function runtime = runAnalyzer(self, options)
% RUNTIME = RUN()
% Run the solver
if nargin<2
    options = self.getOptions;
end
% Wall-clock time-budget launch marker (see options.timeout / lineTimeoutExceeded)
options.timeout_tic = tic;

[pyHandled, options, pyRuntime] = self.runAnalyzerPreamble(options, 'NC');
if pyHandled
    runtime = pyRuntime;
    return
end

iter = NaN;
self.runAnalyzerChecks(options);

sn = self.getStruct();
if isfield(sn,'immfeed') && ~isempty(sn.immfeed) && any(sn.immfeed(:))
    line_warning(mfilename,'SolverNC does not handle immediate feedback (immfeed); the solver will treat self-loops as class-switching with re-queueing.\n');
end

Solver.resetRandomGeneratorSeed(options.seed);


origmethod = options.method;

% With enableChecks=false (SolverLN layer backend) an unknown method silently
% falls through to the default analyzer; see _kb/06-solver-catalog.md (NC section)
validMethods = [self.listValidMethods(), {'comomld'}]; % comomld is auto-selected internally from 'default'
if ~any(strcmpi(origmethod, validMethods))
    line_debug(options, 'NC: unrecognized method ''%s'', falling back to the default normalizing-constant analyzer (nc_analyzer/comom).', origmethod);
end

%options.lang = 'java';

switch options.lang
    case 'java'
        line_debug(options, 'NC: using lang=java, delegating to JLINE');
        sn = self.getStruct;
        jmodel = LINE2JLINE(self.model);
        %M = jmodel.getNumberOfStatefulNodes;
        M = jmodel.getNumberOfStations;
        R = jmodel.getNumberOfClasses;
        jsolver = JLINE.SolverNC(jmodel, options);
        [QN,UN,RN,WN,AN,TN] = JLINE.arrayListToResults(jsolver.getAvgTable);
        runtime = jsolver.result.runtime;
        CN = [];
        XN = [];
        QN = reshape(QN',R,M)';
        UN = reshape(UN',R,M)';
        RN = reshape(RN',R,M)';
        TN = reshape(TN',R,M)';
        WN = reshape(WN',R,M)';
        AN = reshape(AN',R,M)';
        lG = NaN;
        lastiter = NaN;
        for ind = 1:sn.nnodes
            if sn.nodetype(ind) == NodeType.Cache
                jnode = jmodel.getNodeByIndex(ind-1);
                self.model.nodes{ind}.setResultHitProb(JLINE.from_jline_matrix(jnode.getHitRatio()));
                self.model.nodes{ind}.setResultMissProb(JLINE.from_jline_matrix(jnode.getMissRatio()));
                % Retrieval-cache extras (delayed-hit ratio, per-list hit ratio
                % and expected latency) so getAvgCacheTable matches the native path.
                self.model.nodes{ind}.setResultDelayedHitProb(JLINE.from_jline_matrix(jnode.getDelayedHitRatio()));
                self.model.nodes{ind}.setResultHitProbList(JLINE.from_jline_matrix(jnode.getHitRatioByList()));
                self.model.nodes{ind}.setResultItemProb(JLINE.from_jline_matrix(jnode.getItemProb()));
                self.model.nodes{ind}.setResultResidT(JLINE.from_jline_matrix(jnode.getResidT()));
            end
        end
        %self.model.refreshChains();
        self.model.refreshStruct(true);
        self.setAvgResults(QN,UN,RN,TN,AN,WN,CN,XN,runtime,options.method,lastiter);
        self.result.Prob.logNormConstAggr = lG;
        return
    case 'matlab'
        line_debug(options, 'NC: using lang=matlab');
        sn = getStruct(self); % doesn't need initial state

        % Fork-join: solver-agnostic fjFixedPoint with ncDispatch as inner solve,
        % intercepted first; see _kb/06-solver-catalog.md (MVA section)
        if self.model.hasFork
            fjres = self.fjFixedPoint(options, @(sn_, opt_) self.ncDispatch(sn_, opt_));
            QN = fjres.QN; UN = fjres.UN; RN = fjres.RN; TN = fjres.TN;
            CN = fjres.CN; XN = fjres.XN; lG = fjres.lG;
            runtime = fjres.runtime; iter = fjres.iter;
            actualmethod = fjres.actualmethod;
            sn = self.model.getStruct();
            AN = sn_get_arvr_from_tput(sn, TN, self.getAvgTputHandles());
            if strcmp(origmethod,'default') && ~isempty(actualmethod) && ~strcmp(actualmethod,'default')
                self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,['default/' actualmethod],iter);
            else
                self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,options.method,iter);
            end
            self.result.Prob.logNormConstAggr = real(lG);
            if lineTimeoutExceeded(options)
                self.result.Avg.timedOut = true;
                line_warning(mfilename,'Solver exceeded the wall-clock time budget (options.timeout=%gs).\n', options.timeout);
            end
            return
        end

        % OI closed network -> solver_nc_oi_analyzer (exact pfqn_ncoi), intercepted
        % before the multiserver->lldscaling conversion; see
        % _kb/06-solver-catalog.md (NC section, analyzer routing)
        if nc_is_oi_model(sn) && any(strcmpi(options.method,{'default','exact'}))
            line_debug(options, 'NC: order-independent closed network, routing to solver_nc_oi_analyzer');
            [QN,UN,RN,TN,CN,XN,lG,runtime,iter,actualmethod] = solver_nc_oi_analyzer(sn, options);
            AN = sn_get_arvr_from_tput(sn, TN, self.getAvgTputHandles());
            if strcmp(origmethod,'default') && ~strcmp(actualmethod,'default')
                self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,['default/' actualmethod],iter);
            else
                self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,actualmethod,iter);
            end
            self.result.Prob.logNormConstAggr = real(lG);
            return
        end

        switch options.method
            case 'default'
                if sn.nstations == 2 && ~any(sn.nodetype == NodeType.Cache) && any(sn.nodetype == NodeType.Delay) && any(sn.nservers(isfinite(sn.nservers))>1)
                    % 2-station Delay+multiserver (every SolverLN layer submodel):
                    % exact load-dependent CoMoM when product-form, else Seidmann
                    % (comom); see _kb/06-solver-catalog.md (NC section)
                    if self.model.hasProductFormSolution && isempty(sn.lldscaling)
                        Nt = sum(sn.njobs);
                        if isfinite(Nt)
                            sn.lldscaling = ones(sn.nstations,Nt);
                            for i=1:sn.nstations
                                if sn.nservers(i) > 1 && isfinite(sn.nservers(i))
                                    % keep server count c for U normalization; see
                                    % _kb/06-solver-catalog.md (NC section)
                                    sn.lldscaling(i,:) = min(1:Nt,sn.nservers(i));
                                end
                            end
                        end
                        line_debug(options, 'Default method: 2-station Delay+multiserver product-form network, using exact load-dependent comomld');
                    else
                        options.method = 'comom'; % non-product-form (e.g. LN layer submodels): Seidmann approximation
                        line_debug(options, 'Default method: 2-station Delay+multiserver non-product-form network, using comom');
                    end
                end
            case {'exact','is','panaceald'}
                % 'is' and 'panaceald' need the same model as 'exact' (same
                % multiserver conversion),
                % except OI/P&S which route to solver_nc_analyzer (pfqn_pas_is /
                % pfqn_oi_is); see _kb/06-solver-catalog.md (NC section)
                if strcmpi(options.method,'is') && nc_is_pas_model(sn)
                    % no-op: handled by solver_nc_analyzer (pfqn_pas_is)
                elseif ~self.model.hasProductFormSolution
                    line_error(mfilename,'The %s method requires the model to have a product-form solution. This model does not have one. You can use Network.hasProductFormSolution() to check before running the solver.', options.method);
                elseif isempty(sn.lldscaling) && any(sn.nservers(isfinite(sn.nservers)) > 1)
                    % convert multiserver to lld ONLY when a genuine multiserver is
                    % present (forcing it on an all-single-server model mishandles a
                    % self-looping closed chain); see _kb/06-solver-catalog.md (NC section)
                    Nt = sum(sn.njobs);
                    if isfinite(Nt)
                        sn.lldscaling = ones(sn.nstations,Nt);
                        for i=1:sn.nstations
                            if sn.nservers(i) > 1 && isfinite(sn.nservers(i))
                                sn.lldscaling(i,:) = min(1:Nt,sn.nservers(i));
                            end
                        end
                    end
                    line_debug(options, 'Exact method: converted multiserver stations to load-dependent (Nt=%d)', Nt);
                end
        end

        % The feature gate ran in runAnalyzerChecks above, through
        % supportsModelMethod, and reports which features are unsupported. A
        % coarse supports(self.model) repeat here would only re-reject the same
        % models with a message that names nothing.

        Solver.resetRandomGeneratorSeed(options.seed);

        ci_cache = find(sn.nodetype == NodeType.Cache, 1);
        hasRetrieval = ~isempty(ci_cache) && isfield(sn.nodeparam{ci_cache}, 'retrievalSystemCapacity') ...
            && sn.nodeparam{ci_cache}.retrievalSystemCapacity > 0;
        if hasRetrieval % delayed-hit (retrieval-system) cache
            if any(sn.nodetype == NodeType.Source)
                line_debug(options, 'Open delayed-hit retrieval cache, routing to nc_retrieval_analyzer');
                [QN,UN,RN,TN,CN,XN,lG,hitprob,missprob,delayedprob,hitproblist,itemprob,latency,runtime,actualmethod] = solver_nc_retrieval_analyzer(sn, options);
            else
                line_debug(options, 'Closed integrated delayed-hit retrieval cache, routing to nc_cacheqn_retrieval_analyzer');
                [QN,UN,RN,TN,CN,XN,lG,hitprob,missprob,delayedprob,hitproblist,itemprob,latency,runtime,actualmethod] = solver_nc_cacheqn_retrieval_analyzer(sn, options);
            end
            iter = NaN;
            for ind = 1:sn.nnodes
                if sn.nodetype(ind) == NodeType.Cache
                    self.model.nodes{ind}.setResultHitProb(hitprob);
                    self.model.nodes{ind}.setResultMissProb(missprob);
                    self.model.nodes{ind}.setResultDelayedHitProb(delayedprob);
                    self.model.nodes{ind}.setResultHitProbList(hitproblist);
                    self.model.nodes{ind}.setResultItemProb(itemprob);
                    self.model.nodes{ind}.setResultResidT(latency);
                end
            end
            self.model.refreshStruct(true);
        elseif sn.nclosedjobs == 0 && length(sn.nodetype)==3 && all(sort(sn.nodetype)' == sort([NodeType.Source,NodeType.Cache,NodeType.Sink])) % is a non-rentrant cache
            line_debug(options, 'Non-reentrant cache (Source-Cache-Sink), routing to nc_cache_analyzer');
            % random initialization
            for ind = 1:sn.nnodes
                if sn.nodetype(ind) == NodeType.Cache
                    prob = self.model.nodes{ind}.server.hitClass;
                    prob(prob>0) = 0.5;
                    self.model.nodes{ind}.setResultHitProb(prob);
                    self.model.nodes{ind}.setResultMissProb(1-prob);
                end
            end
            self.model.refreshChains();
            % start iteration
            [QN,UN,RN,TN,CN,XN,lG,pij,runtime,actualmethod,hitproblist,itemprob] = solver_nc_cache_analyzer(sn, options);
            self.result.Prob.itemProb = pij;
            for ind = 1:sn.nnodes
                if sn.nodetype(ind) == NodeType.Cache
                    self.model.nodes{ind}.setResultHitProbList(hitproblist);
                    self.model.nodes{ind}.setResultItemProb(itemprob);
                    %prob = self.model.nodes{ind}.server.hitClass;
                    %prob(prob>0) = 0.5;
                    hitClass = self.model.nodes{ind}.getHitClass;
                    missClass = self.model.nodes{ind}.getMissClass;
                    hitprob = zeros(1,length(hitClass));
                    for k=1:length(self.model.nodes{ind}.getHitClass)
                        %                for k=1:length(self.model.nodes{ind}.server.hitClass)
                        chain_k = sn.chains(:,k)>0;
                        inchain = sn.chains(chain_k,:)>0;
                        h = hitClass(k);
                        m = missClass(k);
                        if h>0 && m>0
                            hitprob(k) = XN(h) / sum(XN(inchain),"omitnan");
                        end
                    end
                    self.model.nodes{ind}.setResultHitProb(hitprob);
                    self.model.nodes{ind}.setResultMissProb(1-hitprob);
                end
            end
            self.model.refreshChains;
        else % queueing network
            if any(sn.nodetype == NodeType.Cache) % if integrated caching-queueing
                line_debug(options, 'Integrated caching-queueing network, routing to nc_cacheqn_analyzer');
                [QN,UN,RN,TN,CN,XN,lG,hitprob,missprob,runtime,iter,actualmethod] = solver_nc_cacheqn_analyzer(self, options);
                for ind = 1:sn.nnodes
                    if sn.nodetype(ind) == NodeType.Cache
                        self.model.nodes{ind}.setResultHitProb(hitprob(ind,:));
                        self.model.nodes{ind}.setResultMissProb(missprob(ind,:));
                    end
                end
                self.model.refreshStruct(true);  % Force refresh to get updated actualhitprob/actualmissprob
                sn = self.model.sn;
            elseif sn_is_mm1k_loss(sn) % single-station M/M/1/K with tail drop
                % Exact probability-based loss analysis (M/M/1/K stationary
                % distribution); see qsys_mm1k_loss.
                queue_ist = sn.nodeToStation(sn.nodetype == NodeType.Queue);
                source_ist = sn.nodeToStation(sn.nodetype == NodeType.Source);
                Kcap = sn.cap(queue_ist);
                lambda = sn.rates(source_ist)*sn.visits{source_ist}(sn.stationToStateful(queue_ist));
                mu = sn.rates(queue_ist);
                rho = lambda/mu;
                [Ploss, ~] = qsys_mm1k_loss(lambda, mu, Kcap);
                Tq = lambda*(1-Ploss);          % carried throughput
                if abs(rho-1) < 1e-10
                    Lsys = Kcap/2;              % L'Hopital limit at rho=1
                else
                    Lsys = rho/(1-rho) - (Kcap+1)*rho^(Kcap+1)/(1-rho^(Kcap+1));
                end
                Vq = sn.visits{1}(sn.stationToStateful(queue_ist));
                M = sn.nstations;
                QN = zeros(M,1); UN = QN; RN = QN; TN = QN; XN = QN; CN = 0;
                RN(queue_ist) = Lsys/Tq;        % per-visit response time (Little)
                QN(queue_ist) = Lsys;
                UN(queue_ist) = Tq/mu;          % single-server utilization
                TN(queue_ist) = Tq;             % carried (effective) rate
                TN(source_ist) = lambda;        % offered arrival rate
                XN(queue_ist) = Tq;             % system throughput = carried rate
                CN = RN(queue_ist)*Vq;
                lG = 0; iter = 1; actualmethod = 'mm1k.loss';
                AN = sn_get_arvr_from_tput(sn, TN, self.getAvgTputHandles());
                if strcmp(origmethod,'default')
                    self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,0,['default/' actualmethod],iter);
                else
                    self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,0,actualmethod,iter);
                end
                self.result.Prob.logNormConstAggr = real(lG);
                return
            else % ordinary queueing network
                % Check for open model with single FCR containing single Delay (loss network)
                if ~sn_has_closed_classes(sn) && sn.nregions == 1
                    regionMatrix = sn.region{1};
                    % Find stations in FCR (those with non-negative constraints)
                    stationsInFCR = find(any(regionMatrix(:,1:end-1) >= 0, 2) | regionMatrix(:,end) >= 0);
                    if length(stationsInFCR) == 1 && isinf(sn.nservers(stationsInFCR(1)))
                        % Single delay node in FCR - check drop rule
                        if all(sn.regionrule(1,:) == DropStrategy.DROP)
                            % Use loss network solver
                            line_debug(options, 'Open model with single FCR + Delay node (DROP), routing to nc_lossn_analyzer');
                            [QN,UN,RN,TN,CN,XN,lG,runtime,iter,actualmethod] = solver_nc_lossn_analyzer(sn, options);
                            AN = sn_get_arvr_from_tput(sn, TN, self.getAvgTputHandles());
                            if strcmp(origmethod,'default') && exist('actualmethod','var') && ~strcmp(actualmethod,'default')
                                self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,['default/' actualmethod],iter);
                            else
                                self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,actualmethod,iter);
                            end
                            self.result.Prob.logNormConstAggr = real(lG);
                            return
                        else
                            % WAITQ (blocking) not supported - error and stop
                            line_error(mfilename, 'SolverNC does not support finite capacity regions with WAITQ (blocking) policy. Use DROP policy instead.');
                        end
                    end
                end
                % Residual FCR (not the single-Delay loss-network case dispatched
                % above): NC does not enforce the aggregate region limit; reject
                % instead of silently returning the unconstrained answer.
                if isfield(sn,'nregions') && sn.nregions > 0
                    line_error(mfilename,'This model uses a Finite Capacity Region (addRegion) on queueing stations, which is not supported by SolverNC. Use SolverJMT, or setCapacity for a single-station limit.');
                end
                if ~isempty(sn.lldscaling) || ~isempty(sn.cdscaling) || ~isempty(sn.jdscaling)
                    line_debug(options, 'Load-dependent scaling detected, routing to ncld_analyzer');
                    [QN,UN,RN,TN,CN,XN,lG,runtime,iter,actualmethod] = solver_ncld_analyzer(sn, options);
                else
                    switch options.method
                        case 'exact'
                            if ~sn_has_open_classes(sn)
                                % empty lldscaling here means an all-single-server closed
                                % model -> exact normalizing-constant path (ncld would
                                % mishandle a self-looping chain); see
                                % _kb/06-solver-catalog.md (NC section)
                                line_debug(options, 'NC method=exact, single-server closed model, routing to nc_analyzer');
                                [QN,UN,RN,TN,CN,XN,lG,runtime,iter,actualmethod] = solver_nc_analyzer(sn, options);
                            else%if ~snHasClosedClasses(sn)
                                line_debug(options, 'NC method=exact, open model, routing to nc_analyzer');
                                [QN,UN,RN,TN,CN,XN,lG,runtime,iter,actualmethod] = solver_nc_analyzer(sn, options);
                            end
                        case {'rd','nrp','nrl','comomld','panaceald'}
                            line_debug(options, 'NC method=%s, routing to ncld_analyzer', options.method);
                            [QN,UN,RN,TN,CN,XN,lG,runtime,iter,actualmethod] = solver_ncld_analyzer(sn, options);
                        otherwise
                            line_debug(options, 'NC method=%s, routing to nc_analyzer', options.method);
                            [QN,UN,RN,TN,CN,XN,lG,runtime,iter,actualmethod] = solver_nc_analyzer(sn, options);
                    end
                end
            end
        end
        % Compute average arrival rate at steady-state
        AN = sn_get_arvr_from_tput(sn, TN, self.getAvgTputHandles());
end
if strcmp(origmethod,'default') && exist('actualmethod','var') && ~strcmp(actualmethod,'default')
    self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,['default/' actualmethod],iter);
else
    self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,actualmethod,iter);
end

self.result.Prob.logNormConstAggr = real(lG);
if lineTimeoutExceeded(options)
    self.result.Avg.timedOut = true;
    line_warning(mfilename,'Solver exceeded the wall-clock time budget (options.timeout=%gs).\n', options.timeout);
end
end