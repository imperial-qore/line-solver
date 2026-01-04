function runtime = runAnalyzer(self, options)
% RUNTIME = RUN()
% Run the solver
if nargin<2
    options = self.getOptions;
end

iter = NaN;
self.runAnalyzerChecks(options);
Solver.resetRandomGeneratorSeed(options.seed);

% Show library attribution if verbose and not yet shown
if options.verbose ~= VerboseLevel.SILENT && ~GlobalConstants.isLibraryAttributionShown()
    sn = self.getStruct();
    libs = SolverNC.getLibrariesUsed(sn, options);
    if ~isempty(libs)
        line_printf('The solver will leverage %s.\n', strjoin(libs, ', '));
        GlobalConstants.setLibraryAttributionShown(true);
    end
end

origmethod = options.method;

%options.lang = 'java';

switch options.lang
    case 'java'
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
                self.model.nodes{ind}.setResultHitProb(JLINE.from_jline_matrix(jmodel.getNodeByIndex(ind-1).getHitRatio()));
                self.model.nodes{ind}.setResultMissProb(JLINE.from_jline_matrix(jmodel.getNodeByIndex(ind-1).getMissRatio()));
            end
        end
        %self.model.refreshChains();
        self.model.refreshStruct(true);
        self.setAvgResults(QN,UN,RN,TN,AN,WN,CN,XN,runtime,options.method,lastiter);
        self.result.Prob.logNormConstAggr = lG;
        return
    case 'matlab'
        sn = getStruct(self); % doesn't need initial state

        switch options.method
            case 'default'
                if sn.nstations == 2 && ~any(sn.nodetype == NodeType.Cache) && any(sn.nodetype == NodeType.Delay) && any(sn.nservers(isfinite(sn.nservers))>1)
                    options.method = 'comomld'; % default for multi-server models
                end
            case 'exact'
                if ~self.model.hasProductFormSolution
                    line_error(mfilename,'The exact method requires the model to have a product-form solution. This model does not have one. You can use Network.hasProductFormSolution() to check before running the solver.');
                elseif isempty(sn.lldscaling)
                    % if exact is requested and does not override a lldscaling assigment
                    Nt = sum(sn.njobs);
                    if isfinite(Nt)
                        % trasform multi-server nodes into lld nodes
                        sn.lldscaling = ones(sn.nstations,Nt);
                        for i=1:sn.nstations
                            if sn.nservers(i) > 1 && isfinite(sn.nservers(i))
                                sn.lldscaling(i,:) = min(1:Nt,sn.nservers(i));
                                sn.nservers(i) = 1;
                            end
                        end
                    end
                end
        end

        if self.enableChecks && ~self.supports(self.model)
            line_error(mfilename,'This model contains features not supported by the solver.');
        end

        Solver.resetRandomGeneratorSeed(options.seed);

        if sn.nclosedjobs == 0 && length(sn.nodetype)==3 && all(sort(sn.nodetype)' == sort([NodeType.Source,NodeType.Cache,NodeType.Sink])) % is a non-rentrant cache
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
            [QN,UN,RN,TN,CN,XN,lG,pij,runtime,actualmethod] = solver_nc_cache_analyzer(sn, options);
            self.result.Prob.itemProb = pij;
            for ind = 1:sn.nnodes
                if sn.nodetype(ind) == NodeType.Cache
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
                [QN,UN,RN,TN,CN,XN,lG,hitprob,missprob,runtime,iter,actualmethod] = solver_nc_cacheqn_analyzer(self, options);
                for ind = 1:sn.nnodes
                    if sn.nodetype(ind) == NodeType.Cache
                        self.model.nodes{ind}.setResultHitProb(hitprob(ind,:));
                        self.model.nodes{ind}.setResultMissProb(missprob(ind,:));
                    end
                end
                self.model.refreshStruct(true);  % Force refresh to get updated actualhitprob/actualmissprob
                sn = self.model.sn;
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
                if ~isempty(sn.lldscaling) || ~isempty(sn.cdscaling)
                    [QN,UN,RN,TN,CN,XN,lG,runtime,iter,actualmethod] = solver_ncld_analyzer(sn, options);
                else
                    switch options.method
                        case 'exact'
                            if ~sn_has_open_classes(sn)
                                % multi-servers have already been transformed before
                                [QN,UN,RN,TN,CN,XN,lG,runtime,iter,actualmethod] = solver_ncld_analyzer(sn, options);
                            else%if ~snHasClosedClasses(sn)
                                [QN,UN,RN,TN,CN,XN,lG,runtime,iter,actualmethod] = solver_nc_analyzer(sn, options);
                            end
                        case {'rd','nrp','nrl','comomld'}
                            [QN,UN,RN,TN,CN,XN,lG,runtime,iter,actualmethod] = solver_ncld_analyzer(sn, options);
                        otherwise
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
end