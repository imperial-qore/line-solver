function runtime = runAnalyzer(self, options)
% RUNTIME = RUN()
% Run the solver

sn = getStruct(self); % doesn't need initial state

T0=tic;
iter = NaN;
if nargin<2
    options = self.getOptions;
end

self.runAnalyzerChecks(options);
Solver.resetRandomGeneratorSeed(options.seed);
method = options.method;
wasDefault = false;

switch options.method
    case 'default'
        if sn.nstations == 2 && ~any(sn.nodetype == NodeType.ID_CACHE) && any(sn.nodetype == NodeType.ID_DELAY) && any(sn.nservers(isfinite(sn.nservers))>1)
            options.method = 'comomld'; % default for multi-server models
            wasDefault = true;
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

if sn.nclosedjobs == 0 && length(sn.nodetype)==3 && all(sort(sn.nodetype)' == sort([NodeType.Source,NodeType.ID_CACHE,NodeType.Sink])) % is a non-rentrant cache
    % random initialization
    for ind = 1:sn.nnodes
        if sn.nodetype(ind) == NodeType.ID_CACHE
            prob = self.model.nodes{ind}.server.hitClass;
            prob(prob>0) = 0.5;
            self.model.nodes{ind}.setResultHitProb(prob);
            self.model.nodes{ind}.setResultMissProb(1-prob);
        end
    end
    self.model.refreshChains();
    % start iteration
    [QN,UN,RN,TN,CN,XN,lG,pij,runtime] = solver_nc_cache_analyzer(sn, options);
    self.result.Prob.itemProb = pij;
    for ind = 1:sn.nnodes
        if sn.nodetype(ind) == NodeType.ID_CACHE
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
    if any(sn.nodetype == NodeType.ID_CACHE) % if integrated caching-queueing
        [QN,UN,RN,TN,CN,XN,lG,hitprob,missprob,runtime,iter] = solver_nc_cacheqn_analyzer(self, options);
        for ind = 1:sn.nnodes
            if sn.nodetype(ind) == NodeType.ID_CACHE
                self.model.nodes{ind}.setResultHitProb(hitprob(ind,:));
                self.model.nodes{ind}.setResultMissProb(missprob(ind,:));
            end
        end
        self.model.refreshChains();
    else % ordinary queueing network
        if ~isempty(sn.lldscaling) || ~isempty(sn.cdscaling)
            [QN,UN,RN,TN,CN,XN,lG,runtime,iter,method] = solver_ncld_analyzer(sn, options);
        else
            switch options.method
                case 'exact'
                    if ~snHasOpenClasses(sn)
                        % multi-servers have already been transformed before
                        [QN,UN,RN,TN,CN,XN,lG,runtime,iter,method] = solver_ncld_analyzer(sn, options);
                    else%if ~snHasClosedClasses(sn)
                        [QN,UN,RN,TN,CN,XN,lG,runtime,iter,method] = solver_nc_analyzer(sn, options);
                    end
                case {'rd','nrp','nr.probit','nrl','nr.logit','comomld'}
                    [QN,UN,RN,TN,CN,XN,lG,runtime,iter,method] = solver_ncld_analyzer(sn, options);
                otherwise
                    [QN,UN,RN,TN,CN,XN,lG,runtime,iter,method] = solver_nc_analyzer(sn, options);
            end
        end
    end
end
% Compute average arrival rate at steady-state
AN = getAvgArvRFromTput(sn, TN, self.getAvgTputHandles());

switch method
    case 'default'
        self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,['default(',options.method,')'],iter);
    otherwise
        if wasDefault
            self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,['default(',options.method,')'],iter);
        else
            self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,method,iter);
        end
end
self.result.Prob.logNormConstAggr = real(lG);
end