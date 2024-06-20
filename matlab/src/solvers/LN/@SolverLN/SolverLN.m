classdef SolverLN < LayeredNetworkSolver & EnsembleSolver
    % LINE native solver for layered networks.
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.

    properties (Hidden) % registries of quantities to update at every iteration
        nlayers; % number of model layers
        lqn; % lqn data structure
        hasconverged; % true if last iteration converged, false otherwise
        averagingstart; % iteration at which result averaging started
        idxhash; % ensemble model associated to host or task
        servtmatrix; % auxiliary matrix to determine entry servt
        ptaskcallers; % probability that a task is called by a given task, directly or indirectly (remotely)
        ptaskcallers_step; % probability that a task is called by a given task, directly or indirectly (remotely) up to a given step distance
        ilscaling; % interlock scalings
        njobs; % number of jobs for each caller in a given submodel
        njobsorig; % number of jobs for each caller at layer build time
        routereset; % models that require hard reset of service chains
        svcreset; % models that require hard reset of service process
        maxitererr; % maximum error at current iteration over all layers
        unique_route_prob_updmap; % auxiliary cache of unique route_prob_updmap rows
    end

    properties (Hidden) % performance metrics and related processes
        util;
        util_ilock; % interlock matrix (ntask x ntask), element (i,j) says how much the utilization of task i is imputed to task j
        tput;
        tputproc;
        servt; % this is the mean service time of an activity, which is the residence time at the lower layer (if applicable)
        servtproc; % this is the service time process with mean fitted to the servt value
        servtcdf; % this is the cdf of the service time process
        thinkt;
        thinkproc;
        thinktproc;
        entryproc;
        entrycdfrespt;
        callresidt;
        callresidtproc;
        callresidtcdf;
    end

    properties (Access = protected, Hidden) % registries of quantities to update at every iteration
        arvproc_classes_updmap; % [modelidx, actidx, node, class]
        thinkt_classes_updmap; % [modelidx, actidx, node, class]
        servt_classes_updmap; % [modelidx, actidx, node, class]
        call_classes_updmap;  % [modelidx, callidx, node, class]
        route_prob_updmap; % [modelidx, actidxfrom, actidxto, nodefrom, nodeto, classfrom, classto]
    end

    methods
        function self = SolverLN(lqnmodel, solverFactory, varargin)
            % SELF = SOLVERLN(MODEL,SOLVERFACTORY,VARARGIN)

            self@LayeredNetworkSolver(lqnmodel, mfilename);
            self@EnsembleSolver(lqnmodel, mfilename);

            if any(cellfun(@(s) strcmpi(s,'java'), varargin))
                self.obj = JLINE.SolverLN(JLINE.from_line_layered_network(lqnmodel));
                self.obj.options.verbose = jline.lang.constant.VerboseLevel.SILENT;
            else
                if nargin == 1 %case SolverLN(model)
                    solverFactory = @(m) SolverAuto(m,'verbose',false);
                    self.setOptions(SolverLN.defaultOptions);
                elseif nargin>1 && isstruct(solverFactory)
                    options = solverFactory;
                    self.setOptions(options);
                    solverFactory = @(m) SolverAuto(m,'verbose',false);
                elseif nargin>2 % case SolverLN(model,'opt1',...)
                    if ischar(solverFactory)
                        inputvar = {solverFactory,varargin{:}}; %#ok<CCAT>
                        solverFactory = @(m) SolverAuto(m,'verbose',false);
                    else % case SolverLN(model, solverFactory, 'opt1',...)
                        inputvar = varargin;
                    end
                    self.setOptions(Solver.parseOptions(inputvar, SolverLN.defaultOptions));
                else %case SolverLN(model,solverFactory)
                    self.setOptions(SolverLN.defaultOptions);
                end
                self.lqn = lqnmodel.getStruct();
                self.construct();
                for e=1:self.getNumberOfModels
                    if numel(find(self.lqn.isfunction == 1))
                        if ~isempty(self.ensemble{e}.stations{2}.setupTime)
                            solverFactory = @(m) SolverMAM(m,'verbose',false,'method','dec.poisson');
                        else
                            solverFactory = @(m) SolverMVA(m,'verbose',false);
                        end
                        self.setSolver(solverFactory(self.ensemble{e}),e);
                    else
                        self.setSolver(solverFactory(self.ensemble{e}),e);
                    end
                end
            end
        end

        function construct(self)
            % initialize internal data structures
            self.entrycdfrespt = cell(length(self.lqn.nentries),1);
            self.hasconverged = false;

            % initialize svc and think times
            self.servtproc = self.lqn.hostdem;
            self.thinkproc = self.lqn.think;
            self.callresidtproc = cell(self.lqn.ncalls,1);
            for cidx = 1:self.lqn.ncalls
                self.callresidtproc{cidx} = self.lqn.hostdem{self.lqn.callpair(cidx,2)};
            end

            % perform layering
            self.njobs = zeros(self.lqn.tshift + self.lqn.ntasks, self.lqn.tshift + self.lqn.ntasks);
            buildLayers(self);
            self.njobsorig = self.njobs;
            self.nlayers = length(self.ensemble);

            % initialize data structures for interlock correction
            self.ptaskcallers = zeros(self.lqn.nhosts+self.lqn.ntasks, self.lqn.nhosts+self.lqn.ntasks);
            self.ptaskcallers_step = cell(1,self.nlayers);
            for e=1:self.nlayers
                self.ptaskcallers_step{e} = zeros(self.lqn.nhosts+self.lqn.ntasks, self.lqn.nhosts+self.lqn.ntasks);
            end

            % layering generates update maps that we use here to cache the elements that need reset
            self.routereset = unique(self.idxhash(self.route_prob_updmap(:,1)))';
            self.svcreset = unique(self.idxhash(self.thinkt_classes_updmap(:,1)))';
            self.svcreset = union(self.svcreset,unique(self.idxhash(self.call_classes_updmap(:,1)))');
        end


        function initFromRawAvgTables(self, NodeAvgTable, CallAvgTable)
            line_error(mfilename,'initFromRawAvgTables not yet available');
        end

        function paramFromRawAvgTables(self, NodeAvgTable, CallAvgTable)
            line_error(mfilename,'paramFromRawAvgTables not yet available');
        end

        function self = reset(self)
            % no-op
        end

        bool = converged(self, it); % convergence test at iteration it

        function init(self) % operations before starting to iterate
            % INIT() % OPERATIONS BEFORE STARTING TO ITERATE
            self.unique_route_prob_updmap = unique(self.route_prob_updmap(:,1))';
            self.tput = zeros(self.lqn.nidx,1);
            self.util = zeros(self.lqn.nidx,1);
            self.servt = zeros(self.lqn.nidx,1);
            self.servtmatrix = getEntryServiceMatrix(self);

            for e= 1:self.nlayers
                self.solvers{e}.enableChecks=false;
            end
        end


        function pre(self, it) % operations before an iteration
            % PRE(IT) % OPERATIONS BEFORE AN ITERATION
            % no-op
        end

        function [result, runtime] = analyze(self, it, e)
            % [RESULT, RUNTIME] = ANALYZE(IT, E)
            T0 = tic;
            result = struct();
            jresult = struct();
            %it
            %if it>1%2*length(self.model.ensemble)
            [result.QN, result.UN, result.RN, result.TN, result.AN, result.WN] = self.solvers{e}.getAvg();
            %else
            %    [result.QN, result.UN, result.RN, result.TN, result.AN, result.WN] = SolverMVA(self.model.ensemble{e},' ','verbose',self.solvers{e}.options.verbose).getAvg();
            %end
            runtime = toc(T0);
        end

        function post(self, it) % operations after an iteration
            % POST(IT) % OPERATIONS AFTER AN ITERATION
            % convert the results of QNs into layer metrics
            self.updateMetrics(it);

            % recompute think times
            self.updateThinkTimes(it);

            if self.options.config.interlocking
                % recompute layer populations
                self.updatePopulations(it);
            end

            % update the model parameters
            self.updateLayers(it);

            % update entry selection and cache routing probabilities within callers
            self.updateRoutingProbabilities(it);

            % reset all layers with routing probability changes
            for e= self.routereset
                self.ensemble{e}.refreshChains();
                self.solvers{e}.reset();
            end

            % refresh visits and network model parameters
            for e= self.svcreset
                switch self.solvers{e}.name
                    case {'SolverMVA', 'SolverNC'} %leaner than refreshService, no need to refresh phases
                        % note: this does not refresh the sn.proc field, only sn.rates and sn.scv
                        switch self.options.method
                            case 'default'
                                refreshRates(self.ensemble{e});
                            case 'moment3'
                                refreshService(self.ensemble{e});
                        end
                    otherwise
                        refreshService(self.ensemble{e});
                end
                self.solvers{e}.reset(); % commenting this out des not seem to produce a problem, but it goes faster with it
            end

            % this is required to handle population changes due to interlocking
            if self.options.config.interlocking
                for e=1:self.nlayers
                    self.ensemble{e}.refreshJobs();
                end
            end

            if it==1
                % now disable all solver support checks for future iterations
                for e=1:length(self.ensemble)
                    self.solvers{e}.setDoChecks(false);
                end
            end
        end


        function finish(self) % operations after iterations are completed
            % FINISH() % OPERATIONS AFTER INTERATIONS ARE COMPLETED
            if self.options.verbose
                line_printf('\n');
            end
            E = size(self.results,2);
            for e=1:E
                s = self.solvers{e};
                s.getAvgTable();
                self.solvers{e} = s;
            end
            self.model.ensemble = self.ensemble;
        end

        function [QNlqn_t, UNlqn_t, TNlqn_t] = getTranAvg(self)
            self.getAvg;
            QNclass_t = {};
            UNclass_t = {};
            TNclass_t = {};
            QNlqn_t = cell(0,0);
            for e=1:self.nlayers
                [crows, ccols] = size(QNlqn_t);
                [QNclass_t{e}, UNclass_t{e}, TNclass_t{e}] = self.solvers{e}.getTranAvg();
                QNlqn_t(crows+1:crows+size(QNclass_t{e},1),ccols+1:ccols+size(QNclass_t{e},2)) = QNclass_t{e};
                UNlqn_t(crows+1:crows+size(UNclass_t{e},1),ccols+1:ccols+size(UNclass_t{e},2)) = UNclass_t{e};
                TNlqn_t(crows+1:crows+size(TNclass_t{e},1),ccols+1:ccols+size(TNclass_t{e},2)) = TNclass_t{e};
            end
        end

        function varargout = getAvg(varargin)
            % [QN,UN,RN,TN,AN,WN] = GETAVG(SELF,~,~,~,~,USELQNSNAMING)
            [varargout{1:nargout}] = getEnsembleAvg( varargin{:} );
        end

        %        function [NodeAvgTable, CallAvgTable] = getRawAvgTables(self)
        %            line_error(mfilename,'getRawAvgTables not yet available');
        %        end

        function [cdfRespT] = getCdfRespT(self)
            if isempty(self.entrycdfrespt{1})
                % save user-specified method to temporary variable
                curMethod = self.getOptions.method;
                % run with moment 3
                self.options.method = 'moment3';
                self.getAvgTable();
                % restore user-specified method
                self.options.method = curMethod;
            end
            cdfRespT = self.entrycdfrespt;
        end

    end

    methods
        [QN,UN,RN,TN,AN,WN] = getEnsembleAvg(self,~,~,~,~, useLQNSnaming);
    end

    methods (Hidden)
        buildLayers(self, lqn, resptproc, callresidtproc);
        buildLayersRecursive(self, idx, callers, ishostlayer);
        updateLayers(self, it);
        updatePopulations(self, it);
        updateThinkTimes(self, it);
        updateMetrics(self, it);
        updateRoutingProbabilities(self, it);
        svcmatrix = getEntryServiceMatrix(self)
    end

    methods (Static)
        function [allMethods] = listValidMethods()
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            allMethods = {'default','moment3'};
        end

        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)
            % Todo: if this is static it cannot access self.solvers{e}
            ensemble = model.getEnsemble;
            featSupported = cell(length(ensemble),1);
            bool = true;
            % for e = 1:length(ensemble)
            %     [solverSupports,featSupported{e}] = self.solvers{e}.supports(ensemble{e});
            %     bool = bool && solverSupports;
            % end
        end
    end

    methods (Static)
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = lineDefaults('LN');
        end
    end
end
