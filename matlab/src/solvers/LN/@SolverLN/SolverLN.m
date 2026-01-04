classdef SolverLN < EnsembleSolver
    % SolverLN Layered network solver for hierarchical performance models
    %
    % SolverLN implements analysis of layered queueing networks (LQNs) which model
    % hierarchical software systems with clients, application servers, and resource
    % layers. It uses iterative decomposition to solve multi-layer models by
    % analyzing each layer separately and propagating service demands between layers.
    %
    % @brief Layered network solver for hierarchical software performance models
    % 
    % Example:
    % @code
    % solver = SolverLN(layered_model, 'maxIter', 100);
    % solver.runAnalyzer();       % Iterative layer analysis
    % metrics = solver.getEnsembleAvg(); % Layer performance metrics
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties %(Hidden) % registries of quantities to update at every iteration
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
        % Under-relaxation state for convergence improvement
        relax_omega;        % Current relaxation factor
        relax_err_history;  % Error history for adaptive mode
        servt_prev;         % Previous service times for relaxation
        residt_prev;        % Previous residence times for relaxation
        tput_prev;          % Previous throughputs for relaxation
        thinkt_prev;        % Previous think times for relaxation
        callservt_prev;     % Previous call service times for relaxation
        % MOL (Method of Layers) properties for hierarchical iteration
        hostLayerIndices;   % Indices of host (processor) layers in ensemble
        taskLayerIndices;   % Indices of task layers in ensemble
        util_prev_host;     % Previous processor utilizations (for delta computation)
        util_prev_task;     % Previous task utilizations (for delta computation)
        mol_it_host_outer;  % MOL host outer iteration counter
        mol_it_task_inner;  % MOL task inner iteration counter
        % Phase-2 support properties
        hasPhase2;          % Flag: model has phase-2 activities
        servt_ph1;          % Phase-1 service time per activity (nidx x 1)
        servt_ph2;          % Phase-2 service time per activity (nidx x 1)
        util_ph1;           % Phase-1 utilization per entry
        util_ph2;           % Phase-2 utilization per entry
        prOvertake;         % Overtaking probability per entry (nentries x 1)
    end

    properties %(Hidden) % performance metrics and related processes
        util;
        util_ilock; % interlock matrix (ntask x ntask), element (i,j) says how much the utilization of task i is imputed to task j
        tput;
        tputproc;
        servt; % this is the mean service time of an activity, which is the response time at the lower layer (if applicable)
        residt; % this is the residence time at the lower layer (if applicable)
        servtproc; % this is the service time process with mean fitted to the servt value
        servtcdf; % this is the cdf of the service time process
        thinkt;
        thinkproc;
        thinktproc;
        entryproc;
        entrycdfrespt;
        callresidt;
        callservt;
        callservtproc;
        callservtcdf;
        ignore; % elements to be ignored (e.g., components disconnected from a REF node)
    end

    properties %(Access = protected, Hidden) % registries of quantities to update at every iteration
        arvproc_classes_updmap; % [modelidx, actidx, node, class]
        thinkt_classes_updmap; % [modelidx, actidx, node, class]
        actthinkt_classes_updmap; % [modelidx, actidx, node, class] for activity think-times
        servt_classes_updmap; % [modelidx, actidx, node, class]
        call_classes_updmap;  % [modelidx, callidx, node, class]
        route_prob_updmap; % [modelidx, actidxfrom, actidxto, nodefrom, nodeto, classfrom, classto]
        unique_route_prob_updmap; % auxiliary cache of unique route_prob_updmap rows
        solverFactory; % function handle to create layer solvers
    end

    methods
        function self = SolverLN(lqnmodel, solverFactory, varargin)
            % SELF = SOLVERLN(MODEL,SOLVERFACTORY,VARARGIN)
            self@EnsembleSolver(lqnmodel, mfilename);

            if any(cellfun(@(s) strcmpi(s,'java'), varargin))
                self.obj = JLINE.SolverLN(JLINE.from_line_layered_network(lqnmodel));
                self.obj.options.verbose = jline.VerboseLevel.SILENT;
            else
                % Default solver factory: Use JMT for open networks, MVA for closed networks
                defaultSolverFactory = @(m) adaptiveSolverFactory(m, self.options);

                if nargin == 1 %case SolverLN(model)
                    solverFactory = defaultSolverFactory;
                    self.setOptions(SolverLN.defaultOptions);
                elseif nargin>1 && isstruct(solverFactory)
                    options = solverFactory;
                    self.setOptions(options);
                    solverFactory = defaultSolverFactory;
                elseif nargin>2 % case SolverLN(model,'opt1',...)
                    if ischar(solverFactory)
                        inputvar = {solverFactory,varargin{:}}; %#ok<CCAT>
                        solverFactory = defaultSolverFactory;
                    else % case SolverLN(model, solverFactory, 'opt1',...)
                        inputvar = varargin;
                    end
                    self.setOptions(Solver.parseOptions(inputvar, SolverLN.defaultOptions));
                else %case SolverLN(model,solverFactory)
                    self.setOptions(SolverLN.defaultOptions);
                end
                self.lqn = lqnmodel.getStruct();

                % Detect and initialize phase-2 support
                if isfield(self.lqn, 'actphase') && any(self.lqn.actphase > 1)
                    self.hasPhase2 = true;
                    self.servt_ph1 = zeros(self.lqn.nidx, 1);
                    self.servt_ph2 = zeros(self.lqn.nidx, 1);
                    self.util_ph1 = zeros(self.lqn.nidx, 1);
                    self.util_ph2 = zeros(self.lqn.nidx, 1);
                    self.prOvertake = zeros(self.lqn.nentries, 1);
                else
                    self.hasPhase2 = false;
                end

                self.construct();
                for e=1:self.getNumberOfModels
                    if numel(find(self.lqn.isfunction == 1))
                        if ~isempty(self.ensemble{e}.stations{2}.setupTime)
                            solverFactory = @(m) SolverMAM(m,'verbose',false,'method','dec.poisson');
                        else
                            solverFactory = @(m) SolverAuto(m,'verbose',false);
                        end
                        self.setSolver(solverFactory(self.ensemble{e}),e);
                    else
                        self.setSolver(solverFactory(self.ensemble{e}),e);
                    end
                end
                self.solverFactory = solverFactory; % Store for later use
            end
        end

        function runtime = runAnalyzer(self, options) %#ok<INUSD> % generic method to run the solver
            line_error(mfilename,'Use getEnsembleAvg instead.');
        end

        function sn = getStruct(self)
            % SN = GETSTRUCT()

            % Get data structure summarizing the model
            sn = self.model.getStruct();
        end

        function construct(self)
            % mark down to ignore unreachable disconnected components
            self.ignore = false(self.lqn.nidx,1);
            [~,wccs] = weaklyconncomp(self.lqn.graph'+self.lqn.graph);
            uwccs = unique(wccs);
            if length(uwccs)>1
                % the model has disconnected submodels
                wccref = false(1,length(uwccs));
                for t=1:self.lqn.ntasks
                    tidx = self.lqn.tshift+t;
                    if self.lqn.sched(tidx) == SchedStrategy.REF
                        wccref(wccs(tidx)) = true;
                    end
                end
                if any(wccref==false)
                    for dw=find(wccref==false) % disconnected component
                        self.ignore(find(wccs==dw)) = true;
                    end
                end
            end

            % initialize internal data structures
            self.entrycdfrespt = cell(length(self.lqn.nentries),1);
            self.hasconverged = false;

            % initialize svc and think times
            self.servtproc = self.lqn.hostdem;
            self.thinkproc = self.lqn.think;
            self.callservtproc = cell(self.lqn.ncalls,1);
            for cidx = 1:self.lqn.ncalls
                self.callservtproc{cidx} = self.lqn.hostdem{self.lqn.callpair(cidx,2)};
            end

            % perform layering
            self.njobs = zeros(self.lqn.tshift + self.lqn.ntasks, self.lqn.tshift + self.lqn.ntasks);
            buildLayers(self); % build layers
            self.njobsorig = self.njobs;
            self.nlayers = length(self.ensemble);

            % initialize data structures for interlock correction
            self.ptaskcallers = zeros(self.lqn.nhosts+self.lqn.ntasks, self.lqn.nhosts+self.lqn.ntasks);
            self.ptaskcallers_step = cell(1,self.nlayers+1);
            for step=1:self.nlayers % upper bound on maximum dag height
                self.ptaskcallers_step{step} = zeros(self.lqn.nhosts+self.lqn.ntasks, self.lqn.nhosts+self.lqn.ntasks);
            end

            % layering generates update maps that we use here to cache the elements that need reset
            self.routereset = unique(self.idxhash(self.route_prob_updmap(:,1)))';
            self.svcreset = unique(self.idxhash(self.thinkt_classes_updmap(:,1)))';
            self.svcreset = union(self.svcreset,unique(self.idxhash(self.call_classes_updmap(:,1)))');
        end

        function self = reset(self)
            % no-op
        end

        bool = converged(self, it); % convergence test at iteration it

        function init(self) % operations before starting to iterate
            % INIT() % OPERATIONS BEFORE STARTING TO ITERATE
            self.unique_route_prob_updmap = unique(self.route_prob_updmap(:,1))';
            self.tput = zeros(self.lqn.nidx,1);
            self.tputproc = cell(self.lqn.nidx,1);
            self.util = zeros(self.lqn.nidx,1);
            self.servt = zeros(self.lqn.nidx,1);
            self.servtmatrix = getEntryServiceMatrix(self);

            for e= 1:self.nlayers
                self.solvers{e}.enableChecks=false;
            end

            % Initialize under-relaxation state
            relax_mode = self.options.config.relax;
            switch relax_mode
                case {'auto'}
                    self.relax_omega = 1.0; % Start without relaxation
                case {'fixed', 'adaptive'}
                    self.relax_omega = self.options.config.relax_factor;
                otherwise % 'none' or unrecognized
                    self.relax_omega = 1.0; % No relaxation
            end
            self.relax_err_history = [];
            self.servt_prev = NaN(self.lqn.nidx, 1);
            self.residt_prev = NaN(self.lqn.nidx, 1);
            self.tput_prev = NaN(self.lqn.nidx, 1);
            self.thinkt_prev = NaN(self.lqn.nidx, 1);
            self.callservt_prev = NaN(self.lqn.ncalls, 1);

            % Initialize MOL-specific state
            self.mol_it_host_outer = 0;
            self.mol_it_task_inner = 0;
            self.util_prev_host = zeros(self.lqn.nhosts, 1);
            self.util_prev_task = zeros(self.lqn.ntasks, 1);
        end


        function pre(self, it) %#ok<INUSD> % operations before an iteration
            % PRE(IT) % OPERATIONS BEFORE AN ITERATION
            % no-op
        end

        function [result, runtime] = analyze(self, it, e) %#ok<INUSD>
            % [RESULT, RUNTIME] = ANALYZE(IT, E)
            T0 = tic;
            result = struct();
            %jresult = struct();
            if e==1 && self.solvers{e}.options.verbose
                line_printf('\n');
            end

            % Protection for unstable queues during LN iterations
            % If a solver fails (e.g., due to queue instability with open arrivals),
            % use results from previous iteration if available and continue
            try
                [result.QN, result.UN, result.RN, result.TN, result.AN, result.WN] = self.solvers{e}.getAvg();
            catch ME
                if it > 1 && ~isempty(self.results) && size(self.results, 1) >= (it-1) && size(self.results, 2) >= e
                    if self.solvers{e}.options.verbose
                        warning('LINE:SolverLN:Instability', ...
                            'Layer %d at iteration %d encountered instability (possibly due to high service demand with open arrivals). Using previous iteration values and continuing.', ...
                            e, it);
                    end
                    % Use results from previous iteration
                    prevResult = self.results{it-1, e};
                    result.QN = prevResult.QN;
                    result.UN = prevResult.UN;
                    result.RN = prevResult.RN;
                    result.TN = prevResult.TN;
                    result.AN = prevResult.AN;
                    result.WN = prevResult.WN;
                else
                    % First iteration or no previous results, re-throw the exception
                    error('LINE:SolverLN:FirstIterationFailure', ...
                        'Layer %d failed at iteration %d with no previous iteration to fall back on: %s', ...
                        e, it, ME.message);
                end
            end
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
                    case {'SolverMVA', 'SolverNC'} %leaner than refreshProcesses, no need to refresh phases
                        % note: this does not refresh the sn.proc field, only sn.rates and sn.scv
                        switch self.options.method
                            case 'default'
                                self.ensemble{e}.refreshRates();
                            case 'moment3'
                                self.ensemble{e}.refreshProcesses();
                        end
                    otherwise
                        self.ensemble{e}.refreshProcesses();
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
                    self.solvers{e}.setChecks(false);
                end
            end
        end


        function finish(self) % operations after iterations are completed
            % FINISH() % OPERATIONS AFTER INTERATIONS ARE COMPLETED
            E = size(self.results,2);
            for e=1:E
                s = self.solvers{e};
                s.getAvg();
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

        function [cdfRespT] = getCdfRespT(self)
            if isempty(self.entrycdfrespt{1})
                % save user-specified method to temporary variable
                curMethod = self.getOptions.method;
                % run with moment 3
                self.options.method = 'moment3';
                self.getAvg();
                % restore user-specified method
                self.options.method = curMethod;
            end
            cdfRespT = self.entrycdfrespt;
        end

        function [AvgTable,QT,UT,RT,WT,AT,TT] = getAvgTable(self)
            % [AVGTABLE,QT,UT,RT,WT,TT] = GETAVGTABLE(USELQNSNAMING)
            if (GlobalConstants.DummyMode)
                [AvgTable, QT, UT, RT, TT, WT] = deal([]);
                return
            end

            if ~isempty(self.obj)
                avgTable = self.obj.getEnsembleAvg();
                [QN,UN,RN,WN,AN,TN] = JLINE.arrayListToResults(avgTable);
            else
                [QN,UN,RN,TN,AN,WN] = getAvg(self);
            end

            % attempt to sanitize small numerical perturbations
            variables = {QN, UN, RN, TN, AN, WN};  % Put all variables in a cell array
            for i = 1:length(variables)
                rVar = round(variables{i} * 10);
                toRound = abs(variables{i} * 10 - rVar) < GlobalConstants.CoarseTol * variables{i} * 10;
                variables{i}(toRound) = rVar(toRound) / 10;
                variables{i}(variables{i}<=GlobalConstants.FineTol) = 0;
            end
            [QN, UN, RN, TN, AN, WN] = deal(variables{:});  % Assign the modified values back to the original variables

            %%
            Node = label(self.lqn.names);
            O = length(Node);
            NodeType = label(O,1);
            for o = 1:O
                switch self.lqn.type(o)
                    case LayeredNetworkElement.PROCESSOR
                        NodeType(o,1) = label({'Processor'});
                    case LayeredNetworkElement.TASK                        
                        if self.lqn.isref(o)
                            NodeType(o,1) = label({'RefTask'});
                        else
                            NodeType(o,1) = label({'Task'});
                        end
                    case LayeredNetworkElement.ENTRY
                        NodeType(o,1) = label({'Entry'});
                    case LayeredNetworkElement.ACTIVITY
                        NodeType(o,1) = label({'Activity'});
                    case LayeredNetworkElement.CALL
                        NodeType(o,1) = label({'Call'});
                end
            end
            QLen = QN;
            QT = Table(Node,QLen);
            Util = UN;
            UT = Table(Node,Util);
            RespT = RN;
            RT = Table(Node,RespT);
            Tput = TN;
            TT = Table(Node,Tput);
            %SvcT = SN;
            %ST = Table(Node,SvcT);
            %ProcUtil = PN;
            %PT = Table(Node,ProcUtil);
            ResidT = WN;
            WT = Table(Node,ResidT);
            ArvR = AN;
            AT = Table(Node,ArvR);
            AvgTable = Table(Node, NodeType, QLen, Util, RespT, ResidT, ArvR, Tput);%, ProcUtil, SvcT);
        end
    end

    methods
        [QN,UN,RN,TN,AN,WN] = getEnsembleAvg(self);

        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)
            % This method cannot be static as otherwise it cannot access self.solvers{e}
            ensemble = model.getEnsemble;
            featSupported = cell(length(ensemble),1);
            bool = true;
            for e = 1:length(ensemble)
                [solverSupports,featSupported{e}] = self.solvers{e}.supports(ensemble{e});
                bool = bool && solverSupports;
            end
        end
    end

    methods (Hidden)
        buildLayers(self, lqn, resptproc, callservtproc);
        buildLayersRecursive(self, idx, callers, ishostlayer);
        updateLayers(self, it);
        updatePopulations(self, it);
        updateThinkTimes(self, it);
        updateMetrics(self, it);
        updateRoutingProbabilities(self, it);
        svcmatrix = getEntryServiceMatrix(self)
        prOt = overtake_prob(self, eidx);  % Phase-2 overtaking probability

        function delta = computeTaskDelta(self)
            % COMPUTETASKDELTA Compute max queue-length change for task layers (MOL inner loop)
            %
            % Returns the maximum queue-length change across task layers,
            % normalized by total jobs, similar to converged.m logic.
            results = self.results;
            it = size(results, 1);
            if it < 2
                delta = Inf;
                return;
            end

            delta = 0;
            for e = self.taskLayerIndices
                metric = results{it, e}.QN;
                metric_1 = results{it-1, e}.QN;
                N = sum(self.ensemble{e}.getNumberOfJobs);
                if N > 0
                    try
                        iterErr = max(abs(metric(:) - metric_1(:))) / N;
                    catch
                        iterErr = 0;
                    end
                    delta = max(delta, iterErr);
                end
            end
        end

        function delta = computeHostDelta(self)
            % COMPUTEHOSTDELTA Compute max queue-length change for host layers (MOL outer loop)
            %
            % Returns the maximum queue-length change across host layers,
            % normalized by total jobs, similar to converged.m logic.
            results = self.results;
            it = size(results, 1);
            if it < 2
                delta = Inf;
                return;
            end

            delta = 0;
            for e = self.hostLayerIndices
                metric = results{it, e}.QN;
                metric_1 = results{it-1, e}.QN;
                N = sum(self.ensemble{e}.getNumberOfJobs);
                if N > 0
                    try
                        iterErr = max(abs(metric(:) - metric_1(:))) / N;
                    catch
                        iterErr = 0;
                    end
                    delta = max(delta, iterErr);
                end
            end
        end

        function postTaskLayers(self, it)
            % POSTTASKLAYERS Post-iteration updates for task layers only (MOL inner loop)
            %
            % Updates think times, layer parameters, and routing probabilities
            % after task layer analysis, then resets task layer solvers.

            self.updateThinkTimes(it);
            if self.options.config.interlocking
                self.updatePopulations(it);
            end
            self.updateLayers(it);
            self.updateRoutingProbabilities(it);

            % Reset task layer solvers
            for e = self.taskLayerIndices
                switch self.solvers{e}.name
                    case {'SolverMVA', 'SolverNC'}
                        switch self.options.method
                            case 'mol'
                                self.ensemble{e}.refreshRates();
                            otherwise
                                self.ensemble{e}.refreshRates();
                        end
                    otherwise
                        self.ensemble{e}.refreshProcesses();
                end
                self.solvers{e}.reset();
            end

            if self.options.config.interlocking
                for e = self.taskLayerIndices
                    self.ensemble{e}.refreshJobs();
                end
            end
        end

        function postHostLayers(self, it)
            % POSTHOSTLAYERS Post-iteration updates for host layers (MOL outer loop)
            %
            % Updates think times, layer parameters, and routing probabilities
            % after host layer analysis, then resets host layer solvers.

            self.updateThinkTimes(it);
            if self.options.config.interlocking
                self.updatePopulations(it);
            end
            self.updateLayers(it);
            self.updateRoutingProbabilities(it);

            % Reset host layer solvers
            for e = self.hostLayerIndices
                switch self.solvers{e}.name
                    case {'SolverMVA', 'SolverNC'}
                        self.ensemble{e}.refreshRates();
                    otherwise
                        self.ensemble{e}.refreshProcesses();
                end
                self.solvers{e}.reset();
            end

            if self.options.config.interlocking
                for e = self.hostLayerIndices
                    self.ensemble{e}.refreshJobs();
                end
            end
        end
    end

    methods
        function state = get_state(self)
            % GET_STATE Export current solver state for continuation
            %
            % STATE = GET_STATE() returns a struct containing the current
            % solution state, which can be used to continue iteration with
            % a different solver via set_state().
            %
            % The exported state includes:
            % - Service time processes (servtproc)
            % - Think time processes (thinktproc)
            % - Call service time processes (callservtproc)
            % - Throughput processes (tputproc)
            % - Performance metrics (util, tput, servt, residt, etc.)
            % - Relaxation state
            % - Last iteration results
            %
            % Example:
            %   solver1 = SolverLN(model, @(m) SolverMVA(m));
            %   solver1.getEnsembleAvg();
            %   state = solver1.get_state();
            %
            %   solver2 = SolverLN(model, @(m) SolverJMT(m));
            %   solver2.set_state(state);
            %   solver2.getEnsembleAvg();  % Continues from MVA solution

            state = struct();

            % Service/think time processes
            state.servtproc = self.servtproc;
            state.thinktproc = self.thinktproc;
            state.callservtproc = self.callservtproc;
            state.tputproc = self.tputproc;
            state.entryproc = self.entryproc;

            % Performance metrics
            state.util = self.util;
            state.tput = self.tput;
            state.servt = self.servt;
            state.residt = self.residt;
            state.thinkt = self.thinkt;
            state.callresidt = self.callresidt;
            state.callservt = self.callservt;

            % Relaxation state
            state.relax_omega = self.relax_omega;
            state.servt_prev = self.servt_prev;
            state.residt_prev = self.residt_prev;
            state.tput_prev = self.tput_prev;
            state.thinkt_prev = self.thinkt_prev;
            state.callservt_prev = self.callservt_prev;

            % Results from last iteration
            state.results = self.results;

            % Interlock data
            state.njobs = self.njobs;
            state.ptaskcallers = self.ptaskcallers;
            state.ilscaling = self.ilscaling;
        end

        function set_state(self, state)
            % SET_STATE Import solution state for continuation
            %
            % SET_STATE(STATE) initializes the solver with a previously
            % exported state, allowing iteration to continue from where
            % a previous solver left off.
            %
            % This enables hybrid solving schemes where fast solvers (MVA)
            % provide initial estimates and accurate solvers (JMT, DES)
            % refine the solution.
            %
            % Example:
            %   solver1 = SolverLN(model, @(m) SolverMVA(m));
            %   solver1.getEnsembleAvg();
            %   state = solver1.get_state();
            %
            %   solver2 = SolverLN(model, @(m) SolverJMT(m));
            %   solver2.set_state(state);
            %   solver2.getEnsembleAvg();  % Continues from MVA solution

            % Service/think time processes
            self.servtproc = state.servtproc;
            self.thinktproc = state.thinktproc;
            self.callservtproc = state.callservtproc;
            self.tputproc = state.tputproc;
            if isfield(state, 'entryproc')
                self.entryproc = state.entryproc;
            end

            % Performance metrics
            self.util = state.util;
            self.tput = state.tput;
            self.servt = state.servt;
            if isfield(state, 'residt')
                self.residt = state.residt;
            end
            if isfield(state, 'thinkt')
                self.thinkt = state.thinkt;
            end
            if isfield(state, 'callresidt')
                self.callresidt = state.callresidt;
            end
            if isfield(state, 'callservt')
                self.callservt = state.callservt;
            end

            % Relaxation state
            if isfield(state, 'relax_omega')
                self.relax_omega = state.relax_omega;
            end
            if isfield(state, 'servt_prev')
                self.servt_prev = state.servt_prev;
            end
            if isfield(state, 'residt_prev')
                self.residt_prev = state.residt_prev;
            end
            if isfield(state, 'tput_prev')
                self.tput_prev = state.tput_prev;
            end
            if isfield(state, 'thinkt_prev')
                self.thinkt_prev = state.thinkt_prev;
            end
            if isfield(state, 'callservt_prev')
                self.callservt_prev = state.callservt_prev;
            end

            % Results
            if isfield(state, 'results')
                self.results = state.results;
            end

            % Interlock data
            if isfield(state, 'njobs')
                self.njobs = state.njobs;
            end
            if isfield(state, 'ptaskcallers')
                self.ptaskcallers = state.ptaskcallers;
            end
            if isfield(state, 'ilscaling')
                self.ilscaling = state.ilscaling;
            end

            % Update layer models with imported state
            it = 1;
            if ~isempty(self.results)
                it = size(self.results, 1);
            end
            self.updateLayers(it);

            % Refresh all layer solvers with new parameters
            for e = 1:self.nlayers
                self.ensemble{e}.refreshChains();
                switch self.solvers{e}.name
                    case {'SolverMVA', 'SolverNC'}
                        self.ensemble{e}.refreshRates();
                    otherwise
                        self.ensemble{e}.refreshProcesses();
                end
                self.solvers{e}.reset();
            end
        end

        function update_solver(self, solverFactory)
            % UPDATE_SOLVER Change the solver for all layers
            %
            % UPDATE_SOLVER(FACTORY) replaces all layer solvers with
            % new solvers created by the given factory function.
            %
            % This allows switching between different solving methods
            % (e.g., from MVA to JMT/DES) while preserving the current
            % solution state.
            %
            % Example:
            %   solver = SolverLN(model, @(m) SolverMVA(m));
            %   solver.getEnsembleAvg();  % Fast initial solution
            %
            %   solver.update_solver(@(m) SolverJMT(m, 'samples', 1e5));
            %   solver.getEnsembleAvg();  % Refine with simulation

            self.solverFactory = solverFactory;

            % Replace all layer solvers
            for e = 1:self.nlayers
                self.setSolver(solverFactory(self.ensemble{e}), e);
            end
        end

        function [allMethods] = listValidMethods(self)
            sn = self.model.getStruct();
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            allMethods = {'default','moment3','mol'};
        end
    end

    methods (Static)
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = SolverOptions('LN');
        end

        function libs = getLibrariesUsed(sn, options)
            % GETLIBRARIESUSED Get list of external libraries used by LN solver
            % LN uses internal algorithms, no external library attribution needed
            libs = {};
        end
    end
end

function solver = adaptiveSolverFactory(model, parentOptions)
    % ADAPTIVESOLVERFACTORY - Select appropriate solver based on model characteristics
    % Use JMT for models with open classes, MVA for pure closed networks
    if nargin < 2
        verbose = false;
    else
        verbose = parentOptions.verbose;
    end

    solver = SolverMVA(model, 'verbose', verbose);
end
