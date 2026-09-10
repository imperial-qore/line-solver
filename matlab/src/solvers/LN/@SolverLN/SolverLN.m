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
        momentPassDone; % method='moment3': true once the moment-based entry-law pass has run
        averagingstart; % iteration at which result averaging started
        idxhash; % ensemble model associated to host or task
        servtmatrix; % auxiliary matrix to determine entry servt
        ilscaling; % interlock scalings
        % LQNS V5-style interlock data structures (built once at init)
        il_table_all;    % (nentries x nentries) reachability probability, all phases
        il_table_ph1;    % (nentries x nentries) reachability probability, phase-1 only
        il_common_entries;  % cell(nhosts+ntasks,1) common parent entry indices per server
        il_source_tasks_all; % cell(nhosts+ntasks,1) all-phase source tasks per server
        il_source_tasks_ph2; % cell(nhosts+ntasks,1) phase-2 source tasks per server
        il_num_sources;  % (nhosts+ntasks,1) total source multiplicity per server
        njobs; % number of jobs for each caller in a given submodel
        njobsorig; % number of jobs for each caller at layer build time
        routereset; % models that require hard reset of service chains
        svcreset; % models that require hard reset of service process
        maxitererr; % maximum error at current iteration over all layers
        % Under-relaxation state for convergence improvement
        relax_omega;        % Current relaxation factor
        relax_err_history;  % Error history for adaptive mode
        % Stochastic iteration (Robbins-Monro / Polyak-Ruppert) state,
        % used when one or more layer solvers return noisy estimates
        stochiter_mode;     % resolved mode: 'rm' | 'crn' | 'off'
        stochiter_auto;     % true if mode was resolved from 'auto'
        stochiter_start;    % iteration at which RM averaging started
        stochlayers;        % logical(1,nlayers): layer solver is stochastic
        stoch_avg;          % cell(1,nlayers) Polyak-Ruppert averages of layer results
        stoch_avg_count;    % iterations accumulated into stoch_avg
        stoch_servt_avg;    % Polyak-Ruppert average of the servt iterate
        stoch_residt_avg;   % Polyak-Ruppert average of the residt iterate
        servt_prev;         % Previous service times for relaxation
        residt_prev;        % Previous residence times for relaxation
        tput_prev;          % Previous throughputs for relaxation
        thinkt_prev;        % Previous think times for relaxation
        callservt_prev;     % Previous call service times for relaxation
        callresidt_prev;    % Previous call residence times for growth rate capping
        singleReplicaTasks; % Task indices modeled as single representative replica (fan-out)
        % MOL (Method of Layers) properties for hierarchical iteration
        hostLayerIndices;   % Indices of host (processor) layers in ensemble
        taskLayerIndices;   % Indices of task layers in ensemble
        util_prev_host;     % Previous processor utilizations (for delta computation)
        util_prev_task;     % Previous task utilizations (for delta computation)
        pyMode;             % Flag: delegate to native python SolverLN (lang='python')
        cppMode;            % Flag: delegate to the C++ line-cli layered solver (lang='cpp')
        % Phase-2 support properties
        hasPhase2;          % Flag: model has phase-2 activities
        servt_ph1;          % Phase-1 service time per activity (nidx x 1)
        servt_ph2;          % Phase-2 service time per activity (nidx x 1)
        util_ph1;           % Phase-1 utilization per entry
        util_ph2;           % Phase-2 utilization per entry
        prOvertake;         % Overtaking probability per entry (nentries x 1)
        % the PH encodings ('srvn.ph', 'flat.ph'): per-entry Workflow objects,
        % their composed laws and the layer maps that split a per-caller layer
        % result back over the entries, activities and calls -- see buildLayersPH
        ph;
        % The method the layers were actually BUILT for: 'srvn.ph', 'srvn.cs',
        % 'flat.cs', 'flat.ph' or 'moment3'. Resolved once in buildLayers,
        % because the alias 'srvn' may fall back; every dispatch reads this and
        % not options.method, so a reconstruction can never disagree with the
        % layers it is reading.
        lnmethod;
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
        joint; % join times at AND-Join activities (synchronization delay)
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
        layerHasRegion; % logical per ensemble index, true if the layer carries an admission constraint
        layerChains; % cell per ensemble index, cached sn.chains of a constrained layer
    end

    methods
        function self = SolverLN(lqnmodel, solverFactory, varargin)
            % SELF = SOLVERLN(MODEL,SOLVERFACTORY,VARARGIN)
            self@EnsembleSolver(lqnmodel, mfilename);

            % Solver console: constructing the ensemble compiles one auxiliary
            % network per layer. Each keeps its headline; the stage-by-stage
            % detail of a layer would bury the layered narration itself.
            consoleQuiet = LineConsole.pushQuiet(); %#ok<NASGU>

            % Collect all trailing args (solverFactory may itself be the lang
            % string or an options struct) to detect lang='python'/'java'.
            allArgs = varargin;
            if nargin > 1
                allArgs = [{solverFactory}, varargin];
            end
            % THE LANGUAGE IS NOT ONLY AN ARGUMENT. options.lang is what every
            % other solver dispatches on -- runAnalyzerPreamble for the eight
            % Network solvers, an explicit switch in AUTO and JMT, and
            % self.options.lang at solve time in ENV -- and SolverOptions fills
            % it from the LINEDefaultLang global, itself seeded from
            % LINE_DEFAULT_LANG. A caller that selects the engine GLOBALLY, which
            % is how the cross-language parity rows do it, therefore puts NOTHING
            % in allArgs, and reading the arguments alone left this solver in
            % MATLAB mode on exactly those runs.
            %
            % What that cost was not an error but a silent downgrade: the
            % ensemble delegation below was skipped and the MATLAB fixed point
            % ran instead, dispatching each layer to the foreign engine on its
            % own. Under lang='cpp' that is one line-cli PROCESS per layer per
            % iteration -- several hundred on a small LQN, each paying the
            % binary's multi-second static initialization -- and the engine's own
            % LAYERED solver was never exercised at all, so the row measured a
            % different thing than its name says.
            %
            % An explicit argument still WINS over the ambient setting, and the
            % last one given wins, so SolverLN(model,f,'matlab') under a global
            % 'cpp' stays MATLAB.
            global LINEDefaultLang
            lang = '';
            for argIdx = 1:numel(allArgs)
                thisArg = allArgs{argIdx};
                if ischar(thisArg) && any(strcmpi(thisArg, {'matlab','java','python','cpp'}))
                    lang = lower(thisArg);
                elseif isstruct(thisArg) && isfield(thisArg,'lang') && ~isempty(thisArg.lang)
                    lang = lower(char(thisArg.lang));
                end
            end
            if isempty(lang)
                if isempty(LINEDefaultLang)
                    LINEDefaultLang = getenv('LINE_DEFAULT_LANG');
                end
                lang = lower(char(LINEDefaultLang));
            end
            wantsPython = strcmpi(lang, 'python');
            wantsCpp = strcmpi(lang, 'cpp');

            if strcmpi(lang, 'java')
                % THE OPTIONS ARRIVE IN THE THIRD SLOT and THE LAYER FACTORY IN
                % THE SECOND, on the documented form
                % SolverLN(model, solverFactory, options); this branch read
                % NEITHER. It built the JAR solver from the model alone, so the
                % JAR fell back to its DefaultSolverFactory -- SolverMVA at
                % every layer -- and LN(model, @(l)NC(l,...)) under lang='java'
                % answered with MVA layers. On lcq_threehosts that is the cache
                % hit probability 0.5 against the exact 0.48331, because the
                % cache layer is what the substituted solver was substituted in.
                % The python and cpp branches below already resolve both; this
                % is the same resolution.
                % ALL FOUR CALL SHAPES, as the native branch below resolves
                % them: on SolverLN(model,'verbose',false) the option NAME sits
                % in solverFactory and only its VALUE in varargin, so reading
                % varargin alone hands parseOptions the list {false} and
                % `Solver.m @ line 299: Invalid parameter` kills the row.
                if nargin == 1
                    self.setOptions(SolverLN.defaultOptions);
                elseif isstruct(solverFactory)
                    self.setOptions(solverFactory);
                elseif nargin > 2
                    if ischar(solverFactory)
                        inputvar = [{solverFactory}, varargin];
                    else
                        inputvar = varargin;
                    end
                    self.setOptions(Solver.parseOptions(inputvar, SolverLN.defaultOptions));
                else
                    self.setOptions(SolverLN.defaultOptions);
                end
                self.options.lang = 'java';
                if isa(solverFactory,'function_handle')
                    % Kept unconstructed, exactly as the python and cpp
                    % branches do: it is read only to resolve which layer
                    % engine the JAR SolverLN must be given.
                    self.solverFactory = solverFactory;
                end
                self.obj = JLINE.SolverLN(JLINE.from_line_layered_network(lqnmodel), ...
                    self.options, JLINE.lnLayerSolverType(self));
                self.obj.options.verbose = jline.VerboseLevel.SILENT;
                % see _kb/06-solver-catalog.md (LN section) for rationale
                self.lqn = lqnmodel.getStruct();
            elseif wantsPython
                % see _kb/06-solver-catalog.md (LN section) for rationale
                self.pyMode = true;
                % THE OPTIONS ARRIVE IN THE THIRD SLOT on the documented form
                % SolverLN(model, solverFactory, options), and reading only the
                % SECOND one discarded them whole: `options.config.layering`,
                % `options.method` and the convergence knobs never reached the
                % bridge at all, so widening PYLINE.acceptedKwargs alone was
                % INERT. Resolve them the way the native branch below does.
                % ALL FOUR CALL SHAPES, as the native branch below resolves
                % them: on SolverLN(model,'verbose',false) the option NAME sits
                % in solverFactory and only its VALUE in varargin, so reading
                % varargin alone hands parseOptions the list {false} and
                % `Solver.m @ line 299: Invalid parameter` kills the row.
                if nargin == 1
                    self.setOptions(SolverLN.defaultOptions);
                elseif isstruct(solverFactory)
                    self.setOptions(solverFactory);
                elseif nargin > 2
                    if ischar(solverFactory)
                        inputvar = [{solverFactory}, varargin];
                    else
                        inputvar = varargin;
                    end
                    self.setOptions(Solver.parseOptions(inputvar, SolverLN.defaultOptions));
                else
                    self.setOptions(SolverLN.defaultOptions);
                end
                if isa(solverFactory,'function_handle')
                    % Kept unconstructed, exactly as the cpp branch below does:
                    % it is read only to resolve which layer engine the native
                    % SolverLN must be given. Dropping it made
                    % LN(model,@(l)NC(l,...)) answer with the native default
                    % (AMVA) layers.
                    self.solverFactory = solverFactory;
                end
                self.options.lang = 'python';
                self.lqn = lqnmodel.getStruct();
            elseif wantsCpp
                % ONE subprocess solves the whole ensemble, so no layers are
                % built here: letting the MATLAB fixed point run and dispatching
                % each layer separately would spawn one process per layer per
                % iteration. see _kb/06-solver-catalog.md (LN section)
                self.cppMode = true;
                % THE OPTIONS ARRIVE IN THE THIRD SLOT on the documented form
                % SolverLN(model, solverFactory, options), exactly as in the
                % python branch above, and reading only the SECOND one replaced
                % the whole struct with the defaults: config.layering never
                % reached CPPLINE.lnKnobs, so lang='cpp' silently answered the
                % srvn model whatever was asked for.
                % ALL FOUR CALL SHAPES, as the native branch below resolves
                % them: on SolverLN(model,'verbose',false) the option NAME sits
                % in solverFactory and only its VALUE in varargin, so reading
                % varargin alone hands parseOptions the list {false} and
                % `Solver.m @ line 299: Invalid parameter` kills the row.
                if nargin == 1
                    self.setOptions(SolverLN.defaultOptions);
                elseif isstruct(solverFactory)
                    self.setOptions(solverFactory);
                elseif nargin > 2
                    if ischar(solverFactory)
                        inputvar = [{solverFactory}, varargin];
                    else
                        inputvar = varargin;
                    end
                    self.setOptions(Solver.parseOptions(inputvar, SolverLN.defaultOptions));
                else
                    self.setOptions(SolverLN.defaultOptions);
                end
                self.options.lang = 'cpp';
                if isa(solverFactory,'function_handle')
                    % Kept unconstructed: it is read only to resolve which layer
                    % engine the C++ --layer-solver must name.
                    self.solverFactory = solverFactory;
                end
                self.lqn = lqnmodel.getStruct();
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
                % see _kb/06-solver-catalog.md (LN section) for rationale
                self.lqn = lqn_fwd_rendezvous(self.lqn);

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

                % stored before construct(): buildLayers gates routed call groups
                % on the layer solver, so the handle has to be readable by then
                self.solverFactory = solverFactory;
                self.construct();
                line_debug('LN: solver factory=%s, constructing layers', func2str(solverFactory));
                for e=1:self.getNumberOfModels
                    % see _kb/06-solver-catalog.md (LN section) for rationale
                    % A setup no longer forces the MAM decomposition on the layer:
                    % the cold start is charged to the entry by lqn_setup_charge,
                    % not wired into the station, so the layer is an ordinary one
                    % and the user's own layer solver serves it.
                    layerFactory = solverFactory;
                    layerSolver = layerFactory(self.ensemble{e});
                    self.assertLayerSolverSupportsModel(layerSolver, self.ensemble{e}, e);
                    self.setSolver(layerSolver,e);
                end
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
            % one cell per ENTRY: nentries is a scalar, so length() of it was 1
            % and the table only reached its full size because the moment3 pass
            % grew it on assignment
            self.entrycdfrespt = cell(max(1,self.lqn.nentries),1);
            self.hasconverged = false;
            self.momentPassDone = false;

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
            line_debug('LN construct: built %d layers from LQN model (%d hosts, %d tasks, %d entries, %d activities)', ...
                length(self.ensemble), self.lqn.nhosts, self.lqn.ntasks, self.lqn.nentries, self.lqn.nacts);
            self.njobsorig = self.njobs;
            self.nlayers = length(self.ensemble);

            % interlock data structures are built in init() via initInterlock()

            % layering generates update maps that we use here to cache the elements that need reset
            self.routereset = unique(self.idxhash(self.route_prob_updmap(:,1)))';
            self.svcreset = unique(self.idxhash(self.thinkt_classes_updmap(:,1)))';
            self.svcreset = union(self.svcreset,unique(self.idxhash(self.call_classes_updmap(:,1)))');
        end

        function self = reset(self)
            % no-op
        end

        bool = converged(self, it); % convergence test at iteration it
        bool = convergedStoch(self, it); % convergence test for stochastic layer solvers (Robbins-Monro mode)

        function init(self) % operations before starting to iterate
            % INIT() % OPERATIONS BEFORE STARTING TO ITERATE
            % The moment3 pass is terminal WITHIN ONE SOLVE, so the flag is
            % scoped to one iterate(). Left standing across solves, the terminal
            % test in converged() fires at it=0 on the NEXT one, the loop body
            % never runs, and every metric comes back ZERO -- which a caller
            % meets simply by asking for the table twice, or by calling
            % getCdfRespT and then the table again. See BUGS.md BUG-97.
            self.momentPassDone = false;
            self.unique_route_prob_updmap = unique(self.route_prob_updmap(:,1))';
            self.tput = zeros(self.lqn.nidx,1);
            self.tputproc = cell(self.lqn.nidx,1);
            self.util = zeros(self.lqn.nidx,1);
            self.servt = zeros(self.lqn.nidx,1);
            self.servtmatrix = getEntryServiceMatrix(self);

            % see _kb/06-solver-catalog.md (LN section) for rationale
            for e= 1:self.nlayers
                self.solvers{e}.enableChecks=true;
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
            line_debug('LN init: %d layers, relaxation=%s (omega=%.3f)', ...
                self.nlayers, self.options.config.relax, self.relax_omega);
            self.servt_prev = NaN(self.lqn.nidx, 1);
            self.residt_prev = NaN(self.lqn.nidx, 1);
            self.tput_prev = NaN(self.lqn.nidx, 1);
            self.thinkt_prev = NaN(self.lqn.nidx, 1);
            self.thinkt = zeros(self.lqn.nidx, 1); % Initialize to zeros (Python parity)
            self.callservt_prev = NaN(self.lqn.ncalls, 1);
            self.callresidt_prev = NaN(self.lqn.ncalls, 1);

            % Build interlock tables (LQNS V5 static analysis)
            if self.options.config.interlocking
                self.initInterlock();
            end

            % Initialize MOL-specific state
            self.util_prev_host = zeros(self.lqn.nhosts, 1);
            self.util_prev_task = zeros(self.lqn.ntasks, 1);

            % see _kb/06-solver-catalog.md (LN section) for rationale
            initMode = 'none';
            if isfield(self.options.config, 'layer_init') && ~isempty(self.options.config.layer_init)
                initMode = self.options.config.layer_init;
            end
            if any(strcmpi(initMode, {'bound','boxbound','mwba'}))
                try
                    bnd = lqn_boxbounds(self.lqn);
                    for idx = 1:self.lqn.nidx
                        u = bnd.TN_up(idx); l = bnd.TN_lo(idx);
                        if isfinite(u) && isfinite(l) && u > 0 && l > 0
                            x = sqrt(u*l);
                        elseif isfinite(u)
                            x = u;
                        elseif isfinite(l)
                            x = l;
                        else
                            x = 0;
                        end
                        if x > 0
                            self.tput(idx) = x;
                            self.tputproc{idx} = Exp.fitRate(x);
                        end
                    end
                    line_debug('LN init: throughputs initialized from robust box bounds');
                catch ME
                    line_debug('LN box-bound initialization skipped: %s', ME.message);
                end
            end

            % see _kb/06-solver-catalog.md (LN section) for rationale
            self.stochlayers = false(1, self.nlayers);
            for e = 1:self.nlayers
                self.stochlayers(e) = self.solvers{e}.isStochastic();
            end
            mode = self.options.config.stochiter;
            self.stochiter_auto = strcmpi(mode, 'auto');
            if self.stochiter_auto
                if any(self.stochlayers)
                    mode = 'rm';
                else
                    mode = 'off';
                end
            end
            self.stochiter_mode = lower(mode);
            self.stochiter_start = [];
            self.stoch_avg = cell(1, self.nlayers);
            self.stoch_avg_count = 0;
            self.stoch_servt_avg = [];
            self.stoch_residt_avg = [];
            line_debug('LN init: stochastic iteration mode=%s (%d stochastic layers)', ...
                self.stochiter_mode, sum(self.stochlayers));
        end


        function pre(self, it) % operations before an iteration
            % PRE(IT) % OPERATIONS BEFORE AN ITERATION
            % Seed control for stochastic layer solvers
            if isempty(self.stochiter_mode)
                return
            end
            switch self.stochiter_mode
                case 'rm'
                    % see _kb/06-solver-catalog.md (LN section) for rationale
                    for e = find(self.stochlayers)
                        self.solvers{e}.options.seed = self.options.seed + (it-1)*self.nlayers + e;
                    end
                case 'crn'
                    % see _kb/06-solver-catalog.md (LN section) for rationale
                    for e = find(self.stochlayers)
                        self.solvers{e}.options.seed = self.options.seed + e;
                    end
            end
        end

        function [result, runtime] = analyze(self, it, e)
            % [RESULT, RUNTIME] = ANALYZE(IT, E)
            T0 = tic;
            line_debug('LN analyze: iteration %d, layer %d (%s)', it, e, class(self.solvers{e}));
            result = struct();
            %jresult = struct();
            % A layer solver is SILENT (see setSolver), so it prints no banner
            % of its own and there is nothing left here to separate: the blank
            % line that used to precede layer 1 is gone with the banner.

            % Protection for unstable queues during LN iterations
            % If a solver fails (e.g., due to queue instability with open arrivals),
            % use results from previous iteration if available and continue
            try
                [result.QN, result.UN, result.RN, result.TN, result.AN, result.WN] = self.solvers{e}.getAvg();
            catch ME
                if it > 1 && ~isempty(self.results) && size(self.results, 1) >= (it-1) && size(self.results, 2) >= e
                    % LN'S OWN verbosity decides, not the layer's: the layer is
                    % always SILENT now, and gating on it swallowed the warning
                    if self.options.verbose ~= VerboseLevel.SILENT
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
            % see _kb/06-solver-catalog.md (LN section) for rationale
            if it == 1 && ~isempty(self.stochlayers)
                self.stochlayers(e) = self.solvers{e}.isStochastic();
            end
            % see _kb/06-solver-catalog.md (LN section) for rationale
            if strcmp(self.solvers{e}.name, 'SolverMVA')
                sne = self.ensemble{e}.getStruct(false);
                QNe = result.QN;
                QNe(~isfinite(QNe)) = 0;
                if size(QNe,1) == sne.nstations && size(QNe,2) == sne.nclasses
                    Qch = zeros(sne.nstations, sne.nchains);
                    for c = 1:sne.nchains
                        Qch(:,c) = sum(QNe(:,sne.chains(c,:)>0),2);
                    end
                    self.solvers{e}.options.init_sol = Qch;
                end
            end
            % A fluid layer integrates its ODEs from the cold placement built
            % by SOLVER_FLUID_INITSOL, so without this every outer iteration
            % re-runs the whole transient. Carry the terminal state forward:
            % SOLVER_FLUID_ANALYZER rebuilds it when the layer's phase
            % expansion moved under it.
            if strcmp(self.solvers{e}.name, 'SolverFLD') ...
                    && (~isfield(self.options.config, 'fluid_warmstart') ...
                    || self.options.config.fluid_warmstart)
                lastSol = self.solvers{e}.result;
                if isstruct(lastSol) && isfield(lastSol, 'solverSpecific') ...
                        && isstruct(lastSol.solverSpecific) ...
                        && isfield(lastSol.solverSpecific, 'odeStateVec') ...
                        && ~isempty(lastSol.solverSpecific.odeStateVec)
                    self.solvers{e}.options.init_sol = lastSol.solverSpecific.odeStateVec(:);
                end
            end
            runtime = toc(T0);
        end

        function post(self, it) % operations after an iteration
            % POST(IT) % OPERATIONS AFTER AN ITERATION
            line_debug('LN post: iteration %d, updating metrics and layer parameters', it);
            % convert the results of QNs into layer metrics

            self.updateMetrics(it);

            if self.options.config.interlocking
                % apply interlock correction to call residence times
                self.updatePopulations(it);
            end

            % recompute think times
            self.updateThinkTimes(it);

            % update the model parameters
            self.updateLayers(it);

            % update entry selection and cache routing probabilities within callers
            self.updateRoutingProbabilities(it);

            % reset all layers with routing probability changes
            for e= self.routereset
                self.ensemble{e}.refreshChains();
                % refreshChains can change the chain basis, invalidating the
                % warm-start solution cached by analyze()
                self.solvers{e}.options.init_sol = [];
                % see _kb/06-solver-catalog.md (LN section) for rationale
                self.solvers{e}.reset();
            end

            % refresh visits and network model parameters
            for e= self.svcreset
                switch self.solvers{e}.name
                    case {'SolverMVA', 'SolverNC'} %leaner than refreshProcesses, no need to refresh phases
                        % see _kb/06-solver-catalog.md (LN section) for rationale
                        switch self.lnmethod
                            case {'moment3','srvn.ph','flat.ph'}
                                % both carry a phase-type service law, whose
                                % phases a rate-only refresh would drop
                                self.ensemble{e}.refreshProcesses();
                            otherwise
                                self.ensemble{e}.refreshRates();
                        end
                    otherwise
                        self.ensemble{e}.refreshProcesses();
                end
                self.solvers{e}.reset(); % commenting this out des not seem to produce a problem, but it goes faster with it
            end

            % Note: interlock correction is done via callresidt adjustment
            % in updatePopulations, no population changes needed

            if it==1
                % now disable all solver support checks for future iterations
                for e=1:length(self.ensemble)
                    self.solvers{e}.setChecks(false);
                end
            end
        end


        function finish(self) % operations after iterations are completed
            % FINISH() % OPERATIONS AFTER INTERATIONS ARE COMPLETED
            line_debug('LN finish: final analysis of %d layers', size(self.results,2));
            E = size(self.results,2);
            % In Robbins-Monro mode, report the Polyak-Ruppert averaged
            % results rather than the last (noisy) iterate
            if ~isempty(self.stochiter_mode) && strcmp(self.stochiter_mode,'rm') && self.stoch_avg_count > 0
                for e = 1:E
                    fnames = fieldnames(self.stoch_avg{e});
                    for f = 1:length(fnames)
                        self.results{end,e}.(fnames{f}) = self.stoch_avg{e}.(fnames{f});
                    end
                end
                self.servt = self.stoch_servt_avg;
                self.residt = self.stoch_residt_avg;
            end
            for e=1:E
                s = self.solvers{e};
                s.getAvg();
                self.solvers{e} = s;
            end
            self.model.ensemble = self.ensemble;
        end

        function [QNlqn_t, UNlqn_t, TNlqn_t] = getTranAvg(self, Qt, Ut, Tt)
            % [QNLQN_T,UNLQN_T,TNLQN_T] = GETTRANAVG(SELF,QT,UT,TT)
            % Block-diagonal aggregate transient over the LQN layers.
            %
            % options.config.ln_transient selects the inter-layer coupling of
            % the transient:
            %   'decoupled' - freeze inter-layer demands at the converged fixed
            %       point (getAvg) and run each layer's transient in isolation.
            %   'coupled'   - reconcile the per-layer transients by waveform
            %       relaxation, so layer populations and inter-layer demands
            %       co-evolve in model time (getTranAvgCoupled).
            % Both modes return the SAME block-diagonal layout; iteration 0 of
            % the coupled relaxation is exactly the decoupled result.
            if nargin < 2, Qt = []; end
            if nargin < 3, Ut = []; end
            if nargin < 4, Tt = []; end
            mode = 'coupled'; % default: waveform-relaxation coupled transient
            if isfield(self.options,'config') && isfield(self.options.config,'ln_transient') ...
                    && ~isempty(self.options.config.ln_transient)
                mode = self.options.config.ln_transient;
            end
            switch lower(mode)
                case 'coupled'
                    [QNlqn_t, UNlqn_t, TNlqn_t] = self.getTranAvgCoupled(Qt, Ut, Tt);
                case 'decoupled'
                    [QNlqn_t, UNlqn_t, TNlqn_t] = self.getTranAvgDecoupled(Qt, Ut, Tt);
                otherwise
                    line_error(mfilename, sprintf('Unknown ln_transient mode ''%s'' (use ''coupled'' or ''decoupled'').', mode));
            end
        end

        function [QNlqn_t, UNlqn_t, TNlqn_t] = getTranAvgDecoupled(self, Qt, Ut, Tt)
            % [QNLQN_T,UNLQN_T,TNLQN_T] = GETTRANAVGDECOUPLED(SELF,QT,UT,TT)
            % Decoupled (frozen-demand) layered transient. The optional
            % (Qt,Ut,Tt) handles only fix the aggregate M x K layout, which the
            % per-layer block concatenation already reproduces. %#ok<INUSD>
            self.getAvg;
            QNclass_t = {};
            UNclass_t = {};
            TNclass_t = {};
            QNlqn_t = cell(0,0);
            % see _kb/06-solver-catalog.md (LN section) for rationale
            hasTs = isfield(self.options,'timespan') && numel(self.options.timespan)>=2 ...
                && all(isfinite(self.options.timespan));
            for e=1:self.nlayers
                [crows, ccols] = size(QNlqn_t);
                s = self.solvers{e};
                if hasTs
                    savedTs = s.options.timespan;
                    s.options.timespan = self.options.timespan;
                end
                [QNclass_t{e}, UNclass_t{e}, TNclass_t{e}] = s.getTranAvg();
                if hasTs
                    s.options.timespan = savedTs;
                end
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
                % The distribution pass reads the routing encoding of the
                % activity graph, which srvn.ph layers do not carry: re-running
                % getAvg over them would reconstruct the wrong topology rather
                % than a coarser answer. Refuse by name.
                if self.isPHEncoding()
                    line_error(mfilename, sprintf(['getCdfRespT needs the routing encoding of ' ...
                        'the activity graph, which method=''%s'' does not build. Rebuild the ' ...
                        'solver with method=''srvn.cs'' or method=''moment3''.'], self.lnmethod));
                end
                % save user-specified method to temporary variable
                curMethod = self.getOptions.method;
                curLnMethod = self.lnmethod;
                % Run with moment 3. BOTH the option and the RESOLVED method have
                % to move: updateMetrics dispatches on self.lnmethod, which
                % buildLayers resolved once, so setting options.method alone left
                % the mean-based update in place and returned an EMPTY table. The
                % routing layers already built serve moment3 unchanged, so only
                % the update pass changes.
                self.options.method = 'moment3';
                self.lnmethod = 'moment3';
                self.getAvg();
                % restore user-specified method
                self.options.method = curMethod;
                self.lnmethod = curLnMethod;
            end
            cdfRespT = self.entrycdfrespt;
        end

        function varargout = getAvgTable(self, varargin)
            % [AVGTABLE,QT,UT,RT,WT,TT] = GETAVGTABLE(USELQNSNAMING)
            % The result recorder captures the returned table together with the solver
            % that produced it -- see LineResultRecorder. Recording an ensemble here
            % rather than in the member solver it delegates to is what keeps an
            % AUTO/LN/ENV/UQ answer from being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getAvgTable_impl(varargin{:});
            LineResultRecorder.capture(scope, self, 'avg', varargout{1});
        end

        function [AvgTable,QT,UT,RT,WT,AT,TT] = getAvgTable_impl(self)
            % GETAVGTABLE_IMPL Implementation of GETAVGTABLE; see the wrapper above.
            if (GlobalConstants.DummyMode)
                [AvgTable, QT, UT, RT, TT, WT] = deal([]);
                return
            end

            boundMethod = '';
            if isfield(self.options,'method') && ischar(self.options.method) ...
                    && any(strcmp(self.options.method, {'mwba.upper','mwba.lower'}))
                boundMethod = self.options.method;
            end
            if ~isempty(boundMethod)
                % Majumdar-Woodside robust box bounds for the LQN
                bnd = lqn_boxbounds(self.lqn);
                nidx = self.lqn.nidx;
                if strcmp(boundMethod,'mwba.upper')
                    TN = bnd.TN_up; UN = bnd.UN_up;
                else
                    TN = bnd.TN_lo; UN = bnd.UN_lo;
                end
                TN(isnan(TN)) = 0; UN(isnan(UN)) = 0;
                QN = zeros(nidx,1); RN = zeros(nidx,1);
                WN = zeros(nidx,1); AN = zeros(nidx,1);
            elseif ~isempty(self.obj)
                avgTable = self.obj.getEnsembleAvg();
                [QN,UN,RN,WN,AN,TN] = JLINE.arrayListToResults(avgTable);
            elseif ~isempty(self.pyMode) && self.pyMode
                % the layer solver this ensemble was built with decides the
                % bridge's layer-solver class; see PYLINE.SolverLN. Under
                % pyMode no layer is constructed, so the class is read off the
                % factory the same way CPPLINE.lnLayerSolver reads it.
                layerSolverName = PYLINE.lnLayerSolverName(self);
                [QN,UN,RN,TN,AN,WN] = PYLINE.getEnsembleAvg(self.model, self.options, self.lqn.names, layerSolverName);
            elseif ~isempty(self.cppMode) && self.cppMode
                [QN,UN,RN,TN,AN,WN] = CPPLINE.getEnsembleAvg(self, self.options);
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

        function [AvgTable,QT,UT,RT,WT,AT,TT] = avgTable(self)
            % AVGTABLE Alias for getAvgTable
            [AvgTable,QT,UT,RT,WT,AT,TT] = self.getAvgTable();
        end

        function [AvgTable,QT,UT,RT,WT,AT,TT] = avgT(self)
            % AVGT Short alias for getAvgTable
            [AvgTable,QT,UT,RT,WT,AT,TT] = self.getAvgTable();
        end

        function [AvgTable,QT,UT,RT,WT,AT,TT] = aT(self)
            % AT Short alias for getAvgTable (MATLAB-compatible)
            [AvgTable,QT,UT,RT,WT,AT,TT] = self.getAvgTable();
        end
    end

    methods
        [QN,UN,RN,TN,AN,WN] = getEnsembleAvg(self);
        [QN,UN,RN,TN,AN,WN] = getEnsembleAvgPH(self);
        [QNlqn_t, UNlqn_t, TNlqn_t] = getTranAvgCoupled(self, Qt, Ut, Tt);

        function [bool, featSupported] = supports(self, model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(SELF, MODEL)
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
        buildLayersRecursive(self, idxSet, callers, ishostlayer, flat);
        ok = buildLayersPH(self, mode, flat);  % mode 'probe' answers feasibility only
        phComposeEntryLaws(self);
        initInterlock(self);
        updateLayers(self, it);
        updateLayersPH(self, it);
        updatePopulations(self, it);
        updateThinkTimes(self, it);
        updateThinkTimesPH(self, it);
        updateMetrics(self, it);
        updateMetricsPH(self, it);
        updateRoutingProbabilities(self, it);
        svcmatrix = getEntryServiceMatrix(self)
        prOt = overtake_prob(self, eidx);  % Phase-2 overtaking probability

        function tf = isPHEncoding(self)
            % True when the layers carry the COMPOSED phase-type server law
            % rather than the routing encoding of the activity graph, under
            % either layering. The encoding, not the layering, decides which
            % update and reconstruction passes run, so every such dispatch asks
            % this and not for one method name.
            tf = any(strcmp(self.lnmethod, {'srvn.ph','flat.ph'}));
        end

        function [e, sIdx] = layerOf(self, idx)
            % [E, SIDX] = LAYEROF(IDX) returns the ensemble index of the layer
            % where the LQN element IDX acts as a server, and the index of its
            % station within that layer. Under 'flat' layering every server
            % lives in the same layer, so SIDX is what tells them apart.
            e = self.idxhash(idx);
            if isnan(e)
                sIdx = NaN;
                return
            end
            sIdx = self.stationIdxOf(e, idx);
        end

        function sIdx = stationIdxOf(self, e, idx)
            % Station index of the LQN element IDX inside layer E, falling back
            % to the layer's own server when IDX is not a server there
            attr = self.ensemble{e}.attribute;
            sIdx = attr.serverIdx;
            if isfield(attr,'serverIdxOf') && ~isnan(idx) && idx >= 1 && ...
                    idx <= length(attr.serverIdxOf) && ~isnan(attr.serverIdxOf(idx))
                sIdx = attr.serverIdxOf(idx);
            end
        end

        function sIdx = stationIdxOfClass(self, e, c)
            % Station of layer E that class C is served at: the processor of an
            % activity, the called task of a call, the layer's server otherwise
            attr = self.ensemble{e}.classes{c}.attribute;
            elem = NaN;
            if ~isempty(attr)
                switch attr(1)
                    case LayeredNetworkElement.ACTIVITY
                        elem = self.lqn.parent(self.lqn.parent(attr(2)));
                    case LayeredNetworkElement.CALL
                        elem = self.lqn.parent(self.lqn.callpair(attr(2),2));
                end
            end
            sIdx = self.stationIdxOf(e, elem);
        end

        function rows = serverStationsOf(self, e, ishost)
            % Station indices of the host (ISHOST true) or task servers of layer E
            if ishost
                rows = self.ensemble{e}.attribute.hostStations;
            else
                rows = self.ensemble{e}.attribute.taskStations;
            end
            rows = rows(:)';
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
            %   solver2 = SolverLN(model, @(m) SolverNC(m));
            %   solver2.set_state(state);
            %   solver2.getEnsembleAvg();  % Continues from MVA solution
            %
            %   NOTE: a pure-JMT layer factory (@(m) SolverJMT(m)) is
            %   unsupported because LN server layers carry immediate feedback
            %   (sn.immfeed), which SolverJMT rejects; SolverLN raises a clear
            %   error upfront in that case.

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
            state.callresidt_prev = self.callresidt_prev;

            % Results from last iteration
            state.results = self.results;

            % Interlock data
            state.njobs = self.njobs;
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
            % provide initial estimates and accurate solvers (JMT, LDES)
            % refine the solution.
            %
            % Example:
            %   solver1 = SolverLN(model, @(m) SolverMVA(m));
            %   solver1.getEnsembleAvg();
            %   state = solver1.get_state();
            %
            %   solver2 = SolverLN(model, @(m) SolverNC(m));
            %   solver2.set_state(state);
            %   solver2.getEnsembleAvg();  % Continues from MVA solution
            %
            %   NOTE: a pure-JMT layer factory (@(m) SolverJMT(m)) is
            %   unsupported because LN server layers carry immediate feedback
            %   (sn.immfeed), which SolverJMT rejects; SolverLN raises a clear
            %   error upfront in that case.

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
            if isfield(state, 'callresidt_prev')
                self.callresidt_prev = state.callresidt_prev;
            end

            % Results
            if isfield(state, 'results')
                self.results = state.results;
            end

            % Interlock data
            if isfield(state, 'njobs')
                self.njobs = state.njobs;
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
                % Ensure layer model struct is fully built before refresh
                self.ensemble{e}.getStruct(false);
                self.ensemble{e}.refreshChains();
                % refreshChains can change the chain basis, invalidating the
                % warm-start solution cached by analyze()
                self.solvers{e}.options.init_sol = [];
                if isa(self.solvers{e},'SolverMVA')
                    self.solvers{e}.resetForkWarmStart();
                end
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
            % (e.g., from MVA to LDES) while preserving the current
            % solution state.
            %
            % Example:
            %   solver = SolverLN(model, @(m) SolverMVA(m));
            %   solver.getEnsembleAvg();  % Fast initial solution
            %
            %   solver.update_solver(@(m) SolverLDES(m, 'samples', 1e5));
            %   solver.getEnsembleAvg();  % Refine with simulation
            %
            %   NOTE: a pure-JMT layer factory (@(m) SolverJMT(m, ...)) is
            %   unsupported: LN server layers carry immediate feedback
            %   (sn.immfeed), which SolverJMT rejects. update_solver raises a
            %   clear error upfront in that case.

            self.solverFactory = solverFactory;

            % Replace all layer solvers
            for e = 1:self.nlayers
                layerSolver = solverFactory(self.ensemble{e});
                self.assertLayerSolverSupportsModel(layerSolver, self.ensemble{e}, e);
                self.setSolver(layerSolver, e);
            end
        end

        function assertLayerSolverSupportsModel(self, layerSolver, layerModel, e) %#ok<INUSL>
            % ASSERTLAYERSOLVERSUPPORTSMODEL Reject a layer solver that cannot
            % represent its layer model, upfront with a clear message.
            %
            % LN server-layer stations carry immediate feedback (sn.immfeed):
            % successive same-host activities retain the server, modelled as
            % immediate-feedback self-loops so the layer solver does not
            % re-queue the job. SolverJMT rejects any model with immfeed and
            % returns no solution, so a pure-JMT layer factory otherwise fails
            % cryptically at layer 1, iteration 1. Detect it here instead.
            %
            % The guard is CONDITIONAL: it fires only when the specific layer
            % model actually carries immfeed. A SolverJMT layer solver on an
            % immfeed-free layer is allowed, and non-JMT factories are never
            % rejected.
            % A layer carrying an admission constraint needs a Region-capable
            % solver; the check is by feature set, not by solver name.
            lsnRegion = layerModel.getStruct();
            if isfield(lsnRegion,'nregions') && ~isempty(lsnRegion.nregions) && lsnRegion.nregions > 0
                if ~layerSolver.supports(layerModel)
                    line_error(mfilename, '%s cannot solve LN layer %d: the layer carries an admission constraint (finite capacity region) that its feature set does not cover. Use a layer solver that declares Region, such as SolverCTMC, SolverLDES or SolverSSA.', class(layerSolver), e);
                end
            end
            if isa(layerSolver, 'SolverJMT')
                lsn = layerModel.getStruct();
                if isfield(lsn,'immfeed') && ~isempty(lsn.immfeed) && any(lsn.immfeed(:))
                    line_error(mfilename, ['SolverJMT cannot solve LN layer %d: the layer carries immediate feedback (sn.immfeed), which SolverJMT does not support, so LN would fail at the first iteration. Use the default layer factory (MVA/NC) or another layer solver that supports immediate feedback.'], e);
                end
            end
        end

        function solver = setSolver(self, solver, e)
            % SOLVER = SETSOLVER(SOLVER, E) registers a layer solver, silenced
            %
            % A LAYER SOLVER NEVER NARRATES. The fixed point runs every layer
            % once per iteration, so a layer left at the caller's verbosity
            % prints its own banner nlayers*iter_max times and buries the
            % layered narration that the caller actually asked for. The level
            % is stamped HERE rather than in the factory because a factory the
            % user supplied -- SolverLN(model, @(m) SolverNC(m)) -- never sees
            % the LN options at all, and stamping it in the default factory
            % alone left exactly that case loud.
            %
            % SolverLN's own reporting is unaffected: it reads
            % self.options.verbose, not the layer's.
            if nargin < 3
                solver = setSolver@EnsembleSolver(self, SolverLN.silenced(solver));
            else
                solver = setSolver@EnsembleSolver(self, SolverLN.silenced(solver), e);
            end
        end

        function reportCompletion(self, runtime)
            % REPORTCOMPLETION(RUNTIME) writes the closing banner of an LN run
            %
            % NetworkSolver.setAvgResults prints this line for every
            % NetworkSolver, but SolverLN is an EnsembleSolver and never passes
            % through it, so an LN run used to end on its iteration summary
            % without ever saying which method had run or how long it took --
            % the one solver whose table arrived anonymous. SolverENV, the
            % other ensemble solver, already reported its own; this mirrors it.
            %
            % deferPrint, not line_printf: with the console narrating, the
            % banner is held until after the closing DONE line, exactly as the
            % NetworkSolver one is.
            if self.options.verbose == VerboseLevel.SILENT
                return
            end
            % THE RESOLVED LAYERING, not just the token the caller asked for.
            % 'srvn' and 'default' both resolve at layer-build time to
            % 'srvn.ph' or 'srvn.cs' depending on what the model needs, and
            % that choice is what the run actually made -- reporting 'default'
            % told the reader nothing. Printed 'requested/resolved', the same
            % shape NetworkSolver uses for 'default/exact' and 'default/nrm'.
            method = 'default';
            if isfield(self.options,'method') && ~isempty(self.options.method)
                method = char(self.options.method);
            end
            if ~isempty(self.lnmethod) && ~strcmp(char(self.lnmethod), method)
                method = sprintf('%s/%s', method, char(self.lnmethod));
            end
            lang = 'matlab';
            if isfield(self.options,'lang') && ~isempty(self.options.lang)
                lang = char(self.options.lang);
            end
            iter = 0;
            if ~isempty(self.results)
                iter = size(self.results,1);
            end
            LineConsole.deferPrint(['LN analysis [method: %s; type: %s; lang: %s; ' ...
                'env: %s] completed in %fs. Iterations: %d.\n'], ...
                method, line_method_type('LN', method), lang, ...
                version('-release'), runtime, iter);
        end

        function [allMethods] = listValidMethods(self)
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver.
            %
            % The list is the same eight names for every model; which of them can
            % actually encode THIS model is SUPPORTSMODELMETHOD's question, and
            % LN_METHOD_REFUSAL is where the rules live. This used to compute SN
            % and discard it, which read as if the list were being narrowed.
            allMethods = {'srvn','srvn.ph','srvn.cs','flat','flat.cs','flat.ph','moment3','default'};
        end

        function [bool, reason] = supportsModelMethod(self, method)
            % [BOOL, REASON] = SUPPORTSMODELMETHOD(METHOD)
            % The encoding rules the layer builders enforce at solve time,
            % stated here so a CALLER can see them before running.
            %
            % 'srvn.ph' and 'flat.ph' compose each entry into ONE phase-type law,
            % and several constructs have nowhere to go in that law: a forwarding
            % call whose target is not in the caller's activity graph, a routed
            % call group whose dispatch order the composition folds away, a cache
            % task, an admission constraint, a queue-dependent rate on a station
            % the composition replaces. 'flat.ph' additionally squashes every
            % layer into one network, which per-layer state (a replica, a
            % powered-down setup thread) cannot survive.
            %
            % None of these is a feature name, so none can be a feature-set
            % delta: they are properties of what the METHOD does to the model.
            % Left only in BUILDLAYERSPH they were invisible to every gate above
            % it, and model.help() offered all eight names on every layered
            % model.
            % NO SUPER CALL: SolverLN descends from EnsembleSolver, which
            % descends from Solver, and neither declares supportsModelMethod --
            % the base gate lives on NetworkSolver, which is a sibling. A model
            % feature set is a Network notion and a LayeredNetwork does not carry
            % one, so there is nothing above this to consult.
            bool = true;
            reason = '';
            if ~isempty(self.lqn)
                % the layering the run would use, since 'flat' rules differ
                layering = '';
                if isstruct(self.options) && isfield(self.options,'config') ...
                        && isstruct(self.options.config) && isfield(self.options.config,'layering')
                    layering = self.options.config.layering;
                end
                [bool, reason] = ln_method_refusal(self.lqn, method, layering);
            end
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

        function solver = silenced(solver)
            % SOLVER = SILENCED(SOLVER) sets a layer solver to SILENT
            %
            % Accepts the cell form of setSolver as well as a single solver.
            if iscell(solver)
                for k = 1:numel(solver)
                    solver{k} = SolverLN.silenced(solver{k});
                end
                return
            end
            if ~isempty(solver) && isprop(solver,'options') && isstruct(solver.options) ...
                    && isfield(solver.options,'verbose')
                solver.options.verbose = VerboseLevel.SILENT;
            end
        end
    end
end

function solver = adaptiveSolverFactory(model, parentOptions) %#ok<INUSD>
    % ADAPTIVESOLVERFACTORY - Select appropriate solver based on model characteristics
    % Use JMT for models with open classes, MVA for pure closed networks
    %
    % The caller's verbosity is NOT inherited: SolverLN.setSolver silences
    % every layer solver, this one included, so reading parentOptions.verbose
    % here only decided what to build and immediately discard.
    verbose = VerboseLevel.SILENT;

    % A layer carrying an admission constraint needs a Region-capable solver
    if layerHasRegion(model)
        solver = regionCapableLayerSolver(model, verbose);
        return
    end

    % Create MVA solver with reduced iter_max for sublayer stability (Python parity)
    solver = SolverMVA(model, 'verbose', verbose);
    solver.options.iter_max = 1000; % Cap sublayer MVA iterations
end

function tf = layerHasRegion(model)
    % TF = LAYERHASREGION(MODEL) - true if the layer carries a finite capacity region
    lsn = model.getStruct();
    tf = isfield(lsn,'nregions') && ~isempty(lsn.nregions) && lsn.nregions > 0;
end

function solver = regionCapableLayerSolver(model, verbose)
    % SOLVER = REGIONCAPABLELAYERSOLVER(MODEL, VERBOSE)
    %
    % Picks the first solver whose feature set covers the layer. The order is
    % by decreasing accuracy: CTMC is exact but state-space bound, LDES and SSA
    % simulate. Selection is by SolverX.supports so it self-corrects if another
    % solver later declares Region.
    ctors = {@SolverCTMC, @SolverLDES, @SolverSSA};
    checks = {@SolverCTMC.supports, @SolverLDES.supports, @SolverSSA.supports};
    names = {'SolverCTMC','SolverLDES','SolverSSA'};
    for k = 1:length(ctors)
        if checks{k}(model)
            solver = ctors{k}(model, 'verbose', verbose);
            return
        end
    end
    line_error(mfilename,'LN layer %s carries an admission constraint but none of %s supports it. Supply a layer solver factory explicitly.', model.getName(), strjoin(names,', '));
end
