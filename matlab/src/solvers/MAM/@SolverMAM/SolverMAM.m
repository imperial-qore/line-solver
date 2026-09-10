classdef SolverMAM < NetworkSolver
    % Matrix-Analytic and RCAT methods solver
    %
    % Implements matrix-analytic methods and RCAT for structured Markov chain analysis.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        function self = SolverMAM(model,varargin)
            % SOLVERMAM Create a Matrix-Analytic Methods solver instance
            %
            % @brief Creates a MAM solver for structured Markov chain analysis
            % @param model Network model to be analyzed via matrix-analytic methods
            % @param varargin Optional parameters (method, tolerance, etc.)
            % @return self SolverMAM instance configured for MAM analysis
            
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            self.setLang();
        end
                
        function sn = getStruct(self)
            % QN = GETSTRUCT()
            
            % Get data structure summarizing the model
            sn = self.model.getStruct(true);
        end
        
        runtime = runAnalyzer(self, options);
        RD = getCdfRespT(self, R);
        RD = getCdfPassT(self, R);

            function [allMethods] = listValidMethods(self)
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            sn = self.model.getStruct();
            % Note: Method order must match expected values in test files.
            % New methods should be added at the end to preserve index alignment.
            % SolverMAM.m's own list, MINUS the RCAT names
            % ('inap','inapplus','inapinf','exact') which moved to SolverAG;
            % supportsModelMethod redirects a caller that still asks. 'retrial'
            % names the retrial analyzer that 'default' also resolves to on a
            % retrial topology; it is dispatched, so it is advertised.
            allMethods = {'default','dec.source','dec.mmap','dec.poisson','mna','ldqbd','dec.source.mmap','bgchain','retrial'};
        end

        function featSupported = getMethodFeatureSet(self, method) %#ok<INUSL>
            % Every MAM method shares the solver-level feature envelope, apart
            % from the round-robin split that only 'mna' resolves; the genuine
            % per-method restrictions ('mna', 'ldqbd') are structural and are
            % applied in supportsModelMethod below.
            %
            % Defining this is what lets NetworkSolver.supportsModelMethod name
            % the offending features: with no method feature set it falls back
            % to the coarse supports(model) and returns an empty reason, so the
            % gate could only report "features not supported" without saying
            % which ones.
            featSupported = SolverMAM.getFeatureSet();
            % SETUP/DELAY-OFF IS SOLVER_MAM_BASIC'S ALONE (the hassetup branch
            % of dec.source / dec.poisson): no other MAM algorithm reads
            % sn.hassetup, so every other method drops it here and 'default',
            % which may route to any of them, refuses the combinations that
            % would in supportsModelMethod. FORK-JOIN is served by the
            % synchronizing analyzers only (solver_mam_fj, the mmap
            % decomposition); solver_mam_basic treats a Join as a pass-through,
            % and the mna sweep raises on a Fork, so the plain methods drop it.
            fjNames = {'Fork','Join','Forker','Joiner'};
            switch method
                case 'mna'
                    % deterministic traffic split, open models only; the
                    % closed branch has no counterpart and is rejected in
                    % supportsModelMethod. The flow sweep updates INF, PS and
                    % FCFS stations only: a priority station keeps a zero
                    % queue length (MAM_MNA_APPLICABLE says so by name).
                    featSupported.setTrue({'RoutingStrategy_RROBIN'});
                    featSupported.setFalse([{'SchedStrategy_HOL','SchedStrategy_FCFSPRPRIO', ...
                        'SetupDelayOff'}, fjNames]);
                case 'dec.poisson'
                    featSupported.setFalse(fjNames);
                case 'dec.source.mmap'
                    featSupported.setFalse({'SetupDelayOff'});
                case 'bgchain'
                    featSupported.setFalse([{'SetupDelayOff'}, fjNames]);
                case 'retrial'
                    % The MAP/M/s+G analyzer of solver_mam_retrial, which
                    % 'default' also resolves to; the shape it needs is
                    % structural (MAM_RENEGING_APPLICABLE).
                    featSupported.setTrue({'Reneging'});
                    featSupported.setFalse([{'SetupDelayOff'}, fjNames]);
                case 'dec.mmap'
                    % SOLVER_MAM is an OPEN-network departure-process fixed
                    % point: it iterates on arrival streams a closed population
                    % does not have, and its station ladder serves EXT, FCFS,
                    % HOL, FCFSPRPRIO and PS only. Both restrictions are things
                    % the model HAS, so both belong here rather than in
                    % supportsModelMethod; solver_mam.m raises the matching
                    % message when the method is named by hand. Until this
                    % delta existed the gate offered dec.mmap on a closed model
                    % and the analyzer answered with a table of zeros.
                    % Fork-join goes too: the sweep uses SOLVER_MAM_TRAFFIC, the
                    % plain traffic step, which has no synchronization -- which
                    % is why solver_mam_analyzer routes a fork-join topology
                    % from 'default'/'dec.source' to SOLVER_MAM_FJ and never
                    % here. Left declared, the gate offered dec.mmap on an open
                    % fork-join model and the departure process it then built
                    % had no recurrent state.
                    featSupported.setFalse({'ClosedClass','SelfLoopingClass', ...
                        'SchedStrategy_INF','Fork','Join','Forker','Joiner','SetupDelayOff'});
                case 'default'
                    % LoadDependence and Reneging both reach 'default' through
                    % the shapes it routes to (ldqbd, the reneging analyzer);
                    % outside those shapes supportsModelMethod refuses by name.
                    featSupported.setTrue({'LoadDependence','Reneging'});
                case 'ldqbd'
                    % SOLVER_MAM_LDQBD is the only MAM algorithm that reads
                    % sn.lldscaling: it applies the level-dependent factor in
                    % each departure block. 'default' declares it because the
                    % single-class closed Delay+Queue shape routes there, and
                    % supportsModelMethod refuses a load-dependent model outside
                    % that shape -- a featset cannot see topology, and
                    % solver_mam_basic would otherwise solve every level at the
                    % nominal rate. Every other method keeps refusing.
                    featSupported.setTrue({'LoadDependence'});
                    featSupported.setFalse([{'SetupDelayOff'}, fjNames]);
            end
            % RETRIAL (recorded since 2026-09-05) is served by solver_mam_retrial
            % alone, which 'default' and 'dec.source' route to on the BMAP/PH/N/N
            % shape (solver_mam_analyzer's shared arm) and 'retrial' names.
            % Every other algorithm reads no sn.retrial* field and would answer
            % with the refused jobs lost, so the base grant is withdrawn there;
            % supportsModelMethod refuses an orbit OUTSIDE that shape for the
            % two routing names, the same way it refuses reneging outside it.
            if ~any(strcmp(method, {'default','dec.source','retrial'}))
                featSupported.setFalse({'Retrial'});
            end
            % FINITECAPACITY (registry name since 2026-09-05): the open-network
            % analyzers carry a finite buffer as a loss buffer (solver_mam_basic
            % M/M/c/K and MMAP[K]/G/1/K, solver_mam and solver_mna_open
            % truncate-and-renormalize, solver_mam_retrial's bufferless N/N
            % station), so the base envelope declares it. The two chains that
            % read no sn.cap withdraw it: the LD-QBD levels run to the
            % population or the cutoff, and the background chain to the state
            % cap. A buffer a CLOSED class can fill is refused for every method
            % by MAM_BUFFER_REFUSAL (a closed job blocks, a loss formula does not).
            if any(strcmp(method, {'ldqbd','bgchain'}))
                featSupported.setFalse({'FiniteCapacity'});
            end
        end

        function reason = unsupportedMethodReason(self, method) %#ok<INUSL>
            % REASON = UNSUPPORTEDMETHODREASON(METHOD)
            %
            % The forwarding address for the RCAT names, which are SolverAG's
            % now. Asks nothing of the model, so runAnalyzerChecks can call it
            % before the struct is built; supportsModelMethod returns the same
            % string, so a caller gets one answer whichever gate it meets first.
            reason = '';
            if any(strcmp(method, {'inap','inapplus','inapinf','exact'}))
                reason = sprintf(['The %s method moved to SolverAG: RCAT decomposes the ' ...
                    'model into cooperating agents rather than decomposing traffic, and no ' ...
                    'MAM algorithm shares its machinery. Use SolverAG(model,''%s'').'], ...
                    method, method);
            end
        end

        function [bool, reason] = supportsModelMethod(self, method)
            % Method-aware gate for the genuine per-method restrictions of the
            % MAM analyzer (mirrors the inline guards in solver_mam_analyzer):
            % 'mna' does not support mixed open/closed models, and 'ldqbd'
            % requires a single-class model. All other methods rely on the
            % coarse MAM feature set and the analyzer's topology-based routing
            % (e.g. a Fork-Join model on 'default'/'dec.source' is routed to the
            % FJ solver, not rejected).
            sn = self.model.getStruct();
            % The RCAT methods moved to SolverAG. Name them here rather than
            % letting listValidMethods report an unknown method, so a caller
            % carrying an old options.method is told where they went. The text
            % comes from unsupportedMethodReason, which runAnalyzerChecks also
            % asks, so the gate above the dispatcher and this one cannot drift
            % into two answers.
            moved = self.unsupportedMethodReason(method);
            if ~isempty(moved)
                bool = false;
                reason = moved;
                return;
            end
            % G-network signals belong to the RCAT analyzer alone: no MAM
            % algorithm reads sn.issignal, so every one of them would solve the
            % model with the signals turned into ordinary customers and report
            % that as the answer. Refuse by name so the message says where to go.
            if isfield(sn, 'issignal') && ~isempty(sn.issignal) && any(sn.issignal)
                bool = false;
                reason = sprintf(['The %s method does not support G-network signals: no MAM ' ...
                    'algorithm reads sn.issignal, so this method would solve the model with ' ...
                    'every signal turned into an ordinary customer. Use SolverAG, whose ' ...
                    'RCAT methods are the only ones that model signals.'], method);
                return;
            end
            % A slotted model is routed to SOLVER_MAM_DT before the method name
            % is read, so its scope binds every method; the predicate is the
            % one solver_mam_dt raises with.
            isDiscreteTime = sn_is_discrete_time(sn, self.getOptions());
            if isDiscreteTime
                [dtOk, dtWhy] = mam_dt_supports(sn);
                if ~dtOk
                    bool = false;
                    reason = dtWhy;
                    return;
                end
            end
            % A finite buffer a closed class can fill: no MAM analyzer blocks,
            % so every method refuses it (the open loss case is the one the
            % envelope declares). Same predicate solver_mam_analyzer raises.
            bufferWhy = mam_buffer_refusal(sn);
            if ~isempty(bufferWhy)
                bool = false;
                reason = bufferWhy;
                return;
            end
            % A retrial orbit reaches 'default' and 'dec.source' only through
            % the BMAP/PH/N/N shape their shared arm routes to solver_mam_retrial
            % on (QSYS_IS_RETRIAL); outside it the arm falls through to
            % solver_mam_basic, which reads no sn.retrial* field and would
            % answer with the refused jobs lost. Refuse by name instead, as the
            % reneging-outside-the-shape rule below does for 'default'.
            if any(strcmp(method, {'default','dec.source'})) && mam_has_retrial(sn)
                [isRetrialShape, retInfo] = qsys_is_retrial(sn);
                if ~isRetrialShape
                    bool = false;
                    reason = sprintf(['This model declares a retrial orbit outside the BMAP/PH/N/N ' ...
                        'bufferless shape (one open class, one bufferless Queue with a RETRIAL drop ' ...
                        'rule) that the %s method solves through the retrial analyzer; dec.source ' ...
                        'would ignore the orbit. %s Use SolverCTMC, SolverSSA, SolverJMT or ' ...
                        'SolverLDES.'], method, retInfo.errorMsg);
                    return;
                end
            end
            hasSetup = isfield(sn, 'hassetup') && any(sn.hassetup);
            switch method
                case 'default'
                    % getMethodFeatureSet declares LoadDependence for 'default'
                    % because the single-class closed Delay+Queue shape routes
                    % to solver_mam_ldqbd, which reads sn.lldscaling. Outside
                    % that shape 'default' falls through to solver_mam_basic,
                    % which never reads the field and would solve every level at
                    % the nominal rate, so refuse by name here. The shape is
                    % MAM_LDQBD_APPLICABLE's closed regime, the predicate the
                    % analyzer routes on.
                    if sn_has_load_dependence(sn) && any(sn.lldscaling(:) ~= 1)
                        [ldOk, ~, ldRegime] = mam_ldqbd_applicable(sn);
                        if ~ldOk || ~strcmp(ldRegime, 'closed') || hasSetup
                            bool = false;
                            reason = ['This model uses load-dependent service rates outside the ' ...
                                'single-class closed Delay+Queue shape that the default method routes ' ...
                                'to the level-dependent QBD (solver_mam_ldqbd); the dec.source ' ...
                                'decomposition it would otherwise use does not read the scaling and ' ...
                                'would solve every level at the nominal rate. Call method ''ldqbd'' ' ...
                                'directly, which refuses by name if the shape still does not fit.'];
                            return;
                        end
                    end
                    % Reneging reaches 'default' through the MAP/M/s+G analyzer
                    % of solver_mam_retrial and nowhere else: dec.source would
                    % serve the patience as if it were absent. Same predicate
                    % the analyzer decides on.
                    if mam_has_reneging_patience(sn)
                        [renOk, renInfo] = mam_reneging_applicable(sn);
                        if ~renOk
                            bool = false;
                            reason = sprintf(['This model declares a reneging patience outside the ' ...
                                'MAP/M/s+G shape (one Source, one FCFS Queue with exponential service, ' ...
                                'one open class) that the default method solves through the reneging ' ...
                                'analyzer; dec.source would ignore the patience. %s Use SolverCTMC, ' ...
                                'SolverSSA, SolverJMT or SolverLDES.'], renInfo.errorMsg);
                            return;
                        end
                    end
                    % Setup/delay-off is served by solver_mam_basic only, and a
                    % fork-join, retrial or reneging topology routes 'default'
                    % elsewhere, where the setup would be dropped on the floor.
                    if hasSetup && (sn_has_fork_join(sn) || qsys_is_retrial(sn) ...
                            || mam_has_reneging_patience(sn))
                        bool = false;
                        reason = ['This model combines setup/delay-off times with a fork-join, ' ...
                            'retrial or reneging topology; only the dec.source decomposition ' ...
                            'honours setup times, and the default method routes such a topology ' ...
                            'to an analyzer that does not read them. Use SolverJMT or SolverLDES.'];
                        return;
                    end
                case 'mna'
                    % The correctness rules of the two mna sweeps, which the
                    % analyzer's 'mna' arm raises with; see MAM_MNA_APPLICABLE.
                    [mnaOk, mnaWhy] = mam_mna_applicable(sn);
                    if ~mnaOk
                        bool = false;
                        reason = mnaWhy;
                        return;
                    end
                case 'ldqbd'
                    % The two-station shape solver_mam_ldqbd builds its chain on,
                    % asked of the predicate it raises with.
                    [ldOk, ldWhy] = mam_ldqbd_applicable(sn);
                    if ~ldOk
                        bool = false;
                        reason = ldWhy;
                        return;
                    end
                case 'retrial'
                    % A "must be present" rule, which a feature set cannot
                    % state: it says which constructs are ACCEPTED, so it can
                    % refuse a model for having something and never for lacking
                    % it. solver_mam_retrial needs an impatience configuration
                    % to analyze; MAM_RETRIAL_APPLICABLE is the same predicate
                    % solver_mam_analyzer asks before running it.
                    [retrialOk, retrialWhy] = mam_retrial_applicable(sn);
                    if ~retrialOk
                        bool = false;
                        reason = retrialWhy;
                        return;
                    end
                case 'bgchain'
                    % The closed classes ARE the background chain (a purely
                    % closed model is its degenerate, exact case), no class
                    % priorities, no fork-join, and a chain within the state
                    % cap: MAM_BGCHAIN_APPLICABLE, which the analyzer's default
                    % chooser and solver_mam_bgchain itself also ask.
                    [bgOk, bgWhy] = mam_bgchain_applicable(sn, self.getOptions());
                    if ~bgOk
                        bool = false;
                        reason = bgWhy;
                        return;
                    end
            end
            [bool, reason] = supportsModelMethod@NetworkSolver(self, method);
        end
end
    
    methods (Static)

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()

            featSupported = SolverFeatureSet;
            % MAM features
            featSupported.setTrue({'Sink','Source',...
                'Fork','Join','Forker','Joiner',... % Fork-Join support (via FJ_codes)
                'Delay','DelayStation','Queue',...
                'APH','Coxian','Erlang','Exp','HyperExp','MMPP2','MAP','MMAP','DMAP','ME','RAP',...
                'Det','Gamma','Lognormal','Pareto','Uniform','Weibull',...
                'StatelessClassSwitcher','InfiniteServer',...
                'ClassSwitch', ...
                'SharedServer','Buffer','Dispatcher',...
                'Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS','SchedStrategy_HOL',...
                'SchedStrategy_FCFSPRPRIO',... % solver_mam_basic: MMAPPH1PRPR
                'SchedStrategy_FCFS',...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'ClosedClass','SelfLoopingClass',...
                'OpenClass'});
            % Add BMAP/PH/N/N retrial queue features. 'Retrial' is kept by the
            % names that route to solver_mam_retrial only (getMethodFeatureSet).
            featSupported.setTrue({'Retrial', 'BMAP', 'PH'});
            % c-server stations (every analyzer reads sn.nservers; the slotted
            % path's single-server rule stays structural in mam_dt_supports)
            % and finite buffers as LOSS buffers, see getMethodFeatureSet and
            % mam_buffer_refusal for the closed-class refusal.
            featSupported.setTrue({'MultiServer', 'FiniteCapacity'});
            % Discrete-time (slotted) lattice laws, solved by the Q-MAM
            % discrete-time queues; see solver_mam_dt
            featSupported.setTrue({'Geometric', 'DiscreteUniform'});
            % Setup/delay-off: open stations are solved exactly by
            % qbd_setupdelayoff, closed stations by qbd_setupdelayoff_closed on
            % the finite level-dependent chain (BUG-78, 2026-09; it replaced a
            % per-instance cold-start race that carried no queueing term). The
            % closed analysis is exact at a single server with exponential
            % service and no load dependence, which is the regime
            % isClosedDelayQueue routes to ldqbd; outside it the dec.source
            % decomposition carries the same chain approximately.
            featSupported.setTrue({'SetupDelayOff'});
        end
        
        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverMAM.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end
        
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = SolverOptions('MAM');
        end

        function libs = getLibrariesUsed(sn, options)
            % GETLIBRARIESUSED Get list of external libraries used by MAM solver
            % Detect libraries used by MAM solver based on topology and method
            libs = {};

            % MAMSolver used for matrix-analytic methods (M/G/1, GI/M/1 types)
            % This includes default and decomposition methods
            if ismember(options.method, {'default', 'dec.source', 'dec.mmap', 'dec.poisson', 'dec.source.mmap'})
                libs{end+1} = 'MAMSolver';
            end

            % Q-MAM used by the traffic-split method; the RCAT methods that
            % also used it moved to SolverAG.
            if ismember(options.method, {'mna'})
                libs{end+1} = 'Q-MAM';
            end

            % SMCSolver used for QBD (Quasi-Birth-Death) analysis
            % Currently available but not actively used in default paths
            % Uncomment when QBD methods are activated:
            % if ismember(options.method, {'qbd'})
            %     libs{end+1} = 'SMCSolver';
            % end

            % BUTools usage detection (via KPCToolbox MAP functions)
            % BUTools is used when analyzing MAP/PH distributions
            if ~isempty(sn) && isfield(sn, 'proc') && ~isempty(sn.proc) && any(~cellfun(@isempty, sn.proc(:)))
                libs{end+1} = 'BUTools';
            end

            % Remove duplicates and maintain order
            libs = unique(libs, 'stable');
        end
    end
end
