classdef SolverFLD < NetworkSolver
    % FLD - Fluid/Mean-Field Approximation solver for large-scale network analysis
    %
    % Implements fluid approximation by replacing discrete populations with continuous levels.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        function self = SolverFLD(model,varargin)
            % SOLVERFLUID Create a Fluid solver instance
            %
            % @brief Creates a Fluid solver for continuous approximation analysis
            % @param model Network model to be analyzed via fluid approximation
            % @param varargin Optional parameters (method, tolerances, etc.)
            % @return self SolverFLD instance configured for fluid analysis

            % An auxiliary solver passed as first optional argument requests a
            % warm start: its steady-state solution decides the initial state
            % of the ODE integration (see NetworkSolver.initFromSolver).
            initSolver = [];
            if ~isempty(varargin) && isa(varargin{1}, 'NetworkSolver')
                initSolver = varargin{1};
                varargin(1) = [];
            end
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            self.setLang();
            if ~isempty(initSolver)
                self.initFromSolver(initSolver);
            end
        end
    end
    
    methods
        RD = getTranCdfPassT(self, R);
        [Pnir,logPnir] = getProbAggr(self, ist);
        RD = getCdfRespT(self, R);
        RD = getCdfPassT(self, R);
        RD = getCdfPT(self, R);
        [AoI, PAoI, aoiTable] = getAvgAoI(self);
        [AoI_cdf, PAoI_cdf] = getCdfAoI(self, t_values);

        function sn = getStruct(self)
            % QN = GETSTRUCT()
            
            % Get data structure summarizing the model
            sn = self.model.getStruct(true);
        end
        
        [QNt,UNt,TNt] = getTranAvg(self,Qt,Ut,Tt);

        % export the mean-field ODE system to LaTeX (scalar or matrix notation)
        [tex, sys] = exportODEs(self, filename, notation);

        % solve method is supplied by Solver superclass
        runtime = runAnalyzer(self, options);

        function bool = supportsTransientAnalysis(self) %#ok<MANU>
            % Transient averages are available (fluid ODE integrated over options.timespan).
            bool = true;
        end

        function [allMethods] = listValidMethods(self)
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            % Every method the dispatch in runAnalyzer accepts, INCLUDING the
            % 'fluid.'-qualified spelling of each: the dispatch already took
            % 'fluid.tbi', 'fluid.dae' and six others while this list named none
            % of them, so the gate refused names the solver served. 'butools'
            % (the MFQ backend) and 'aoi' (its age-of-information reading) are
            % aliases of 'mfq', which is how native python spells them.
            % 'fluid.default', which the JAR, native python and C++ take, is
            % NOT listed: the qn sanity goldens baseline every name on this
            % list, so a new one needs its block generated first
            % (genMissingSanityBaselines), and 'default' is the same run.
            allMethods = {'default', ...
                'matrix','fluid.matrix','pnorm','fluid.pnorm', ...
                'softmin','fluid.softmin', ...
                'statedep','fluid.statedep', ...
                'closing','fluid.closing', ...
                'minnormal','fluid.minnormal', ...
                'refined','fluid.refined', ...
                'tbi','fluid.tbi', ...
                'diffusion','fluid.diffusion', ...
                'mfq','fluid.mfq','butools', ...
                'rmf','fluid.rmf', ...
                'aoi','fluid.aoi', ...
                'kp','fluid.kp', ...
                'dae','fluid.dae'};
            % The single-station fluid limits, for the Source -> Queue -> Sink
            % shape they are stated for. Listing them on any other model would
            % name a method that cannot run.
            if isa(self.model,'Network')
                sn = self.model.getStruct();
                if sn.nclasses == 1 && sn.nclosedjobs == 0 && numel(sn.nodetype) == 3 && ...
                        any(sn.nodetype == NodeType.Source) && any(sn.nodetype == NodeType.Sink)
                    % the 'fluid.'-qualified spelling too: runAnalyzer.m:270
                    % and :311 both accept it, and the JAR, python and C++ all
                    % advertise it, so omitting it refused a name every
                    % codebase serves
                    allMethods = {allMethods{:}, 'mtginf','fluid.mtginf','mol','fluid.mol'}; %#ok<CCAT>
                    queue_ist = sn.nodeToStation(sn.nodetype == NodeType.Queue);
                    if ~isempty(queue_ist) && ~isempty(sn_patience_handles(sn, queue_ist, 1))
                        % 'ggisgi' and 'tga' are the SHORT spellings the C++
                        % fluid_qsys_canonical maps onto the two primary names;
                        % accepted and advertised here too so the four lists are
                        % the same set
                        allMethods = {allMethods{:}, 'ggisgi.fluid','fluid.ggisgi','ggisgi', ...
                            'ggingi.tga','fluid.tga','tga','tvms','fluid.tvms'}; %#ok<CCAT>
                    end
                end
            end
        end

        function moments = getMoments(self)
            % MOMENTS = GETMOMENTS()
            % Second-order results of the moment-closure methods
            % ('minnormal', 'refined'): state-level covariance Sigma,
            % station-class queue-length variance QVar and standard deviation
            % QStd, per-station population variance sigma2, and the 1/N
            % refinement. Empty for every first-order method, which computes
            % no second moment at all.
            %
            % On a STOCHASTIC PETRI NET (method 'dae', SOLVER_FLUID_PETRI) the
            % same call also returns a `petri` block: the mean marking and its
            % variance as (nnodes x nclasses) matrices, the per-mode firing
            % flows, the immediate firing flows, the P-invariants with the
            % residual each one closed to, and which place capacities bound.
            if isempty(self.result)
                % PRE-EXISTING BUG, found while verifying the dae route: there
                % is no `run` method on SolverFLD or on any of its bases, so
                % getMoments on a solver that had not already been analysed
                % threw "Unrecognized method ... 'run'" rather than analysing.
                % Every other FLD getter that has to trigger the analysis uses
                % getAvg (getAvgAoI, getCdfAoI) or runAnalyzer (getProbAggr);
                % getAvg is the right one here because it is what populates
                % result.solverSpecific, which is the field read just below.
                self.getAvg();
            end
            moments = [];
            if isfield(self.result,'solverSpecific') && isstruct(self.result.solverSpecific) ...
                    && isfield(self.result.solverSpecific,'moments')
                moments = self.result.solverSpecific.moments;
            end
        end
        function featSupported = getMethodFeatureSet(self, method)
            % FLD methods share the solver-level feature envelope except for
            % load dependence, which only the closing family evaluates.
            %
            % Defining this is what lets NetworkSolver.supportsModelMethod name
            % the offending features: with no method feature set it falls back
            % to the coarse supports(model), which returns an empty reason, so
            % the gate could only report "features not supported" without
            % saying which ones.
            %
            % A non-Network model (e.g. a LayeredNetwork) has no
            % getUsedLangFeatures, so it keeps the coarse path and any
            % structural checks or redirects that operate on such models.
            if ~isa(self.model, 'Network')
                featSupported = [];
                return;
            end
            featSupported = SolverFLD.getFeatureSet();
            % EVERY TEST BELOW IS ON THE CANONICAL NAME. Enumerating each
            % spelling by hand is what let the lists drift apart: the bare
            % 'ggisgi' and 'tga' that LISTVALIDMETHODS advertises were missing
            % from the Reneging branch, so the gate refused the two names on
            % exactly the reneging models they are offered for, and 'butools'
            % and 'aoi' lost the HOL grant their alias target 'mfq' keeps.
            % Canonicalize once and the four codebases cannot disagree about a
            % spelling again.
            m = SolverFLD.canonicalMethod(method);
            % A NAME THAT RESOLVES IS GATED AS ITS RESOLUTION (see RESOLVEMETHOD).
            % RUNANALYZER resolves 'default' before the gate, so the run was
            % always gated as the concrete method; a report asking about the
            % literal name got the envelope below with none of the closing
            % family's grants, so findSolver refused 'fluid.default' on a
            % load-dependent, GPS or Petri-net model that SolverFLD(model) then
            % solved. 'mfq' and 'rmf' resolve to 'matrix' off their shape (the
            % documented fallback), so their own deltas below apply only when
            % they run as themselves.
            if any(strcmp(m, {'default','mfq','rmf'}))
                opts = self.getOptions();
                opts.method = m;
                m = SolverFLD.canonicalMethod(self.resolveMethod(opts));
            end
            % The single-station fluid limits. Three of them are stated for a
            % queue customers ABANDON, and Reneging stays out of the base FLD
            % envelope: the drift carries no abandonment flow, so every other
            % method would integrate the model as if nobody left.
            if any(strcmp(m, {'ggisgi.fluid','ggingi.tga','tvms'}))
                featSupported.setTrue({'Reneging'});
            end
            if any(strcmp(m, {'ggisgi.fluid','ggingi.tga','tvms','mtginf','mol'}))
                % Every one of them is stated for a single open station; the
                % base envelope's closed classes have no meaning there, and a
                % fork-join model is answered by a fixed point over a network,
                % which a closed form for one station is not.
                featSupported.setFalse({'ClosedClass','SelfLoopingClass', ...
                    'Fork','Join','Forker','Joiner','JoinPartial'});
            end
            if strcmp(m, 'mfq')
                % SOLVER_MFQ skips a closed class with a warning and reports
                % zeros for it. A MAPt/PHt slot is a schedule, not a {D0,D1}
                % pair, and the parameter reader took its first entry for a
                % one-phase generator: a modulated schedule ran as a Poisson
                % stream of the same mean, without a word.
                featSupported.setFalse({'ClosedClass','SelfLoopingClass','MAPt','PHt'});
            end
            % DPS is a WEIGHTED share, and only the closing family reads
            % sn.schedparam: the matrix drift is the PS aggregate, which is
            % what 'matrix' and 'pnorm' were refusing by name in RUNANALYZER
            % while this envelope admitted it; 'rmf' solves its network step
            % with that same drift, 'kp' shares min(n,c)/n equally among the
            % classes, and SOLVER_FLUID_DIFFUSION refuses DPS outright.
            if any(strcmp(m, {'matrix','pnorm','rmf','kp','diffusion'}))
                featSupported.setFalse({'SchedStrategy_DPS'});
            end
            % ODE_SOFTMIN and ODE_STATEDEP have no branch for an EXT station
            % and raise on one ('does not support open models'); the envelope
            % said otherwise and the run stopped inside the ODE builder.
            if any(strcmp(m, {'softmin','statedep'}))
                featSupported.setFalse({'OpenClass', ...
                    'Source','Sink','RandomSource','JobSink'});
            end
            % A cache model is a DECOMPOSITION (SOLVER_FLD_CACHEQN_ANALYZER),
            % whose network step carries the closure for 'minnormal' and the
            % first-order drift for 'rmf' and nothing else: 'refined' has no
            % single state to correct and 'dae' no single drift to constrain.
            % RUNANALYZER refused both by name; stated here so a report does.
            if any(strcmp(m, {'refined','dae'}))
                featSupported.setFalse({'Cache','CacheClassSwitcher', ...
                    'ReplacementStrategy_RR','ReplacementStrategy_FIFO','ReplacementStrategy_SFIFO'});
            end
            if strcmp(m, 'kp')
                % The Ko-Pender limits are proved for an OPEN network of
                % stations fed by external arrival processes: a closed class has
                % no arrival process to modulate and no source phase to carry,
                % and the cache and class-switch machinery has no counterpart in
                % the paper's event set. Narrow the envelope rather than return
                % the drift of a model the limits do not describe. Mirrors
                % native python's getMethodFeatureSet.
                featSupported.setFalse({'ClosedClass','SelfLoopingClass', ...
                    'Cache','CacheClassSwitcher','ClassSwitch','StatelessClassSwitcher', ...
                    'ReplacementStrategy_RR','ReplacementStrategy_FIFO', ...
                    'ReplacementStrategy_SFIFO'});
            end
            if strcmp(m, 'refined')
                % CLOSED MODELS ONLY, which RUNANALYZER has always enforced by
                % name and the featset never stated: FLUID_REFINE_MEANFIELD
                % solves its 1/N correction on orth(D) over the FULL state, so
                % on an open model it adds a perturbation to the SOURCE POOL
                % mass, a normalisation constant rather than a population. Only
                % 'minnormal' was validated open. Stating it here is what lets a
                % report withdraw the pair instead of offering a run that stops
                % -- on an open fork-join model the same restriction surfaced as
                % a failure inside the MMT fixed point rather than as a refusal.
                featSupported.setFalse({'OpenClass', ...
                    'Source','Sink','RandomSource','JobSink'});
            end
            if strcmp(m, 'diffusion')
                % SOLVER_FLUID_DIFFUSION integrates a stochastic differential
                % equation whose drift PROJECTS each class back onto its own
                % fixed population, which is the closed-network constraint
                % itself: an open class has no population to project onto, and
                % a Source is not a station the SDE has a coordinate for.
                featSupported.setFalse({'OpenClass', ...
                    'Source','Sink','RandomSource','JobSink'});
            end
            if any(strcmp(m, {'diffusion','kp'}))
                % NEITHER OF THESE TWO INTEGRATES A FORK-JOIN MODEL, and each
                % says so by answering rather than by refusing, which is the
                % reason to state it here. Measured on a SYMMETRIC closed
                % fork-join (Delay -> Fork -> two identical FCFS queues -> Join,
                % N = 2) whose exact chain is Q1 = Q2 = 0.664, J = 0.624,
                % D = 1.024: 'diffusion' returns the whole population on ONE
                % station and zero elsewhere -- a different station on a rerun,
                % so the SDE is not integrating this model at all -- and 'kp'
                % returns an ALL-ZERO table on a symmetric OPEN fork-join fed at
                % rate 0.5, an empty network where jobs are arriving. The C++
                % featset has always withheld the names; MATLAB and native
                % python offered them and mis-answered.
                featSupported.setFalse({'Fork','Join','Forker','Joiner','JoinPartial'});
            end
            if strcmp(m, 'tbi')
                % Trajectory-based iteration decomposes the CLOSED population
                % into cells and relaxes the waveforms between them; there is
                % no cell for an unbounded open stream. A cache model is solved
                % by decomposition rather than by one drift, so the cell
                % partition has nothing to partition -- use 'rmf'.
                featSupported.setFalse({'OpenClass', ...
                    'Source','Sink','RandomSource','JobSink', ...
                    'Cache','CacheClassSwitcher', ...
                    'ReplacementStrategy_RR','ReplacementStrategy_FIFO','ReplacementStrategy_SFIFO'});
            end
            % Limited load dependence composes with the closure as a rate
            % multiplier alpha(n_i) on the scheduling share, which only the
            % closing family evaluates (ODE_RATES_CLOSING_FACTORS). The
            % matrix, pnorm, softmin, statedep, tbi, diffusion, mfq and rmf
            % paths build their drift independently and would silently ignore
            % alpha, so they must keep rejecting it.
            if ~any(strcmp(m, {'closing','minnormal','refined','dae'}))
                featSupported.setFalse({'LoadDependence'});
            end
            % Scheduling disciplines with no branch in the drift. A station
            % without a case in ODE_RATES_CLOSING_FACTORS keeps rates = x, i.e.
            % it is integrated as an INFINITE SERVER, and the answer is wrong
            % without any warning: on Delay(Z=1) -> Queue(c=1), N=4, exact
            % Q2 = 3.0154, the fall-through returns 2.0000.
            %
            % SOLVER_FLUID_CLOSING caught LCFS/LCFSPR/HOL in its metric reader
            % ('Unsupported scheduling policy'), but SOLVER_FLUID_MOMENTS reads
            % the event representation directly and never reached that guard,
            % so `minnormal`/`refined` answered all of them silently. SIRO was
            % worse still: that reader accepts it as FCFS, so the ODE
            % integrated it as INF while the metrics were read as if it shared
            % the server, and EVERY closing-family method returned 2.0000
            % silently. Reject at the featset gate instead, where the message
            % names the offending discipline.
            %
            % Only the closing family is affected: matrix/pnorm build a PS
            % drift for every queueing station, which is the right aggregate
            % for any work-conserving discipline, and SOLVER_FLUID_DIFFUSION
            % lists SIRO explicitly.
            if any(strcmp(m, {'closing','statedep','softmin','tbi','minnormal','refined','dae'}))
                featSupported.setFalse({'SchedStrategy_SIRO', ...
                    'SchedStrategy_LCFS','SchedStrategy_LCFSPR'});
            end
            % SOLVER_FLUID_DIFFUSION lists PS, FCFS, INF and SIRO and refuses
            % every other discipline by name; the envelope admitted LCFS and
            % LCFSPR (DPS is withdrawn above), so the run stopped after the gate.
            if strcmp(m, 'diffusion')
                featSupported.setFalse({'SchedStrategy_LCFS','SchedStrategy_LCFSPR'});
            end
            % DPS closes on the covariance BETWEEN the class coordinates of a
            % station, not on the station total. `minnormal` carries those
            % blocks through its outer iteration; the DAE form has no unknown
            % for them, since adding a matrix block per station restores the
            % quartic cost that keeping Sigma out of the Newton vector avoids.
            % GPS is already excluded above, for `minnormal` only.
            % A stochastic Petri net has no drift outside the DAE form: its
            % conserved quantities are P-invariants rather than chain
            % populations, and an immediate transition is an algebraic FLOW
            % rather than an event with a rate. Every other fluid method builds
            % its drift from the station/class/phase encoding, where a Place
            % contributes no coordinate at all, so it would integrate the net as
            % an empty model and report zeros without a warning.
            if ~strcmp(m, 'dae')
                featSupported.setFalse({'Place','Transition','Enabling', ...
                    'Inhibiting','Timing','Firing','Storage','Linkage'});
            end
            if strcmp(m, 'dae')
                featSupported.setFalse({'SchedStrategy_DPS'});
                % A finite capacity region is a linear inequality on the state,
                % which the DAE form can carry as an algebraic equation beside
                % the drift and no ODE method can carry at all. The gate is
                % where this has to be declared: RUNANALYZER's own refusal sits
                % downstream of NETWORKSOLVER.RUNANALYZERCHECKS, so without this
                % the model is rejected as an unsupported feature before the
                % method is ever consulted. FLUID_REGION_CONSTRAINTS still
                % refuses the region forms that are not constraints on this
                % drift, by name.
                featSupported.setTrue({'Region'});
            end
            % GPS divides the server by weight among the BACKLOGGED classes, so
            % its share is a function of the backlog INDICATOR. A first-order
            % closure cannot express it at all: with continuous x_k > 0 every
            % class is always backlogged and the share collapses to the constant
            % w_k/sum_j w_j, the heavy-traffic limit, regardless of load. Only
            % `minnormal` supplies the P(X_k >= 1) that the closure needs.
            % `refined` is excluded too: it reads its rates at the MEAN-FIELD
            % variance, which is exactly the degenerate case. matrix/pnorm build
            % a PS drift, and GPS is not PS -- equal-weight GPS gives each
            % backlogged CLASS an equal share, PS each JOB.
            if ~strcmp(m, 'minnormal')
                featSupported.setFalse({'SchedStrategy_GPS'});
            end
            % HOL allocates capacity in PRIORITY order, not in proportion to
            % population. No fluid drift anywhere reads sn.classprio except
            % SOLVER_MFQ_PRIO on the single-queue mfq path, so every other
            % method would answer a priority model as if the classes shared the
            % server proportionally.
            if ~strcmp(m, 'mfq')
                featSupported.setFalse({'SchedStrategy_HOL'});
            end
            % MULTISERVER (registry name since 2026-09-05): the drifts carry
            % min(n,c) except SOLVER_FLUID_DIFFUSION, whose SDE is written for
            % one or infinitely many servers, and SOLVER_MFQ, a single-queue
            % model on the same server counts (off them 'mfq' resolves to
            % 'matrix' above, so the delta binds only where it runs as itself).
            % FLUID_METHOD_REFUSAL keeps wording the diffusion refusal.
            if any(strcmp(m, {'diffusion','mfq'}))
                featSupported.setFalse({'MultiServer'});
            end
            % FINITECAPACITY IS DECLARED IN THE BASE ENVELOPE AND WITHDRAWN FROM
            % NO METHOD, on purpose: which methods carry a buffer is decided by
            % the structural predicate FLUID_METHOD_REFUSAL already asks
            % (NetworkSolver.checkBindingCapacity, over MNetwork.findBindingCapacity
            % -- the same answer the recorder marks the name on), and RUNANALYZER
            % stops on that same call. A per-method delta here would be a SECOND
            % rule for one fact, and the worse of the two would win: it runs
            % first, so it replaced "Finite station capacity (setCapacity=2) at
            % station 'Q1' ... Use options.method='dae', which carries the buffer
            % as an algebraic constraint on the drift" -- which names the
            % station, the value and the fix -- with a bare "(feature:
            % FiniteCapacity)", and testsFLD/test_solver_fld_dae pins the former.
            % It also refused the age-of-information arm of 'mfq' (a bufferless
            % or single-buffer queue by definition, AOI_IS_AOI) through the
            % COARSE gate in RUNANALYZER, which exempts 'dae' and the qsys names
            % but not 'mfq'. The three arms that do serve a buffer are 'dae' (an
            % algebraic constraint on the drift), 'mol' (the Mt/G/s/0 loss
            % system, where the server count IS the buffer) and that AoI arm;
            % 'default' is gated as its resolution, which is 'dae' on a binding
            % buffer.
        end

        function [bool, reason] = supportsModelMethod(self, method)
            % The rules the feature registry cannot name, stated here so that
            % a CALLER sees them before running: the single-queue shape of
            % 'mfq', the single-station shape and horizon of the fluid limits,
            % a Cache node for 'rmf', the moment-closure applicability of
            % 'minnormal', 'refined' and 'dae', the server count of
            % 'diffusion', the fork-join route, and the binding buffer nothing
            % in the FLD tree but 'dae' reads. Every one of them is a
            % predicate FLUID_METHOD_REFUSAL delegates to and RUNANALYZER stops
            % on, so the report and the run cannot answer differently.
            %
            % Left only in the analyzers those rules were invisible to every
            % gate above them: findSolver offered all 31 fluid methods on the
            % BAS-blocking model of cqn_bas_blocking, 'fluid.mfq' on a tandem
            % and 'fluid.mol' on a Source -> Delay -> Sink, and each run then
            % stopped, or answered as the matrix method under the caller's
            % label.
            %
            % A NAME THAT RESOLVES IS ASKED THROUGH ITS RESOLUTION, not as a
            % name of its own: the featset and every structural rule are those
            % of the method that will run (see RESOLVEMETHOD), so a Petri net,
            % a capped model or a GPS model is gated as 'dae' or 'minnormal',
            % which is what SolverFLD(model) goes on to run, and 'mfq' on a
            % tandem or 'rmf' with no cache is gated as 'matrix', which is what
            % answers. SUPPORTSRESOLVEDMETHOD resolves before calling here; a
            % direct call resolves the same way, so the two cannot disagree.
            m = SolverFLD.canonicalMethod(method);
            isNet = isa(self.model, 'Network');
            if isNet && any(strcmp(m, {'default','mfq','rmf'}))
                opts = self.getOptions();
                opts.method = m;
                m = SolverFLD.canonicalMethod(self.resolveMethod(opts));
            end
            [bool, reason] = supportsModelMethod@NetworkSolver(self, m);
            if ~bool || ~isNet
                return
            end
            opts = self.getOptions();
            [bool, reason] = fluid_method_refusal(self.gateStruct(opts), m, opts, self.model);
        end

        function method = resolveMethod(self, options)
            % METHOD = RESOLVEMETHOD(OPTIONS)
            % The concrete method OPTIONS.METHOD runs as on this model, so
            % every gate above the dispatch (runAnalyzerChecks through
            % supportsResolvedMethod, needsMapEnv, findSolver) asks about the
            % method that answers:
            %   'default'   FLUID_RESOLVE_DEFAULT_METHOD, exactly as RUNANALYZER
            %               resolves it; lang='matlab' only, since the java
            %               path hands 'default' to JLINE.SolverFluid, which
            %               resolves it itself
            %   'mfq' and its aliases   'matrix' off the single-queue shape
            %               FLUID_MFQ_ADMITS decides. That is the documented
            %               fallback the analyzer arm keeps (it warns and
            %               re-enters the matrix method), stated here so the
            %               pair is gated and labelled as the method that
            %               answers rather than refused, or offered as 'mfq'
            %   'rmf'       'matrix' with no Cache node to decompose: the
            %               decomposition's network step IS the matrix method
            method = options.method;
            if ~isa(self.model, 'Network')
                return
            end
            switch SolverFLD.canonicalMethod(method)
                case 'default'
                    if isfield(options, 'lang') && ~strcmp(options.lang, 'matlab')
                        return
                    end
                    method = fluid_resolve_default_method(self.gateStruct(options), options, self.model);
                case 'mfq'
                    if ~fluid_mfq_admits(self.gateStruct(options))
                        method = 'matrix';
                    end
                case 'rmf'
                    if ~any(self.model.getStruct().nodetype == NodeType.Cache)
                        method = 'matrix';
                    end
            end
        end

        function sn = gateStruct(self, options)
            % SN = GATESTRUCT(OPTIONS)
            % The struct the gate reasons on, prepared as RUNANALYZER prepares
            % the one it solves: non-Markovian laws already fitted to
            % phase-type, so a phase count or an arrival-process test answers
            % about the model the run sees rather than about the one the user
            % typed. SSA and the fluid ODEs read mu*phi as a flow, so the
            % surrogate must be a genuine phase-type ('ph'), not a matrix
            % exponential.
            if nargin < 2
                options = self.getOptions();
            end
            sn = self.model.getStruct();
            options.config.phfit = 'ph';
            sn = sn_nonmarkov_toph(sn, options);
        end


    end
    
    methods (Static)

        function m = canonicalMethod(method)
            % M = CANONICALMETHOD(METHOD)
            % The one spelling of a fluid method that every gate tests against.
            %
            % Three families of alias reach this solver and they used to be
            % expanded by hand at each branch, which is why the four codebases
            % drifted apart: a 'fluid.' qualifier RUNANALYZER accepts on every
            % name, the MFQ backend aliases 'butools' and 'aoi', and the short
            % spellings 'ggisgi' and 'tga' of the two single-station limits.
            % Canonicalizing once is what makes an alias carry the same feature
            % envelope as the name it resolves to; the JAR, native python and
            % C++ apply the same three rules in the same order.
            m = method;
            if ~ischar(m) && ~isstring(m)
                return
            end
            m = char(m);
            if strncmp(m, 'fluid.', 6)
                m = m(7:end);
            end
            % The MFQ backend and its age-of-information reading are the same
            % drift as 'mfq': RUNANALYZER dispatches all three to SOLVER_MFQ.
            if any(strcmp(m, {'butools','aoi'}))
                m = 'mfq';
                return
            end
            % The short spellings of the two single-station limits, as
            % fluid_qsys_canonical maps them.
            if strcmp(m, 'ggisgi')
                m = 'ggisgi.fluid';
            elseif strcmp(m, 'tga')
                m = 'ggingi.tga';
            end
        end

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()
            
            featSupported = SolverFeatureSet;
            featSupported.setTrue({
                'ClassSwitch','Delay','DelayStation','Queue',...
                'Cache','CacheClassSwitcher',...
                ... % 'CacheRetrieval' is deliberately NOT declared: no fluid
                ... % code anywhere implements delayed-hit retrieval, and on
                ... % examples/basic/cacheModel/retrieval_simple the ODE
                ... % returned zero QLen, Util and Tput on every row while jobs
                ... % arrived at rate 1, i.e. flow was not conserved. Refusing
                ... % the model is the honest answer; the JAR SolverFluid does
                ... % the same.
                'Cox2','Coxian','Erlang','Exp','HyperExp',...
                'APH', 'Det','MAP','MMPP2','NHPP','MAPt','PHt',...
                ... % Non-Markovian renewal distributions: converted to acyclic PH
                ... % by sn_nonmarkov_toph in runAnalyzer, so the fluid ODE solves them.
                'Gamma','Lognormal','Pareto','Uniform','Weibull',...
                'StatelessClassSwitcher','InfiniteServer','SharedServer','Buffer','Dispatcher',...
                'Server','ServiceTunnel',...
                'LoadDependence',... % closing family only, see getMethodFeatureSet
                'SchedStrategy_INF','SchedStrategy_PS',...
                'SchedStrategy_DPS','SchedStrategy_GPS',... % GPS: minnormal only, see getMethodFeatureSet
                'SchedStrategy_FCFS','SchedStrategy_HOL',...
                'SchedStrategy_SIRO','SchedStrategy_LCFS','SchedStrategy_LCFSPR',...
                'ReplacementStrategy_RR','ReplacementStrategy_FIFO',... % RAND(m)/FIFO(m) refined mean field
                'ReplacementStrategy_SFIFO',... % strict FIFO(m) position-resolved mean field; LRU/HLRU/CLIMB/QLRU rejected at runtime
                ...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'Fork','Join','Forker','Joiner',... % fork-join via the MMT transformation (fjFixedPoint), see fldDispatch
                'JoinPartial',... % quorum join: the MMT fixed point charges the k-th branch completion (fj_ordstat_exp)
                'ClosedClass','SelfLoopingClass','Replayer', ...
                ... % Stochastic Petri nets: the 'dae' method only, see
                ... % getMethodFeatureSet. A Transition node routes the model to
                ... % SOLVER_FLUID_PETRI, which solves the marking as the same
                ... % min-normal closure with the P-invariants as constraints and
                ... % the immediate firing flows as algebraic unknowns.
                ... % 'Storage'/'Linkage' ride along with any Place, as they do
                ... % in the SSA and CTMC sets: GETUSEDLANGFEATURES emits all
                ... % three for one node, so declaring Place without them refuses
                ... % every Petri net at the gate.
                ... % 'Inhibiting' is declared, but an inhibitor arc on an
                ... % IMMEDIATE mode is answered wrongly when the inhibitor place's
                ... % mean sits at its threshold -- the mode is switched off and its
                ... % branch carries no flow. Closed for a TIMED mode only; see
                ... % _kb/06-solver-catalog.md and FLUID_PETRI_IMMEDIATE.
                'Place','Transition','Enabling','Inhibiting','Timing','Firing',...
                'Storage','Linkage',...
                'RandomSource','Sink','Source','OpenClass','JobSink',...
                ... % c-server stations: the drifts carry min(n,c); withdrawn
                ... % from 'diffusion' and 'mfq' in getMethodFeatureSet.
                'MultiServer',...
                ... % A binding buffer: 'dae' carries it as an algebraic
                ... % constraint, 'mol' IS the Mt/G/s/0 loss system and the AoI
                ... % arm of 'mfq' is a bufferless or single-buffer queue. WHICH
                ... % method serves one is the structural rule FLUID_METHOD_REFUSAL
                ... % asks and RUNANALYZER stops on, so no per-method delta
                ... % duplicates it here; see getMethodFeatureSet.
                'FiniteCapacity'});
            %SolverFLD has very weak performance on open models
        end
        
        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverFLD.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end        
        
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = SolverOptions('Fluid');
        end

        function libs = getLibrariesUsed(sn, options)
            % GETLIBRARIESUSED Get list of external libraries used by FLD solver
            libs = {};
            if nargin >= 2 && isfield(struct(options), 'method')
                if any(strcmp(options.method, {'rmf', 'fluid.rmf'}))
                    libs{end+1} = 'rmf_tool';
                end
            end
        end
    end
end
