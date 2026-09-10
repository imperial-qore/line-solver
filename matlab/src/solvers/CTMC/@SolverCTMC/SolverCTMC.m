classdef SolverCTMC < NetworkSolver
    % Continuous-Time Markov Chain solver for exact state-space analysis
    %
    % Implements exact analysis of queueing networks via CTMC formulation.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Access = public)
        % User-supplied MarkovProcess (CTMC) or MarkovChain (DTMC) when the
        % solver runs in chain mode, empty when it analyzes a Network.
        chainModel = [];
    end

    methods
        function self = SolverCTMC(model,varargin)
            % SOLVERCTMC Create a CTMC solver instance
            %
            % @brief Creates a CTMC solver for exact Markov chain analysis
            % @param model Network model, or a MarkovProcess/MarkovChain to be
            %              solved directly (chain mode)
            % @param varargin Optional parameters (cutoff, method, etc.)
            % @return self SolverCTMC instance configured for exact analysis

            % Redirect a LayeredNetwork to the solvers that accept one. This
            % must precede the NetworkSolver constructor: that constructor calls
            % model.getAvgHandles(), which a LayeredNetwork does not implement,
            % so the redirect that used to live in runAnalyzer was unreachable
            % and the user saw "Unrecognized method 'getAvgHandles'" instead.
            if isa(model, 'LayeredNetwork')
                line_error(mfilename, ['This model is a LayeredNetwork. Use SolverLN ', ...
                    '(iterative) or SolverLQNS (external LQNS) instead of SolverCTMC.\n']);
            end
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            if isa(model,'MarkovProcess') || isa(model,'MarkovChain')
                % Chain mode: the state space and the generator are given, so
                % there is nothing to translate to another language backend.
                self.chainModel = model;
                return
            end
            self.setLang();
        end

        function bool = isChainSolver(self)
            % BOOL = ISCHAINSOLVER()
            % True when the solver was constructed from a MarkovProcess or a
            % MarkovChain instead of a Network.
            bool = ~isempty(self.chainModel);
        end

        function bool = isDiscreteChain(self)
            % BOOL = ISDISCRETECHAIN()
            % True in chain mode when the user supplied a DTMC (MarkovChain).
            bool = isa(self.chainModel,'MarkovChain');
        end

        function P = getTransMat(self)
            % P = GETTRANSMAT()
            % Transition matrix of the user-supplied DTMC (chain mode only).
            if ~self.isDiscreteChain()
                line_error(mfilename,'getTransMat requires a SolverCTMC built from a MarkovChain.');
            end
            P = self.chainModel.getTransMat();
        end

        function assertChainMode(self, callerName, wantChain)
            % ASSERTCHAINMODE(CALLERNAME, WANTCHAIN)
            % Guard that keeps the network-oriented and the chain-oriented
            % method sets from being mixed on the same solver instance.
            if wantChain && ~self.isChainSolver()
                line_error(mfilename,'%s requires a SolverCTMC built from a MarkovProcess or a MarkovChain.', callerName);
            elseif ~wantChain && self.isChainSolver()
                line_error(mfilename,['%s requires a Network model. This solver was built from a %s, ', ...
                    'which has no stations or classes: use getProbSys, getGenerator, getStateSpace, ', ...
                    'getTranProbSys or sampleSys instead.'], callerName, class(self.chainModel));
            end
        end


        runtime = runAnalyzer(self, options)
        Pnir = getProb(self, node, state)
        Pn = getProbSys(self)
        Pnir = getProbAggr(self, ist)
        
        Pn = getProbSysAggr(self)
        [Pi_t, SSsysa] = getTranProbSysAggr(self)   
        [Pi_t, SSnode_a] = getTranProbAggr(self, node)
        [Pi_t, SSsys] = getTranProbSys(self)     
        [Pi_t, SSnode] = getTranProb(self, node)
        RD = getCdfRespT(self, R)
        RD = getCdfSysRespT(self)
        [RD, out] = getCdfFirstPassT(self, A, B)
        [m, mall] = getFirstPassTMoments(self, A, B, nmax)
        
        [stateSpace,nodeStateSpace] = getStateSpace(self, options)
        [sn, options, isFJ] = fjAugment(self, sn, options)
        [QN,UN,RN,TN,CN,XN] = ctmcFesAggregation(self, sn, options)
        stateSpaceAggr = getStateSpaceAggr(self)

        % Reward computation methods
        [Rt, t, names] = getTranReward(self, rewardName)
        [R, names] = getAvgReward(self)
        [V, t, names, stateSpace] = runRewardAnalyzer(self)
        
        function [state_space, local_states] = stateSpace(self)
            % STATESPACE Alias for getStateSpace
            if nargout <= 1
                state_space = self.getStateSpace();
            else
                [state_space, local_states] = self.getStateSpace();
            end
        end
        
        function Q = generator(self)
            % GENERATOR Alias for getGenerator
            Q = self.getGenerator();
        end
            
        [infGen, eventFilt, synchInfo, stateSpace, nodeStateSpace] = getSymbolicGenerator(self, invertSymbol)
        [infGen, eventFilt, synchInfo] = getInfGen(self, options)        
        [infGen, eventFilt, synchInfo] = getGenerator(self, options)
        
        tstate = sampleSys(self, numevents)
        sampleAggr = sampleAggr(self, node, numEvents)
        
        function MCTMC = getMarkedCTMC(self, options)        
            % MCTMC = GETMARKEDCTMC(options)

            if nargin < 2
                [infGen, eventFilt, synchInfo] = self.getInfGen();    
            else
                [infGen, eventFilt, synchInfo] = getInfGen(self, options);    
            end
            
            MCTMC = MarkedMarkovProcess(infGen, eventFilt, synchInfo);
        end

        function sn = getStruct(self)
            % QN = GETSTRUCT()

            % Get data structure summarizing the model
            self.assertChainMode('getStruct', false);
            sn = self.model.getStruct(true);
        end
        

        function bool = supportsTransientAnalysis(self) %#ok<MANU>
            % Transient averages are available (uniformization of the generator over options.timespan).
            bool = true;
        end

        function reason = unsupportedMethodReason(self, method) %#ok<INUSL>
            % REASON = UNSUPPORTEDMETHODREASON(METHOD)
            %
            % The forwarding address for the QRF reduction bounds, which are
            % SolverBA's now. runAnalyzer carried this text, which is DOWNSTREAM
            % of runAnalyzerChecks and so never ran: listValidMethods no longer
            % names qrf.*, and the flat "unsupported by this solver" fired
            % first. Asking nothing of the model is what lets the name gate call
            % it.
            reason = '';
            if ischar(method) && startsWith(method, 'qrf')
                reason = sprintf(['The QRF reduction bounds (method ''%s'') moved to ' ...
                    'SolverBA. Use SolverBA(model,''method'',''%s'') (or the ''qr''/''lr'' ' ...
                    'aliases).'], method, method);
            end
        end

        function [allMethods] = listValidMethods(self)
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            if self.isChainSolver()
                % Chain mode solves the given generator directly, there is no
                % state-space generation method to choose.
                allMethods = {'default'};
                return
            end
            sn = self.model.getStruct();
            % QRF (Quadratic/Linear Reduction) LP-based bounds were moved to
            % SolverBA (methods qr/lr/qrf.*); SolverCTMC no longer serves them.
            % 'exact' is an explicit alias for the default state-space path: it
            % pins the intent at the call site so an example or test cannot be
            % re-baselined by a later change of what 'default' selects. It must
            % stay behaviourally identical to 'default' -- neither
            % solver_ctmc_analyzer nor ctmc_solve branches on the name (in
            % ctmc_solve it is neither 'gmres', 'direct' nor 'gpu', so it takes
            % the same path), and that equivalence is the point of the alias.
            %
            % 'gpu' IS THE SAME ALIAS IN MATLAB. No file under matlab/src reads
            % the name (ctmc_stationary has no gpuArray path), so it runs the
            % CPU solve; it is kept because the reference lists carry it and
            % the JAR selects a backend on it.
            %
            % 'mdd' holds the reachable set in a decision diagram and solves
            % K coupled level-CTMCs instead of the |S|-state generator; it is
            % exact on product-form models and approximate otherwise, and is
            % restricted to closed single-class networks (solver_ctmc_mdd_analyzer).
            allMethods = {'default','exact','gpu','mdd','cftp','cftp.approx'};
        end
        function featSupported = getMethodFeatureSet(self, method)
            % Per-method feature deltas applied to the base CTMC envelope.
            % Four of the six methods share it; 'cftp'/'cftp.approx' and 'mdd'
            % narrow it, because neither builds the explicit generator that
            % carries the rest of the envelope.
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
            featSupported = SolverCTMC.getFeatureSet();
            if any(strcmp(method, {'cftp','cftp.approx'}))
                % PERFECT SAMPLING FROM A BALANCE FUNCTION, not from a
                % generator: the sampler encodes the closed single-class
                % product form of Gordon-Newell and nothing else, so every
                % construct outside it has to leave the envelope. The class
                % count and the station count have no registry name and are
                % checked structurally in supportsModelMethod, against the
                % same predicate the analyzer uses.
                featSupported.setFalse({'OpenClass', ...
                    ... % Queue, Delay and Router are the only node kinds the sampler walks
                    'Source','Sink','RandomSource','JobSink', ...
                    'ClassSwitch','StatelessClassSwitcher', ...
                    'Cache','CacheClassSwitcher','CacheRetrieval', ...
                    'ReplacementStrategy_RR','ReplacementStrategy_FIFO','ReplacementStrategy_SFIFO', ...
                    'ReplacementStrategy_LRU','ReplacementStrategy_HLRU', ...
                    'ReplacementStrategy_CLIMB','ReplacementStrategy_QLRU', ...
                    'Fork','Join','Forker','Joiner','ForkFanoutVector', ...
                    'Place','Transition','Linkage','Enabling','Inhibiting','Timing','Firing','Storage', ...
                    ... % disciplines outside INF/PS/FCFS/SIRO/LCFSPR have no product form
                    'SchedStrategy_DPS','SchedStrategy_GPS', ...
                    'SchedStrategy_SEPT','SchedStrategy_LEPT', ...
                    'SchedStrategy_HOL','SchedStrategy_LCFS', ...
                    'SchedStrategy_LCFSPRPRIO','SchedStrategy_FCFSPRPRIO', ...
                    'SchedStrategy_FCFSPR','SchedStrategy_LCFSPI','SchedStrategy_FCFSPI', ...
                    'SchedStrategy_LCFSPIPRIO','SchedStrategy_FCFSPIPRIO', ...
                    'SchedStrategy_PSPRIO','SchedStrategy_DPSPRIO','SchedStrategy_GPSPRIO', ...
                    'SchedStrategy_LPS','SchedStrategy_PAS','SchedStrategy_OI','SchedStrategy_POLLING', ...
                    ... % The one-phase-per-station rule is deliberately NOT
                    ... % spelled as a list of distribution names. The sampler
                    ... % refuses sn.phases(i,1) > 1, and a name is not a phase
                    ... % count: a one-phase Coxian passes and a HyperExp does
                    ... % not, while Det/Gamma/Pareto only acquire their phases
                    ... % in SN_NONMARKOV_TOPH. supportsModelMethod asks the
                    ... % phase count instead, which is also what lets it name
                    ... % the offending station.
                    'Region', ...
                    'LoadDependence','ClassDependence','JointDependence','GlobalDependence', ...
                    ... % a state-dependent decision is not Markovian routing
                    'RoutingStrategy_RROBIN','RoutingStrategy_WRROBIN', ...
                    'RoutingStrategy_JSQ','RoutingStrategy_SQ','RoutingStrategy_SDR', ...
                    ... % the balance function has no orbit, no abandonment, no
                    ... % outage and no trigger class; each is a State construct
                    'Balking','Reneging','Retrial','Breakdown', ...
                    'OpenSignal','ClosedSignal','SignalType_NEGATIVE','SignalType_REPLY', ...
                    'SignalType_CATASTROPHE','SignalBatchRemoval','SignalRemovalPolicy', ...
                    ... % the Gordon-Newell balance function has no buffer:
                    ... % solver_ctmc_cftp_supports refuses a finite one by name
                    'FiniteCapacity'});
            elseif strcmp(method, 'mdd')
                % The decision diagram holds the MARKING of a closed network;
                % an open stream makes it unbounded, so there is no finite
                % diagram to hold. The single-class rule is structural (no
                % registry name for a class count) and lives in
                % supportsModelMethod. A stochastic Petri net keeps the
                % Place/Transition names: SPN_MDD reads the marking directly.
                %
                % A FORK-JOIN MODEL IS NEITHER of the two shapes it serves. The
                % tag augmentation a fork needs adds one auxiliary class per
                % branch, so the struct that reaches the analyzer is never
                % single-class however the model was written, and the level
                % decomposition has no meaning for a firing that does not
                % conserve the per-chain population.
                %
                % THE DESCRIPTOR READS RATES, SERVERS, PHASES, DISCIPLINES AND
                % THE ROUTING CHAIN, AND NOTHING ELSE. Every construct the
                % explicit generator carries through State.afterEventStation
                % (impatience, an orbit, an outage, a cache, a signal, a
                % scaling handle, a region, a polling table, a service-rate
                % function) would be dropped on the floor, so each leaves the
                % envelope here. A Router is stateful and breaks the
                % station-to-station chain the predicate needs stochastic.
                featSupported.setFalse({'OpenClass', ...
                    'Source','Sink','RandomSource','JobSink', ...
                    'Fork','Join','Forker','Joiner','JoinPartial','ForkFanoutVector', ...
                    'Router', ...
                    'Balking','Reneging','Retrial','Breakdown', ...
                    'Cache','CacheClassSwitcher','CacheRetrieval', ...
                    'ReplacementStrategy_RR','ReplacementStrategy_FIFO','ReplacementStrategy_SFIFO', ...
                    'ReplacementStrategy_LRU','ReplacementStrategy_HLRU', ...
                    'ReplacementStrategy_CLIMB','ReplacementStrategy_QLRU', ...
                    'Region', ...
                    'LoadDependence','ClassDependence','JointDependence','GlobalDependence', ...
                    'OpenSignal','ClosedSignal','SignalType_NEGATIVE','SignalType_REPLY', ...
                    'SignalType_CATASTROPHE','SignalBatchRemoval','SignalRemovalPolicy', ...
                    'SchedStrategy_POLLING','SchedStrategy_PAS','SchedStrategy_OI', ...
                    ... % the level decomposition reads rates, servers and phases
                    ... % and no sn.cap/classcap, so a buffer would be dropped
                    'FiniteCapacity'});
            end
        end

        function [bool, reason] = supportsModelMethod(self, method)
            % The per-method rules the feature registry has no name for, asked
            % of the SAME predicates the analyzers use so that the report and
            % the run cannot answer differently.
            %
            % Three of them: the class count and the station count that 'cftp'
            % and 'mdd' need (a class count is not a model feature), and the
            % state-space size that the explicit-generator methods need. The
            % last one is why 'default'/'exact'/'gpu' were offered on models
            % whose chain does not fit memory -- the analyzer priced the state
            % space and refused, and nothing above it had asked.
            %
            % THE TWO STRUCTURAL PREDICATES ARE ASKED BEFORE THE FEATURE GATE,
            % which is the reverse of the usual order and deliberate: each is
            % the analyzer's own assert, so it refuses a strict superset of
            % what the per-method feature deltas above refuse, and its wording
            % names the offending station or class count instead of a feature.
            % Asking the feature gate first would replace 'the cftp method
            % supports closed models only' with '(feature: OpenClass)' on the
            % very run the caller is about to make.
            if isa(self.model, 'Network')
                if any(strcmp(method, {'cftp','cftp.approx'}))
                    [bool, reason] = solver_ctmc_cftp_supports(self.model.getStruct(), self.getOptions());
                    if ~bool
                        return
                    end
                elseif strcmp(method, 'mdd')
                    [bool, reason] = solver_ctmc_mdd_supports(self.model.getStruct(), self.getOptions());
                    if ~bool
                        return
                    end
                end
                % A finite horizon on a fork-join model, which RUNANALYZER
                % refuses after the gate has answered; asked here of the same
                % predicate so the report does not offer what the run refuses.
                [bool, reason] = SolverCTMC.transientSupports(self.model.getStruct(), self.getOptions());
                if ~bool
                    return
                end
                % The fork-join model class, which EVERY method has to clear:
                % the tag augmentation runs before the state space, the decision
                % diagram and the sampler alike, so a model SN_FJ_VALIDATE
                % refuses is refused whichever name was asked for. The featset
                % cannot state it -- Fork and Join are declared and the rules
                % are about how they are wired -- so it is structural, and it is
                % the validator's own body of rules rather than a copy.
                [bool, reason] = sn_fj_supports(self.model.getStruct());
                if ~bool
                    return
                end
            end
            [bool, reason] = supportsModelMethod@NetworkSolver(self, method);
            if ~bool || ~isa(self.model, 'Network')
                return
            end
            if ~any(strcmp(method, {'cftp','cftp.approx','mdd'}))
                % The impatience laws and the server counts the State
                % machinery serves: the rules solver_ctmc.m raises, asked of
                % the same predicate (SOLVER_CTMC_STATE_SUPPORTS) before the
                % run so a phase-type patience is refused here and not there.
                [bool, reason] = solver_ctmc_state_supports(self.model.getStruct(), 'SolverCTMC');
                if ~bool
                    return
                end
                % The explicit state space is what these methods enumerate, and
                % ctmc_memory_gate refuses it above the host budget. Asking the
                % same estimator here costs a combinatorial formula, not a
                % state space, so the report stays cheap.
                [bool, reason] = SolverCTMC.isStateSpaceTractable(self.model, self.getOptions());
            end
        end


    end
           
    methods (Static)

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()

            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Source','Sink',...
                'ClassSwitch','Delay','DelayStation','Queue','Router',...
                'MAP','APH','MMPP2','MMAP','PH','Coxian','Erlang','Exp','HyperExp','ME',...
                'Det','Gamma','Lognormal','Pareto','Uniform','Weibull',...
                'StatelessClassSwitcher','InfiniteServer','SharedServer','Buffer','Dispatcher',...
                ... % Finite capacity regions: CTMC represents them exactly (the
                ... % blocked-job overflow buffer is part of the chain), verified by
                ... % convergence to the simulators as the cutoff grows.
                'Region', ...
                'Cache','CacheClassSwitcher', ...
                'CacheRetrieval', ...
                'Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS',...
                'SchedStrategy_DPS','SchedStrategy_GPS',...
                'SchedStrategy_SIRO','SchedStrategy_SEPT',...
                'SchedStrategy_LEPT','SchedStrategy_FCFS',...
                'SchedStrategy_HOL','SchedStrategy_LCFS',...
                ... % LCFS with priority groups: buffered like HOL in State.fromMarginal and
                ... % served by its own arm of afterEventStation, so it was an under-claim.
                'SchedStrategy_LCFSPRIO',...
                'SchedStrategy_LCFSPR','SchedStrategy_LCFSPRPRIO','SchedStrategy_FCFSPRPRIO',...
                ... % Rest of the preempt family, all handled in afterEventStation.m (restored after a merge dropped it).
                'SchedStrategy_LCFSPI','SchedStrategy_FCFSPR','SchedStrategy_FCFSPI',...
                'SchedStrategy_LCFSPIPRIO','SchedStrategy_FCFSPIPRIO',...
                'SchedStrategy_PSPRIO','SchedStrategy_DPSPRIO','SchedStrategy_GPSPRIO',...
                'SchedStrategy_LPS',...
                'SchedStrategy_PAS',...
                'SchedStrategy_OI',...
                'SchedStrategy_POLLING',...
                'RoutingStrategy_RROBIN',...
                'RoutingStrategy_WRROBIN',...
                'RoutingStrategy_JSQ',...
                'RoutingStrategy_SQ',...
                'RoutingStrategy_SDR',...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO','ReplacementStrategy_SFIFO','ReplacementStrategy_LRU',...
                'ReplacementStrategy_HLRU','ReplacementStrategy_CLIMB','ReplacementStrategy_QLRU',...
                'ClosedClass','SelfLoopingClass','OpenClass','Replayer',...
                'OpenSignal','ClosedSignal',...
                'SignalType_NEGATIVE','SignalType_CATASTROPHE','SignalType_REPLY',...
                'SignalBatchRemoval','SignalRemovalPolicy',...
                'Place','Transition','Linkage','Enabling','Inhibiting','Timing','Firing','Storage',...
                'Fork','Join','Forker','Joiner',...
                ... % A per-destination tasks-per-link vector is carried as one integer
                ... % weight per branch by the tag construction (fjtag, State.afterFJEvent);
                ... % a random count or a branch probability is not (sn_fj_validate).
                'ForkFanoutVector',...
                'Balking','Reneging','Retrial','Breakdown',...
                'LoadDependence','ClassDependence','JointDependence','GlobalDependence',...
                ... % c-server stations and binding buffers are both State
                ... % constructs (State.fromMarginal / afterEventStation): served
                ... % by the explicit generator, withdrawn from cftp and mdd
                'MultiServer','FiniteCapacity'});
        end
        
        function [bool, featSupported, featUsed] = supports(model)
            % [BOOL, FEATSUPPORTED, FEATUSED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverCTMC.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end        
        
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = SolverOptions('CTMC');
        end

        function [bool, reason] = transientSupports(sn, options)
            % [BOOL, REASON] = TRANSIENTSUPPORTS(SN, OPTIONS)
            %
            % Can a finite options.timespan be asked of this model? The tag
            % augmentation of a fork-join model has no transient counterpart
            % (SN_FJ_FOLDBACK folds stationary means only), so RUNANALYZER
            % refuses the pair. One body, two callers: supportsModelMethod
            % answers with it, runAnalyzer raises with it.
            bool = true;
            reason = '';
            if nargin < 2 || isempty(options) || ~isfield(options,'timespan') ...
                    || isempty(options.timespan) || isinf(options.timespan(1))
                return
            end
            if any(sn.nodetype == NodeType.Fork) || any(sn.nodetype == NodeType.Join)
                bool = false;
                reason = 'Transient analysis of fork-join models is not supported by SolverCTMC.';
            end
        end

        function [ok, msg, logNstates] = isStateSpaceTractable(model, options)
            % [OK, MSG, LOGNSTATES] = ISSTATESPACETRACTABLE(MODEL, OPTIONS)
            %
            % True when the worst-case CTMC state space of MODEL fits the
            % host memory budget. Same estimator and gate the analyzer runs,
            % exposed so that a caller (e.g. SolverAUTO) can rank CTMC out
            % before paying for state-space generation.

            if nargin < 2 || isempty(options)
                options = SolverCTMC.defaultOptions();
            end
            ok = true; msg = ''; logNstates = 0;
            if ~isa(model,'Network')
                return
            end
            try
                sn = sn_nonmarkov_toph(model.getStruct(), options);
                logNstates = ctmc_state_space_logsize(sn, options);
            catch ME
                % An estimator failure must not be read as a refusal: the
                % analyzer runs its own gate and reports the real error.
                msg = ME.message;
                return
            end
            forceFlag = isfield(options,'force') && ~isempty(options.force) && options.force;
            if isfield(options,'memorySafetyFraction') && ~isempty(options.memorySafetyFraction)
                safetyFraction = options.memorySafetyFraction;
            else
                safetyFraction = 0.6;
            end
            [ok, msg] = ctmc_memory_gate(logNstates, forceFlag, false, safetyFraction);
        end
        
        function printInfGen(Q,SS)
            % PRINTINFGEN(Q,SS)
            
            SS=full(SS);
            Q=full(Q);
            for s=1:size(SS,1)
                for sp=1:size(SS,1)
                    if Q(s,sp)>0
                        line_printf('\n%s->%s: %f',mat2str(SS(s,:)),mat2str(SS(sp,:)),double(Q(s,sp)));
                    end
                end
            end
            line_printf('\n');
        end
        
        function printEventFilt(sync,D,SS,myevents)
            % PRINTEVENTFILT(SYNC,D,SS,MYEVENTS)
            
            if nargin<4 %~exist('events','var')
                myevents = 1:length(sync);
            end
            SS=full(SS);
            for e=myevents
                D{e}=full(D{e});
                for s=1:size(SS,1)
                    for sp=1:size(SS,1)
                        if D{e}(s,sp)>0
                            line_printf('\n%s-- %d: (%d,%d) => (%d,%d) -->%s: %f',mat2str(SS(s,:)),e,sync{e}.active{1}.node,sync{e}.active{1}.class,sync{e}.passive{1}.node,sync{e}.passive{1}.class,mat2str(SS(sp,:)),double(D{e}(s,sp)));
                        end
                    end
                end
            end
        end

        function libs = getLibrariesUsed(sn, options) %#ok<INUSD>
            % GETLIBRARIESUSED Get list of external libraries used by CTMC solver
            libs = {};
        end

        function backend = symbolicBackend(self)
            % BACKEND = SYMBOLICBACKEND(SELF)
            %
            % Value of options.config.symbolic, or 'auto' when the solver
            % carries no such option (an older options struct, or none at all).
            backend = 'auto';
            try
                if isprop(self, 'options') && isfield(self.options, 'config') && ...
                        isfield(self.options.config, 'symbolic')
                    backend = self.options.config.symbolic;
                end
            catch
            end
        end

        function entries = symbolicEntries(infGenTerms, symbols, invertSymbol, n)
            % ENTRIES = SYMBOLICENTRIES(INFGENTERMS, SYMBOLS, INVERTSYMBOL, N)
            %
            % Renders the per-event coefficient matrices as expression strings,
            % e.g. '2*x1 - 3*x2'. This is the wire form the symbolic backend
            % consumes and the exact format the JAR's getSymbolicEntry emits,
            % which is what lets the two be compared at all.
            entries = cell(n, n);
            for i = 1:n
                for j = 1:n
                    s = '';
                    for e = 1:numel(symbols)
                        if isempty(symbols{e})
                            continue
                        end
                        c = infGenTerms{e}(i, j);
                        if c == 0
                            continue
                        end
                        if isempty(s)
                            if c < 0
                                s = '-';
                            end
                        else
                            if c < 0
                                s = [s, ' - ']; %#ok<AGROW>
                            else
                                s = [s, ' + ']; %#ok<AGROW>
                            end
                        end
                        absC = abs(c);
                        if invertSymbol
                            s = [s, SolverCTMC.symbolicCoeff(absC), '/', symbols{e}]; %#ok<AGROW>
                        else
                            if absC ~= 1
                                s = [s, SolverCTMC.symbolicCoeff(absC), '*']; %#ok<AGROW>
                            end
                            s = [s, symbols{e}]; %#ok<AGROW>
                        end
                    end
                    if isempty(s)
                        s = '0';
                    end
                    entries{i, j} = s;
                end
            end
        end

        function s = symbolicCoeff(c)
            % S = SYMBOLICCOEFF(C)  Coefficient as text the backend reads exactly.
            if c == round(c) && abs(c) < 1e15
                s = sprintf('%d', round(c));
            else
                s = sprintf('%.17g', c);
            end
        end
    end
end
