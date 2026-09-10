classdef SolverSSA < NetworkSolver
    % Stochastic Simulation Analysis solver
    %
    % Implements discrete-event stochastic simulation for queueing network analysis.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    methods
        function self = SolverSSA(model,varargin)
            % SOLVERSSA Create an SSA solver instance
            %
            % @brief Creates a Stochastic Simulation Analysis solver
            % @param model Network model to be analyzed via simulation
            % @param varargin Optional parameters (method, samples, seed, etc.)
            % @return self SolverSSA instance configured for simulation

            % An auxiliary solver passed as first optional argument requests a
            % warm start: its steady-state solution decides the initial state
            % of the simulated trajectory (see NetworkSolver.initFromSolver).
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
        
        function sn = getStruct(self)
            % QN = GETSTRUCT()
            
            % Get data structure summarizing the model
            sn = self.model.getStruct(true);
        end
        
        [runtime, tranSysState, tranSync] = run(self, options);        
        Prob = getProb(self, node, state);
        ProbAggr = getProbAggr(self, node, state);
        ProbSys = getProbSys(self);
        ProbSysAggr = getProbSysAggr(self);
        tranNodeState = sample(self, node, numEvents, markActivePassive);
        tranNodeStateAggr = sampleAggr(self, node, numEvents, markActivePassive);
        tranSysStateAggr = sampleSysAggr(self, numEvents, markActivePassive);
        tranSysState = sampleSys(self, numEvents, markActivePassive);

        function [allMethods] = listValidMethods(self)
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            
            %sn = self.model.getStruct();
            allMethods = {'default','ssa','serial', 'para','parallel','nrm'};
        end

        function bool = isStochasticMethod(self, method) %#ok<INUSD>
            % BOOL = ISSTOCHASTICMETHOD(METHOD)
            % All SSA methods are stochastic simulation.
            bool = true;
        end

        function featSupported = getMethodFeatureSet(self, method) %#ok<INUSD>
            % Every SSA method shares the solver envelope: the serial and the
            % NRM engines differ in what they PREFER, not in what the solver
            % answers, because the explicit 'nrm' arm falls back to the serial
            % engine for everything but a discipline without a reaction form,
            % and that one rule is structural (SSA_NRM_REFUSAL). Defining this
            % is what lets NetworkSolver.supportsModelMethod NAME the offending
            % features: with no method feature set it falls back to the coarse
            % supports(model), whose reason is empty.
            if ~isa(self.model, 'Network')
                featSupported = [];
                return;
            end
            featSupported = SolverSSA.getFeatureSet();
        end

        function [bool, reason] = supportsModelMethod(self, method)
            % The fork-join model class, which EVERY SSA method has to clear.
            %
            % RUNANALYZER tag-augments a fork-join model through
            % MODELADAPTER.FJTAG, whose first act is SN_FJ_VALIDATE, so a model
            % that validator refuses is refused whichever method was asked for.
            % The featset cannot state it -- Fork and Join are declared, and the
            % rules are about how they are WIRED (the pairing, the join
            % strategy, the tasks per link, whether an open class is routed
            % through the fork) -- so it is structural, and it is the
            % validator's own body of rules rather than a copy of them.
            %
            % Without it the report offered every ssa.* row on a fork-join model
            % whose Join names no fork, and each one then errored; SolverCTMC
            % gates on the same predicate for the same reason.
            [bool, reason] = supportsModelMethod@NetworkSolver(self, method);
            if bool && isa(self.model, 'Network')
                [bool, reason] = sn_fj_supports(self.model.getStruct());
            end
            % The impatience laws and the server counts the State machinery
            % serves, shared with SolverCTMC: solver_ssa.m raises them and the
            % NRM reads the same fields, so every SSA method is bound by them.
            if bool && isa(self.model, 'Network')
                [bool, reason] = solver_ctmc_state_supports(self.model.getStruct(), 'SolverSSA');
            end
            % A marking-dependent firing rate is refused by every SSA method
            % (SSA_FIRINGDEP_REFUSAL, which the analyzer raises with).
            if bool && isa(self.model, 'Network')
                [bool, reason] = ssa_firingdep_refusal(self.model.getStruct());
            end
            % The NRM engine is the one SSA method with a model class of its
            % own, and the test for it already existed: SOLVER_SSA_ANALYZER
            % consulted it to PREFER the NRM but nothing consulted it to decide
            % whether 'nrm' could be OFFERED, so the report listed ssa.nrm
            % Runnable on every model and an explicit request then raised
            % UnsupportedPolicy from SOLVER_SSA_ANALYZER_NRM. SSA_NRM_ELIGIBLE
            % is now that one predicate with two callers.
            if bool && isa(self.model, 'Network') && any(strcmpi(method, {'nrm','ssa.nrm'}))
                [bool, reason] = ssa_nrm_refusal(self.model.getStruct());
            end
        end
    end
    
    methods (Static)

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()

            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink','Source','Router',...
                'ClassSwitch','Delay','DelayStation','Queue',...
                'Cache','CacheClassSwitcher',...
                'CacheRetrieval', ...
                'MAP','MMPP2','MMAP', 'APH', 'PH', 'Replayer',...
                'Coxian','Erlang','Exp','HyperExp',...
                'Det','Gamma','Lognormal','Pareto','Uniform','Weibull',...
                'StatelessClassSwitcher','InfiniteServer',...
                'SharedServer','Buffer','Dispatcher',...
                ... % Finite capacity regions: the NRM sample path honours them.
                'Region', ...
                'Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS',...
                'SchedStrategy_DPS','SchedStrategy_FCFS',...
                'SchedStrategy_GPS','SchedStrategy_LPS','SchedStrategy_SIRO',...
                'SchedStrategy_HOL','SchedStrategy_LCFS',...
                'SchedStrategy_SEPT','SchedStrategy_LEPT',...
                'SchedStrategy_LCFSPR',...
                'SchedStrategy_PSPRIO','SchedStrategy_DPSPRIO','SchedStrategy_GPSPRIO',...
                'SchedStrategy_LCFSPRPRIO','SchedStrategy_FCFSPRPRIO',...
                ... % The rest of the preempt family and LCFS-with-priorities: the
                ... % serial engine drives the same State.afterEventStation arms
                ... % SolverCTMC declares; the NRM has no reaction form for them
                ... % and hands them to the serial engine (ssa_nrm_guards.sched).
                'SchedStrategy_LCFSPI','SchedStrategy_FCFSPR','SchedStrategy_FCFSPI',...
                'SchedStrategy_LCFSPIPRIO','SchedStrategy_FCFSPIPRIO','SchedStrategy_LCFSPRIO',...
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
                'SchedStrategy_EXT','ClosedClass','SelfLoopingClass','OpenClass',...
                'OpenSignal','ClosedSignal',...
                'SignalType_NEGATIVE','SignalType_CATASTROPHE',...
                'SignalBatchRemoval','SignalRemovalPolicy',...
                'Fork','Join','Forker','Joiner',...
                ... % A per-destination tasks-per-link vector is carried as one integer
                ... % weight per branch by the tag construction (fjtag, State.afterFJEvent);
                ... % a random count or a branch probability is not (sn_fj_validate).
                'ForkFanoutVector',...
                'Place', 'Transition', 'Linkage', 'Enabling', 'Inhibiting', 'Timing', 'Firing', 'Storage',...
                'Balking','Reneging','Retrial',...
                'LoadDependence','ClassDependence','JointDependence','GlobalDependence',...
                ... % c-server stations and binding buffers: the serial engine
                ... % walks the same State arms SolverCTMC declares and the NRM
                ... % honours both (test_ssa_nrm_closed_capacity)
                'MultiServer','FiniteCapacity'});
        end
        
        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverSSA.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end        
               
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()

            options = SolverOptions('SSA');
        end

        function libs = getLibrariesUsed(sn, options)
            % GETLIBRARIESUSED Get list of external libraries used by SSA solver
            % SSA uses internal simulation, no external libraries needed
            libs = {};
        end

    end
end
