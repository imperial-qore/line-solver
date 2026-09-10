classdef SolverQNS < NetworkSolver
    % Wrapper of LQNS's qnsolver
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = SolverQNS(model,varargin)
            % SELF = SolverQNS(MODEL,VARARGIN)

            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
        end

        runtime = run(self)

        function sn = getStruct(self)
            % QN = GETSTRUCT()

            % Get data structure summarizing the model
            sn = self.model.getStruct(false);
        end

        [runtime, analyzer] = runAnalyzer(self, options);

       function [allMethods] = listValidMethods(self)
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            sn = self.model.getStruct();
            allMethods = {'default','conway','rolia','zhou','suri','reiser','schmidt'};
        end
        function featSupported = getMethodFeatureSet(self, method) %#ok<INUSD>
            % All QNS methods share the solver-level feature envelope.
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
            featSupported = SolverQNS.getFeatureSet();
        end

        function [bool, reason] = supportsModelMethod(self, method)
            % Structural gate, on top of the feature envelope, for what no
            % registry name can state.
            %
            % A BINDING FINITE BUFFER first. NOTHING under the QNS tree reads
            % sn.cap or sn.classcap -- neither the JMVA document qnsolver reads
            % nor the LQN that QN2LQN writes has a buffer -- so a capped station
            % was solved as an unbounded one and the table reported the
            % unconstrained answer under this solver's name. SolverMVA,
            % SolverNC, SolverAG and SolverFLD gate the same way through the
            % same helper. Without it SolverAUTO.listValidMethods offered all
            % eight 'qns' method names on the BAS-blocking model of cqn_bas_blocking.
            %
            % THEN THE PATH. runAnalyzer sends a product-form or an open model
            % to qnsolver through writeJMVA and a closed non-product-form one
            % to SolverLQNS through QN2LQN, and the two carry different things.
            % On the qnsolver path runAnalyzer maps the method name onto
            % options.config.multiserver one to one, and 'qnsolver -m' knows
            % conway, reiser, rolia and zhou only, so 'suri' and 'schmidt' die
            % on a multiserver station: QNS_MULTISERVER_REFUSAL is the predicate
            % SOLVER_QNS asks at run time, so it is asked here too. (The former
            % comment here claimed the config value stayed 'default' under
            % those names and Conway answered; runAnalyzer sets it from the
            % name before SOLVER_QNS runs, so it does not.) The JMVA document
            % also has no fork or join and carries mean demands only, which
            % JMTMETHODREFUSAL states for the engine, exactly as writeJMVA asks
            % it before writing. On the lqns path QN2LQN writes AND forks and
            % joins as activity precedences and passes the full service law, so
            % the same model is served there and no path rule applies.
            %
            % IMMEDIATE FEEDBACK is refused on BOTH paths, in this solver's own
            % words (QNS_IMMFEED_REFUSAL, also asked by runAnalyzer): neither
            % document can keep a self-looping job on its server.
            % STRUCTURAL PREDICATE FIRST, as SolverBA does and for the same
            % reason: it names the station, the cap and the way out, where the
            % feature envelope can only say "(feature: FiniteCapacity)". Once
            % FiniteCapacity became a registry name on 2026-09-05 the base gate
            % started answering first and the useful sentence became unreachable.
            if isa(self.model, 'Network')
                [bool, reason] = NetworkSolver.checkBindingCapacity(self.model, 'SolverQNS');
                if ~bool
                    return
                end
            end
            [bool, reason] = supportsModelMethod@NetworkSolver(self, method);
            if ~bool || ~isa(self.model, 'Network')
                return
            end
            sn = self.model.getStruct();
            reason = qns_immfeed_refusal(sn);
            if ~isempty(reason)
                bool = false;
                return
            end
            if self.model.hasProductFormSolution() || self.model.hasOpenClasses()
                [bool, reason] = qns_multiserver_refusal(sn, method);
                if ~bool
                    return
                end
                reason = jmtMethodRefusal(sn, method, self.getOptions(), 'jmva');
                bool = isempty(reason);
            end
        end


    end

    methods (Static)

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()
            %
            % What BOTH paths of runAnalyzer carry. qnsolver reads the JMVA
            % document writeJMVA emits (a station type, a mean demand and a
            % visit count per chain) and lqns reads the LQN QN2LQN writes (a
            % host per station with its multiplicity and discipline, a task,
            % an entry per class, the service law, OR/AND forks from the
            % routing), so a construct is declared only when neither drops it.
            %
            % WITHDRAWN, and why: the Petri-net names, which writeJMVA ignores
            % and QN2LQN has no arm for, so a Place or Transition vanished; MAP
            % and MMPP2, whose correlation the JMVA document flattens to a rate
            % (undeclared, needsMapEnv now solves them through the random
            % environment image instead); every discipline outside INF, PS and
            % FCFS, since the JMVA document names no discipline and the .lqnx
            % schema knows no 'dps', 'siro', 'sept' or 'lcfs' (lqns rejects the
            % file), while HOL kept its name but lost the class priorities on
            % both paths; and the state-dependent routings (RROBIN, WRROBIN,
            % SQ), which mean visit counts cannot express. Fork/Join stay: the
            % lqns path writes them as precedences, and the qnsolver path,
            % which cannot, is refused by supportsModelMethod. The distributions
            % stay too, by the JMVA policy: a mean is admissible, and the mean-
            % only refusal at a FCFS station is jmtMethodRefusal's to state.
            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink',...
                'Source',...
                'Router',...
                'ClassSwitch',...
                'Delay',...
                'DelayStation',...
                'Queue',...
                'Fork',...
                'Join',...
                'Forker',...
                'Joiner',...
                'Logger',...
                'Coxian',...
                'Cox2',...
                'APH',...
                'Erlang',...
                'Exp',...
                'HyperExp',...
                'Det',...
                'Gamma',...
                'Lognormal',...
                'Normal',...
                'PH',...
                'Pareto',...
                'Weibull',...
                'Replayer',...
                'Trace',...
                'Uniform',...
                'StatelessClassSwitcher',...
                'InfiniteServer',...
                'SharedServer',...
                'Buffer',...
                'Dispatcher',...
                'Server',...
                'JobSink',...
                'RandomSource',...
                'ServiceTunnel',...
                'LogTunnel',...
                'SchedStrategy_INF',...
                'SchedStrategy_PS',...
                'SchedStrategy_FCFS',...
                'RoutingStrategy_PROB',...
                'RoutingStrategy_RAND',...
                'SchedStrategy_EXT',...
                'ClosedClass',...
                'OpenClass',...
                ... % c-server stations: the JMVA document carries the count as an
                ... % <ldstation> and the LQN as a host multiplicity; 'suri' and
                ... % 'schmidt' refuse one on the qnsolver path, which stays
                ... % structural (QNS_MULTISERVER_REFUSAL). FiniteCapacity is NOT
                ... % declared: neither document has a buffer (supportsModelMethod).
                'MultiServer'});
        end

        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)

            featUsed = model.getUsedLangFeatures();
            featSupported = SolverQNS.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end
    end

    methods (Static)
        function bool = hasLocalBinary()
            %HASLOCALBINARY True if the native qnsolver command is on the system path
            if ispc
                [~, ret] = dos('qnsolver -H');
                bool = ~contains(ret, 'not recognized', 'IgnoreCase', true);
            else
                [~, ret] = unix('qnsolver -H');
                bool = ~contains(ret, 'command not found', 'IgnoreCase', true);
            end
        end

        function bool = isAvailable()
            %ISAVAILABLE True if a local qnsolver binary is installed
            %
            % The binary is 'qnsolver' (not 'qns'). qnsolver ships with LQNS,
            % whose licence forbids redistribution, so LINE never runs it from a
            % container image; run-tests.sh --lqns-docker puts a shim on the PATH
            % when a containerised build is what you want to exercise.
            bool = SolverQNS.hasLocalBinary();
        end

        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()

            options = SolverOptions('QNS');
            options.timespan = [Inf,Inf];
        end
    end
end
