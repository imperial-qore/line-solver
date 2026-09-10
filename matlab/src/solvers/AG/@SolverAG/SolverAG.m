classdef SolverAG < NetworkSolver
    % Agent-based (RCAT) solver
    %
    % Solves a network by the Reversed Compound Agent Theorem: every
    % (station, class) pair becomes an isolated agent, and the agents are
    % coupled ONLY through the reversed rates of the synchronizing actions.
    % Agent k carries the generator
    %
    %   Q_k(x) = L_k + sum_{c : passive(c)=k} x_c * Pb_c
    %
    % and publishes, for every action a it is active on, a scalar x_a read off
    % its own stationary vector. The fixed point over that scalar vector is the
    % whole analysis, which is why it parallelises exactly (see the exec option
    % below) -- the per-sweep message is one double per action, not a
    % distribution.
    %
    % Methods: 'inap', 'inapplus', 'inapinf' (and the vestigial 'exact' alias,
    % which warns and falls back to 'inap'). Per Marin, Rota Bulo and Balsamo,
    % "A Numerical Algorithm for the Decomposition of Cooperating Structured
    % Markov Processes", MASCOTS 2012.
    %
    % These methods used to live in SolverMAM; they are the only algorithms in
    % LINE that read sn.issignal, so the G-network feature names belong here
    % and to no other solver.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = SolverAG(model,varargin)
            % SOLVERAG Create an agent-based (RCAT) solver instance
            %
            % @brief Creates an RCAT solver that decomposes the model into
            %        cooperating agents and solves the reversed-rate fixed point
            % @param model Network model to be analyzed
            % @param varargin Optional parameters (method, tolerance, exec, ...)
            % @return self SolverAG instance

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

        function [allMethods] = listValidMethods(self) %#ok<MANU>
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            allMethods = {'default','inap','inapplus','inapinf','exact'};
        end

        function featSupported = getMethodFeatureSet(self, method) %#ok<INUSD>
            % Every AG method is the same RCAT decomposition differing only in
            % how the reversed rate is read off the component, so they share one
            % feature envelope; the genuine restrictions (process type,
            % single-server) are structural and are applied in
            % supportsModelMethod below.
            featSupported = SolverAG.getFeatureSet();
        end

        function [bool, reason] = supportsModelMethod(self, method)
            % Structural gate for the RCAT decomposition. RCAT gives every
            % component a service-phase and an arrival-phase dimension, so any
            % law with a genuine (D0,D1) Markovian representation is admissible;
            % what is not is a law whose matrices are not a generator (ME, RAP),
            % one that is not time-homogeneous (NHPP, MAPt, PHt), one that is not
            % continuous-time (DMAP), and one that arrives in batches (BMAP,
            % MMAP), since a batch moves the level by more than one.
            % See _kb/06-solver-catalog.md for rationale.
            sn = self.model.getStruct();

            unsupp = false(size(sn.procid));
            for ist = 1:size(sn.procid, 1)
                for r = 1:size(sn.procid, 2)
                    unsupp(ist, r) = ~SolverAG.rcatSupportsProcess(sn.procid(ist, r));
                end
            end
            unsupp = unsupp & isfinite(sn.rates) & (sn.rates > 0);
            % A signal is a trigger with no service, so its service entry is
            % never read; only its Source arrival rate is, and that one must stay
            % exponential because the removal is folded into the component as a
            % scalar rate.
            if isfield(sn, 'issignal') && ~isempty(sn.issignal) && any(sn.issignal)
                issource = false(size(unsupp, 1), 1);
                for ist = 1:numel(issource)
                    issource(ist) = sn.nodetype(sn.stationToNode(ist)) == NodeType.Source;
                end
                sigcls = logical(sn.issignal(:))';
                unsupp(~issource, sigcls) = false;
                sigNonExp = (sn.procid ~= ProcessType.EXP) & isfinite(sn.rates) & (sn.rates > 0);
                sigNonExp(~issource, :) = false;
                sigNonExp(:, ~sigcls) = false;
                unsupp = unsupp | sigNonExp;
            end
            if any(unsupp(:))
                [ist, r] = find(unsupp, 1);
                bool = false;
                reason = sprintf(['The %s method supports processes with a Markovian ' ...
                    '(D0,D1) representation only (RCAT builds a phase dimension per ' ...
                    'component out of it), but station %d class %d is %s. Use SolverMAM ' ...
                    '(method ''dec.source'') for such models.'], method, ist, r, ...
                    ProcessType.toText(sn.procid(ist, r)));
                return;
            end
            % RCAT models every station single-server: reject multi-server;
            % see _kb/06-solver-catalog.md for rationale
            multi = isfinite(sn.nservers) & (sn.nservers > 1);
            if any(multi(:))
                ist = find(multi, 1);
                bool = false;
                reason = sprintf(['The %s method supports single-server stations only ' ...
                    '(RCAT does not model sn.nservers, so a multiserver station is driven ' ...
                    'at rho = lambda/mu instead of lambda/(c*mu)), but station %d has %d ' ...
                    'servers. Use SolverMAM (method ''dec.source'') for multiserver ' ...
                    'models.'], method, ist, sn.nservers(ist));
                return;
            end
            % An infinite-server station is decomposed as a single-server
            % component; AG_INF_SUPPORTS says when that is still exact, and
            % solver_ag raises with the same predicate.
            [infOk, infWhy] = ag_inf_supports(sn);
            if ~infOk
                bool = false;
                reason = infWhy;
                return;
            end
            % A FINITE BUFFER IS NOT SOMETHING RCAT CAN CARRY, and it was being
            % answered rather than refused: nothing under solvers/AG reads sn.cap
            % or sn.classcap, so a capped station was decomposed as an unbounded
            % one and the table reported the UNCONSTRAINED figures (a closed
            % 2-job tandem with cap 1 on the second queue returned the same
            % numbers with and without the cap, 1.09 jobs in a buffer of 1).
            %
            % The component's level bound is the class population by
            % construction (solver_ag: nlev = njobs(r)+1), and simply lowering it
            % to the buffer would not fix this: the top-level boundary is a
            % self-loop, which LOSES the arrival, whereas a closed job refused at
            % a full buffer must BLOCK the upstream departure. Blocking couples
            % the two components through more than the reversed rate of the
            % synchronizing action, which is exactly the independence RCAT
            % assumes, so this is a limit of the decomposition and not a bug in
            % its bookkeeping. Refuse, as SolverBA does through sn_has_blocking.
            % see _kb/06-solver-catalog.md
            [capOk, capWhy] = NetworkSolver.checkBindingCapacity(self.model, 'SolverAG');
            if ~capOk
                bool = false;
                reason = capWhy;
                return;
            end
            [bool, reason] = supportsModelMethod@NetworkSolver(self, method);
        end
    end

    methods (Static)

        function bool = rcatSupportsProcess(procid)
            % BOOL = RCATSUPPORTSPROCESS(PROCID)
            %
            % True when the RCAT analyzers can give a process of this type a
            % phase dimension: after sn_nonmarkov_toph (which AG runs with
            % phfit='ph' and preserveDet=false) it must hold a genuine (D0,D1)
            % pair with non-negative off-diagonal rates and a single arrival
            % per epoch.
            %
            % The list is an ALLOW-list on purpose: a process type nobody has
            % checked against this construction must be refused, not answered.
            % Compare BY NAME across codebases, never by the raw integer.
            if isnan(procid)
                bool = true;   % nothing to serve (e.g. a Transition mode)
                return;
            end
            bool = any(procid == [ProcessType.EXP, ProcessType.ERLANG, ...
                ProcessType.HYPEREXP, ProcessType.PH, ProcessType.APH, ...
                ProcessType.COXIAN, ProcessType.COX2, ProcessType.MAP, ...
                ProcessType.MMPP2, ProcessType.DET, ProcessType.UNIFORM, ...
                ProcessType.GAMMA, ProcessType.PARETO, ProcessType.WEIBULL, ...
                ProcessType.LOGNORMAL, ProcessType.REPLAYER, ...
                ProcessType.IMMEDIATE, ProcessType.DISABLED]);
        end

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()

            % NO FORK-JOIN. The decomposition has no synchronizing action for
            % a Join, so a fork-join model would be answered as probabilistic
            % routing through a pass-through Join; the 'Fork-Join support (via
            % FJ_codes)' this list used to carry was SolverMAM's, copied when
            % the RCAT methods moved here. A Delay (INF station) stays
            % declared, under the structural rule of AG_INF_SUPPORTS: every
            % component is a single-server QBD, which is the infinite-server
            % chain only for a closed class of population one.
            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink', 'Source', ...
                'Delay', 'DelayStation', 'Queue', ...
                'APH', 'Coxian', 'Erlang', 'Exp', 'HyperExp', 'MAP', 'MMPP2', ...
                'Det','Gamma','Lognormal','Pareto','Uniform','Weibull',...
                'StatelessClassSwitcher', 'InfiniteServer', ...
                'SharedServer', 'Buffer', 'Dispatcher', ...
                'Server', 'JobSink', 'RandomSource', 'ServiceTunnel', ...
                'SchedStrategy_INF', 'SchedStrategy_PS', ...
                'SchedStrategy_FCFS', ...
                'RoutingStrategy_PROB', 'RoutingStrategy_RAND', ...
                'ClosedClass', ...
                ... % a self-looping class is its own single-station component
                ... % (solver_ag reads sn.isslc); the qn sanity goldens pin the
                ... % inap/inapplus/inapinf rows on the *slc* models. MultiServer
                ... % and FiniteCapacity are deliberately absent: supportsModelMethod
                ... % words both refusals (RCAT drives rho = lambda/mu, no buffer).
                'SelfLoopingClass', ...
                'OpenClass', ...
                'OpenSignal', 'ClosedSignal', ... % G-network signals (solver_ag)
                'SignalType_NEGATIVE', 'SignalType_CATASTROPHE', ...
                'SignalBatchRemoval'}); % AG reads sn.signalremdist
        end

        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)

            featUsed = model.getUsedLangFeatures();
            featSupported = SolverAG.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end

        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = SolverOptions('AG');
        end
    end
end
