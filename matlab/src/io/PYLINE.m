classdef PYLINE
    % PYLINE  MATLAB-to-native-Python bridge (lang='python').
    %
    % Static helpers that construct a native line_solver (pure Python, no JVM)
    % model element-by-element through MATLAB's in-process Python interface
    % (py.*), run a native Python solver, and marshal results back. Mirrors the
    % structure of JLINE.m (line_to_jline/from_line_network/from_line_node/...)
    % with the jline.* Java backend replaced by py.line_solver.* and the Java
    % Matrix class replaced by numpy arrays.
    %
    % Requirements: MATLAB pyenv must point at a CPython where the in-tree
    % python/ line_solver package is importable. No JPype/JVM is used.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods (Static)

        function pynetwork = line_to_pyline(model)
            % LINE_TO_PYLINE  Top-level entry: LINE model -> py.line_solver model.
            switch class(model)
                case {'Network', 'MNetwork'}
                    pynetwork = PYLINE.from_line_network(model);
                case 'LayeredNetwork'
                    pynetwork = PYLINE.from_line_layered_network(model);
                otherwise
                    line_error(mfilename, sprintf('PYLINE (lang=python) does not support model class ''%s'' yet.', class(model)));
            end
        end

        function pynet = from_line_network(model)
            % FROM_LINE_NETWORK  Build a native Python Network from a LINE Network.
            PYLINE.assertPythonReady();
            L = py.importlib.import_module('line_solver');

            line_nodes = model.getNodes;
            line_classes = model.getClasses;
            nnodes = length(line_nodes);
            nclasses = length(line_classes);

            % see CLAUDE.md (lang='python' backend: Model bridge) for rationale.
            % An OI/PAS station joins them because its mu(c) is a MATLAB FUNCTION
            % HANDLE, which has no py.* form; linemodel_save already materializes
            % it as the 'oiServiceRate' macrostate table the native reader reads.
            if any(cellfun(@(nd) isa(nd, 'Cache') || isa(nd, 'Place') || isa(nd, 'Transition') ...
                    || (isa(nd, 'Queue') && ~isempty(nd.svcRateFun)), line_nodes)) ...
                    || PYLINE.needsJsonBridge(model, line_nodes)
                pynet = PYLINE.from_line_via_json(model);
                return
            end

            pynet = L.Network(model.getName);

            % 1) Nodes first (closed classes reference a station node)
            pynodes = cell(1, nnodes);
            for n = 1:nnodes
                if isa(line_nodes{n}, 'ClassSwitch') && line_nodes{n}.autoAdded
                    continue; % Python link() re-adds auto ClassSwitch nodes
                end
                if isa(line_nodes{n}, 'Join')
                    forkNode = pynodes{line_nodes{n}.joinOf.index};
                    pynodes{n} = PYLINE.from_line_node(line_nodes{n}, pynet, L, forkNode);
                else
                    pynodes{n} = PYLINE.from_line_node(line_nodes{n}, pynet, L, []);
                end
            end

            % 2) Classes
            pyclasses = cell(1, nclasses);
            for r = 1:nclasses
                pyclasses{r} = PYLINE.from_line_class(line_classes{r}, pynet, L, pynodes);
            end

            % 3) Service / arrival processes
            for n = 1:nnodes
                if isempty(pynodes{n})
                    continue;
                end
                PYLINE.set_service(line_nodes{n}, pynodes{n}, line_classes, pyclasses, L);
            end

            % 3b) Per-class station attributes that need the classes to exist
            for n = 1:nnodes
                if isempty(pynodes{n})
                    continue;
                end
                PYLINE.set_drop_rules(line_nodes{n}, pynodes{n}, pyclasses, L);
            end

            % 4) Routing
            PYLINE.from_line_links(model, pynet, pynodes, pyclasses, L);

            % 5) THE MODEL'S INITIAL STATE, marshalled LAST because link()
            % rebuilds the struct and would discard it.
            %
            % Nothing carried it before, so every lang='python' solve ran at the
            % NATIVE default initialization. That is invisible for a steady-state
            % mean off an ergodic chain, which is why it survived; it decides a
            % TRANSIENT outright. On init_state_fcfs_nonexp all three priors came
            % back as one number (0.171991) where the reference reports 0.175091,
            % 0.175972 and 0.175821 -- the bridge answering Prior 1's question
            % three times under three names, and only visible once PYLINE.tranAvg
            % existed to ask a transient question at all. This is the py.* twin
            % of the `initialState` key the JSON wire carries for the C++ and JAR
            % routes.
            %
            % A NEGATIVE ENTRY IS NOT A STATE but the "ignore this station" flag
            % of the getProb* family, filtered here as linemodel_save filters it:
            % read as a state it would start the chain where the model never is.
            % A STATE ROW IS AN ENCODING, NOT A VALUE, so the single-row case
            % crosses as the per-class MARGINAL and is re-encoded natively by
            % initFromMarginal. Sending MATLAB's row verbatim sends MATLAB's
            % encoding: on init_state_fcfs_nonexp, whose FCFS station carries
            % PH phases, the row initDefault leaves is read differently on the
            % native side and Prior 1 came back 0.171991 against the reference's
            % 0.175091, while the two initFromMarginal priors -- whose rows
            % happen to agree -- passed. The marginal says WHERE THE JOBS ARE,
            % which is what the reference means by an initial state, and each
            % codebase then builds its own row for it.
            hasPrior = false;
            for n = 1:nnodes
                if isempty(pynodes{n}) || ~isa(line_nodes{n}, 'StatefulNode')
                    continue
                end
                if numel(double(line_nodes{n}.statePrior(:))) > 1
                    hasPrior = true;
                end
            end
            if hasPrior
                % A genuine DISTRIBUTION over several rows: the native analyzer
                % iterates the support and weights the results, so the space and
                % the prior have to cross together and verbatim -- a marginal
                % cannot express a mixture.
                for n = 1:nnodes
                    if isempty(pynodes{n}) || ~isa(line_nodes{n}, 'StatefulNode')
                        continue
                    end
                    prior = double(line_nodes{n}.statePrior(:));
                    space = double(full(line_nodes{n}.space));
                    if isempty(space) || any(space(:) < 0) || size(space, 1) ~= numel(prior)
                        continue
                    end
                    pynodes{n}.setStateSpace(PYLINE.from_line_matrix(space));
                    pynodes{n}.setStatePrior(PYLINE.from_line_matrix(prior'));
                    pynodes{n}.setState(PYLINE.from_line_matrix(space(1, :)));
                end
            else
                snState = model.getStruct(true);  % sync node states into sn.state
                nir = zeros(snState.nstations, snState.nclasses);
                ok = true;
                for ist = 1:snState.nstations
                    isf = snState.stationToStateful(ist);
                    row = snState.state{isf};
                    % A NEGATIVE ENTRY IS NOT A STATE but the "ignore this
                    % station" flag of the getProb* family; nothing is sent when
                    % one is present, and the native side keeps its own default.
                    if isempty(row) || any(row(:) < 0)
                        ok = false;
                        break
                    end
                    [~, marg] = State.toMarginal(snState, snState.stationToNode(ist), row(1, :));
                    nir(ist, :) = double(marg(1, :));
                end
                if ok
                    pynet.initFromMarginal(PYLINE.from_line_matrix(nir));
                end
            end
        end

        function tf = needsJsonBridge(model, line_nodes)
            % TF = NEEDSJSONBRIDGE(MODEL, LINE_NODES)
            % True when a station carries a feature the element-by-element py.*
            % construction below does not marshal.
            %
            % SILENTLY DROPPING ONE IS THE WORST OUTCOME: from_line_node carries
            % the discipline, the server count, the capacity and the
            % load-dependent scaling and NOTHING ELSE, so a setup time, a server
            % breakdown, a switchover, a balking rule or a finite capacity region
            % crossed as a model that simply did not have it -- lqn_setup came
            % back with processor Util 0.75 against 0.436 because the setup layer
            % lost its setup. linemodel_save serializes every one of these, so
            % the whole model goes through the JSON bridge instead, exactly as a
            % Cache, an SPN or an OI/PAS station already does.
            tf = false;
            if isprop(model, 'regions') && ~isempty(model.regions)
                tf = true;
                return
            end
            props = {'setupTime', 'delayoffTime', 'switchoverTime', 'pollingType', ...
                'balkingStrategies', 'balkingThresholds', 'retrialDelays', ...
                'orbitImpatienceDistributions', 'patienceDistributions', ...
                'breakdownFailure', 'breakdownRepair', 'serverTypes', ...
                'serverParallelism', ...
                'immediateFeedback', 'lcdScaling', 'arrivalBatch'};
            for n = 1:numel(line_nodes)
                nd = line_nodes{n};
                for p = 1:numel(props)
                    if ~isprop(nd, props{p})
                        continue
                    end
                    v = nd.(props{p});
                    if iscell(v)
                        if any(~cellfun(@isempty, v(:)))
                            tf = true;
                            return
                        end
                    elseif ~isempty(v)
                        tf = true;
                        return
                    end
                end
            end
        end

        function pynode = from_line_node(line_node, pynet, L, forkNode)
            % FROM_LINE_NODE  LINE node -> py.line_solver node.
            if isa(line_node, 'Source')
                pynode = L.Source(pynet, line_node.getName);
            elseif isa(line_node, 'Sink')
                pynode = L.Sink(pynet, line_node.getName);
            elseif isa(line_node, 'Router')
                pynode = L.Router(pynet, line_node.getName);
            elseif isa(line_node, 'Delay')
                pynode = L.Delay(pynet, line_node.getName);
            elseif isa(line_node, 'Fork')
                pynode = L.Fork(pynet, line_node.getName);
                if ~isempty(line_node.output) && isprop(line_node.output, 'tasksPerLink') && ~isempty(line_node.output.tasksPerLink)
                    tpl = double(line_node.output.tasksPerLink);
                    % Standard forks use tasksPerLink=1 (the native default); only
                    % override when a quorum / task multiplier is actually set.
                    if any(tpl(:) ~= 1)
                        pynode.setTasksPerLink(PYLINE.from_line_matrix(tpl(:)'));
                    end
                end
            elseif isa(line_node, 'Join')
                pynode = L.Join(pynet, line_node.getName, forkNode);
            elseif isa(line_node, 'ClassSwitch')
                % Explicit ClassSwitch: the K x K switch-probability matrix lives
                % in server.csMatrix; the node routes class-preserving to/from its
                % neighbours (switching happens inside the node). Auto-added
                % ClassSwitch nodes are skipped upstream and re-added by link().
                csMatrix = line_node.server.csMatrix;
                pynode = L.ClassSwitch(pynet, line_node.getName, PYLINE.from_line_matrix(csMatrix));
            elseif isa(line_node, 'Queue')
                pysched = PYLINE.to_py_sched(line_node.schedStrategy, L);
                pynode = L.Queue(pynet, line_node.getName, pysched);
                nservers = line_node.getNumberOfServers;
                % INFINITE IS NOT intmax. Materializing it as 2^31-1 servers
                % turned the station into a multiserver whose per-class Util is
                % T*S/c ~ 1e-9 instead of the mean number busy, and SolverLN
                % reads that Util to set every INF task's think time -- which is
                % how a bridged LQN drifted from its native answer (tut11
                % DBProcessor Util 0.83 vs 0.994). The native Queue already
                % holds an infinite server count, so leave it alone.
                if isinf(nservers)
                    % nothing to set: the native default is already infinite
                elseif nservers > 1
                    pynode.setNumberOfServers(int32(nservers));
                end
                if ~isinf(line_node.cap)
                    pynode.setCapacity(int32(line_node.cap));
                end
                if ~isempty(line_node.lldScaling)
                    pynode.setLoadDependence(PYLINE.from_line_matrix(line_node.lldScaling));
                end
            else
                line_error(mfilename, sprintf('PYLINE (lang=python) does not support node ''%s'' (%s) yet.', line_node.getName, class(line_node)));
            end
        end

        function pysched = to_py_sched(schedId, L)
            % TO_PY_SCHED  SchedStrategy id -> py.line_solver.SchedStrategy enum.
            name = char(SchedStrategy.toProperty(SchedStrategy.toText(schedId)));
            try
                pysched = L.SchedStrategy.(name);
            catch
                line_error(mfilename, sprintf('PYLINE (lang=python) does not support the %s scheduling strategy yet.', name));
            end
        end

        function pyclass = from_line_class(line_class, pynet, L, pynodes)
            % FROM_LINE_CLASS  LINE job class -> py.line_solver class.
            % THE SIGNAL BRANCHES MUST COME FIRST. OpenSignal derives from
            % OpenClass and ClosedSignal from ClosedClass, so testing the bases
            % first matches a G-network signal and marshals it as an ORDINARY
            % class: the negative customers stop removing jobs and the bridge
            % silently solves a different model (on ag_gnetwork the signal class
            % picked up a queue length, a utilization and a throughput where
            % MATLAB has zero). A silent downgrade is worse than a refusal, and
            % nothing downstream can detect it.
            if isa(line_class, 'OpenSignal')
                pyclass = L.OpenSignal(pynet, line_class.getName, ...
                    PYLINE.to_py_signal_type(line_class.signalType, L), ...
                    int32(line_class.priority));
                PYLINE.set_signal_removal(line_class, pyclass, L);
            elseif isa(line_class, 'ClosedSignal')
                refnode = pynodes{line_class.refstat.index};
                pyclass = L.ClosedSignal(pynet, line_class.getName, ...
                    PYLINE.to_py_signal_type(line_class.signalType, L), ...
                    refnode, int32(line_class.priority));
                PYLINE.set_signal_removal(line_class, pyclass, L);
            elseif isa(line_class, 'OpenClass')
                pyclass = L.OpenClass(pynet, line_class.getName, int32(line_class.priority));
            elseif isa(line_class, 'SelfLoopingClass')
                refnode = pynodes{line_class.refstat.index};
                pyclass = L.SelfLoopingClass(pynet, line_class.getName, int32(line_class.population), refnode, int32(line_class.priority));
            elseif isa(line_class, 'ClosedClass')
                refnode = pynodes{line_class.refstat.index};
                pyclass = L.ClosedClass(pynet, line_class.getName, int32(line_class.population), refnode, int32(line_class.priority));
            else
                line_error(mfilename, sprintf('PYLINE (lang=python) does not support class type ''%s'' yet.', class(line_class)));
            end
        end

        function pytype = to_py_signal_type(signalType, L)
            % TO_PY_SIGNAL_TYPE  MATLAB SignalType constant -> py SignalType.
            % Bridged through the LOWERCASE TEXT both enums already use on the
            % JSON wire, never through the ordinal: MATLAB numbers REPLY 0 and
            % NEGATIVE 1, while the python member values are the strings, so an
            % ordinal handed over directly would silently pick another signal.
            pytype = L.SignalType(SignalType.toText(signalType));
        end

        function set_signal_removal(line_class, pyclass, L)
            % SET_SIGNAL_REMOVAL  Transfer how many jobs a signal removes, and
            % which. Empty removalDistribution means "exactly one", which is the
            % native default, so it is left alone rather than encoded.
            if ~isempty(line_class.removalDistribution)
                pyclass.setRemovalDistribution(PYLINE.from_line_distribution( ...
                    line_class.removalDistribution, L));
            end
            if ~isempty(line_class.removalPolicy)
                pyclass.setRemovalPolicy(L.RemovalPolicy( ...
                    RemovalPolicy.toText(line_class.removalPolicy)));
            end
        end

        function set_service(line_node, pynode, line_classes, pyclasses, L)
            % SET_SERVICE  Transfer arrival/service processes.
            if isa(line_node, 'Sink') || isa(line_node, 'Router') || isa(line_node, 'ClassSwitch') || ...
                    isa(line_node, 'Fork') || isa(line_node, 'Join') || isa(line_node, 'Cache') || isa(line_node, 'Logger')
                return;
            end
            for r = 1:length(line_classes)
                if isa(line_node, 'Source')
                    matlab_dist = line_node.getArrivalProcess(line_classes{r});
                    if isempty(matlab_dist) || isa(matlab_dist, 'Disabled')
                        continue;
                    end
                    pynode.setArrival(pyclasses{r}, PYLINE.from_line_distribution(matlab_dist, L));
                elseif isa(line_node, 'Queue') || isa(line_node, 'Delay')
                    matlab_dist = line_node.getService(line_classes{r});
                    if isempty(matlab_dist) || isa(matlab_dist, 'Disabled')
                        continue;
                    end
                    % THE SCHEDULING WEIGHT IS THE THIRD ARGUMENT OF setService,
                    % and it is model data with no derivation: under DPS/GPS it
                    % is the share the class is served at. Dropped, the native
                    % Python built the same station with every weight at 1, so
                    % cqn_scheduling_dps solved a DIFFERENT model -- CTMC and
                    % JMT alike disagreed with MATLAB by up to 47% on Queue2,
                    % which no tolerance should absorb.
                    %
                    % LPS IS EXCLUDED because its schedStrategyPar(1) is not a
                    % weight at all -- setLimit stores the concurrency limit
                    % there, and the native side keeps that in its own field.
                    if isa(line_node, 'Queue') && numel(line_node.schedStrategyPar) >= r ...
                            && SchedStrategy.toId(line_node.schedStrategy) ~= SchedStrategy.LPS
                        pynode.setService(pyclasses{r}, ...
                            PYLINE.from_line_distribution(matlab_dist, L), ...
                            line_node.schedStrategyPar(r));
                    else
                        pynode.setService(pyclasses{r}, PYLINE.from_line_distribution(matlab_dist, L));
                    end
                end
            end
        end

        function set_drop_rules(line_node, pynode, pyclasses, L)
            % SET_DROP_RULES  Transfer Station.dropRule, the per-class rule that
            % says what a full buffer does to an arriving job. It is user state
            % with no derivation, so a rule that fails to cross is invisible: on
            % cqn_bas_blocking the native SolverMVA then sees a finite capacity
            % with no blocking policy attached and REFUSES the model, because its
            % finite-capacity gate exempts Blocking-After-Service alone.
            if ~isa(line_node, 'Station') || isa(line_node, 'Source') || isa(line_node, 'Sink')
                return;
            end
            rules = line_node.dropRule;
            for r = 1:min(numel(rules), numel(pyclasses))
                if isempty(pyclasses{r}) || rules(r) == 0
                    continue; % 0 is the "never set" slot MATLAB leaves when growing the array
                end
                pynode.setDropRule(pyclasses{r}, PYLINE.to_py_drop_rule(rules(r), L));
            end
        end

        function pyrule = to_py_drop_rule(dropId, L)
            % TO_PY_DROP_RULE  MATLAB DropStrategy id -> py DropStrategy member.
            % Resolved through the member NAME. The two enums happen to agree
            % numerically today, but DropStrategy.toText returns prose ('BAS
            % blocking') that no python member is keyed by, so an ordinal is the
            % only other option and an ordinal is exactly what silently picks the
            % wrong member the day one side renumbers.
            switch dropId
                case DropStrategy.WAITQ
                    name = 'WAITQ';
                case DropStrategy.DROP
                    name = 'DROP';
                case DropStrategy.BAS
                    name = 'BAS';
                case DropStrategy.BBS
                    name = 'BBS';
                case DropStrategy.RSRD
                    name = 'RSRD';
                case DropStrategy.RETRIAL
                    name = 'RETRIAL';
                case DropStrategy.RETRIAL_WITH_LIMIT
                    name = 'RETRIAL_WITH_LIMIT';
                otherwise
                    line_error(mfilename, sprintf('PYLINE (lang=python) does not support drop rule %d yet.', dropId));
            end
            pyrule = L.DropStrategy.(name);
        end

        function pydist = from_line_distribution(line_dist, L)
            % FROM_LINE_DISTRIBUTION  LINE distribution -> py.line_solver process.
            if isa(line_dist, 'Exp')
                pydist = L.Exp(line_dist.getParam(1).paramValue);
            elseif isa(line_dist, 'Erlang')
                pydist = L.Erlang(line_dist.getParam(1).paramValue, int32(line_dist.getParam(2).paramValue));
            elseif isa(line_dist, 'HyperExp')
                pydist = L.HyperExp(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue, line_dist.getParam(3).paramValue);
            % THE FAMILY MUST SURVIVE THE CROSSING, not just the (alpha,T): a
            % featset names APH and Coxian and does NOT name PH, so collapsing
            % either one to L.PH made every solver reject its own model.
            elseif isa(line_dist, 'APH')
                alpha = line_dist.getParam(1).paramValue;
                T = line_dist.getParam(2).paramValue;
                pydist = L.APH(PYLINE.from_line_matrix(alpha(:)'), PYLINE.from_line_matrix(T));
            elseif isa(line_dist, 'PH')
                alpha = line_dist.getParam(1).paramValue;
                T = line_dist.getParam(2).paramValue;
                pydist = L.PH(PYLINE.from_line_matrix(alpha(:)'), PYLINE.from_line_matrix(T));
            elseif isa(line_dist, 'Coxian') % includes Cox2
                mu = line_dist.getParam(1).paramValue;
                phi = line_dist.getParam(2).paramValue;
                pydist = L.Coxian(PYLINE.from_line_matrix(mu(:)'), PYLINE.from_line_matrix(phi(:)'));
            elseif isa(line_dist, 'Det')
                pydist = L.Det(line_dist.getParam(1).paramValue);
            elseif isa(line_dist, 'Gamma')
                pydist = L.Gamma(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'Pareto')
                pydist = L.Pareto(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'Weibull')
                pydist = L.Weibull(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'Lognormal')
                pydist = L.Lognormal(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'Uniform')
                pydist = L.Uniform(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'Normal')
                pydist = L.Normal(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'MMPP2')
                pydist = L.MMPP2(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue, ...
                    line_dist.getParam(3).paramValue, line_dist.getParam(4).paramValue);
            elseif isa(line_dist, 'MAP')
                pydist = L.MAP(PYLINE.from_line_matrix(line_dist.D(0)), PYLINE.from_line_matrix(line_dist.D(1)));
            elseif isa(line_dist, 'MMPP')
                pydist = L.MAP(PYLINE.from_line_matrix(line_dist.D(0)), PYLINE.from_line_matrix(line_dist.D(1)));
            elseif isa(line_dist, 'ME')
                pydist = L.ME(PYLINE.from_line_matrix(line_dist.getParam(1).paramValue), PYLINE.from_line_matrix(line_dist.getParam(2).paramValue));
            elseif isa(line_dist, 'RAP')
                pydist = L.RAP(PYLINE.from_line_matrix(line_dist.getParam(1).paramValue), PYLINE.from_line_matrix(line_dist.getParam(2).paramValue));
            elseif isa(line_dist, 'NHPP')
                pydist = L.NHPP(PYLINE.from_line_matrix(line_dist.getBreakpoints()), PYLINE.from_line_matrix(line_dist.getRates()), logical(line_dist.isCyclic()));
            elseif isa(line_dist, 'MAPt') || isa(line_dist, 'PHt')
                if isa(line_dist, 'MAPt')
                    segA = line_dist.getD0Segments(); segB = line_dist.getD1Segments();
                else
                    segA = line_dist.getAlphaSegments(); segB = line_dist.getSSegments();
                end
                pyA = py.list(); pyB = py.list();
                for k = 1:numel(segA)
                    pyA.append(PYLINE.from_line_matrix(segA{k}));
                    pyB.append(PYLINE.from_line_matrix(segB{k}));
                end
                if isa(line_dist, 'MAPt')
                    pydist = L.MAPt(PYLINE.from_line_matrix(line_dist.getBreakpoints()), pyA, pyB, logical(line_dist.isCyclic()));
                else
                    pydist = L.PHt(PYLINE.from_line_matrix(line_dist.getBreakpoints()), pyA, pyB, logical(line_dist.isCyclic()));
                end
            elseif isa(line_dist, 'Zipf')
                pydist = L.Zipf(line_dist.getParam(3).paramValue, int32(line_dist.getParam(4).paramValue));
            elseif isa(line_dist, 'Trace') % before Replayer (Trace < Replayer)
                pydist = L.Trace(line_dist.params{1}.paramValue);
            elseif isa(line_dist, 'Replayer')
                pydist = L.Replayer(line_dist.params{1}.paramValue);
            elseif isa(line_dist, 'Bernoulli')
                pydist = L.Bernoulli(line_dist.getParam(1).paramValue);
            elseif isa(line_dist, 'Binomial')
                pydist = L.Binomial(int32(line_dist.getParam(1).paramValue), line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'Geometric')
                pydist = L.Geometric(line_dist.getParam(1).paramValue);
            elseif isa(line_dist, 'Poisson')
                pydist = L.Poisson(line_dist.getParam(1).paramValue);
            elseif isa(line_dist, 'DiscreteUniform')
                pydist = L.DiscreteUniform(int32(line_dist.getParam(1).paramValue), int32(line_dist.getParam(2).paramValue));
            elseif isa(line_dist, 'DiscreteSampler')
                pydist = L.DiscreteSampler(PYLINE.from_line_matrix(line_dist.getParam(1).paramValue), PYLINE.from_line_matrix(line_dist.getParam(2).paramValue));
            elseif isa(line_dist, 'Immediate')
                pydist = L.Immediate();
            elseif isempty(line_dist) || isa(line_dist, 'Disabled')
                pydist = L.Disabled();
            else
                line_error(mfilename, sprintf('PYLINE (lang=python) does not support distribution ''%s'' yet.', class(line_dist)));
            end
        end

        function from_line_links(model, pynet, pynodes, pyclasses, L)
            % FROM_LINE_LINKS  Reconstruct routing on the Python model.
            connections = model.getConnectionMatrix();
            line_nodes = model.getNodes;
            sn = model.getStruct;
            nclasses = length(pyclasses);

            % MATLAB node index -> Python node index (skip auto ClassSwitch)
            m2p = zeros(1, length(line_nodes));
            pidx = 0;
            for i = 1:length(line_nodes)
                if isa(line_nodes{i}, 'ClassSwitch') && line_nodes{i}.autoAdded
                    m2p(i) = -1;
                else
                    m2p(i) = pidx;
                    pidx = pidx + 1;
                end
            end

            useLinkMethod = ~isempty(sn.rtorig);
            hasAutoCS = any(cellfun(@(nd) isa(nd, 'ClassSwitch') && nd.autoAdded, line_nodes));

            if useLinkMethod
                rm = L.RoutingMatrix(pynet);
            end

            if useLinkMethod && hasAutoCS
                % Class-switching routing lives in sn.rtorig (station-indexed,
                % already excluding auto ClassSwitch nodes).
                for r = 1:nclasses
                    for s = 1:nclasses
                        Prs = sn.rtorig{r, s};
                        if isempty(Prs)
                            continue;
                        end
                        [nrows, ncols] = size(Prs);
                        for i = 1:nrows
                            for j = 1:ncols
                                if Prs(i, j) > 0
                                    rm.set(pyclasses{r}, pyclasses{s}, pynodes{i}, pynodes{j}, Prs(i, j));
                                end
                            end
                        end
                    end
                end
                pynet.link(rm);
                return
            end

            for i = 1:size(connections, 1)
                line_node = line_nodes{i};
                if isa(line_node, 'ClassSwitch') && line_node.autoAdded
                    continue;
                end
                pi_idx = m2p(i);
                for k = 1:nclasses
                    output_strat = line_node.output.outputStrategy{k};
                    strat = RoutingStrategy.fromText(output_strat{2});
                    switch strat
                        case RoutingStrategy.DISABLED
                            pynodes{i}.setRouting(pyclasses{k}, L.RoutingStrategy.DISABLED);
                        case RoutingStrategy.RAND
                            pynodes{i}.setRouting(pyclasses{k}, L.RoutingStrategy.RAND);
                            if ~useLinkMethod
                                for j = find(connections(i, :))
                                    if m2p(j) >= 0
                                        pynet.addLink(pynodes{i}, pynodes{j});
                                    end
                                end
                            end
                        case RoutingStrategy.PROB
                            pynodes{i}.setRouting(pyclasses{k}, L.RoutingStrategy.PROB);
                            if length(output_strat) >= 3
                                probs = output_strat{3};
                                for j = 1:length(probs)
                                    dest_idx = probs{j}{1}.index;
                                    if connections(i, dest_idx) ~= 0 && m2p(dest_idx) >= 0
                                        if useLinkMethod
                                            rm.set(pyclasses{k}, pyclasses{k}, pynodes{i}, pynodes{dest_idx}, probs{j}{2});
                                        else
                                            pynodes{i}.setProbRouting(pyclasses{k}, pynodes{dest_idx}, probs{j}{2});
                                        end
                                    end
                                end
                            end
                        case RoutingStrategy.RROBIN
                            if useLinkMethod
                                line_error(mfilename, 'RROBIN cannot be used together with the link() command.');
                            end
                            pynodes{i}.setRouting(pyclasses{k}, L.RoutingStrategy.RROBIN);
                            for j = find(connections(i, :))
                                if m2p(j) >= 0
                                    pynet.addLink(pynodes{i}, pynodes{j});
                                end
                            end
                        case RoutingStrategy.JSQ
                            pynodes{i}.setRouting(pyclasses{k}, L.RoutingStrategy.JSQ);
                            if ~useLinkMethod
                                for j = find(connections(i, :))
                                    if m2p(j) >= 0
                                        pynet.addLink(pynodes{i}, pynodes{j});
                                    end
                                end
                            end
                        case RoutingStrategy.WRROBIN
                            if useLinkMethod
                                line_error(mfilename, 'WRROBIN cannot be used together with the link() command.');
                            end
                            for j = find(connections(i, :))
                                if m2p(j) >= 0
                                    pynet.addLink(pynodes{i}, pynodes{j});
                                end
                            end
                            % output_strat{3} = list of {targetNode, weight} pairs
                            for j = 1:length(output_strat{3})
                                tgt = output_strat{3}{j}{1};
                                weight = output_strat{3}{j}{2};
                                pynodes{i}.setRouting(pyclasses{k}, L.RoutingStrategy.WRROBIN, pynodes{tgt.index}, weight);
                            end
                        case RoutingStrategy.SQ
                            if length(output_strat) >= 3 && ~isempty(output_strat{3})
                                dparam = output_strat{3}{1};
                                pynodes{i}.setRouting(pyclasses{k}, L.RoutingStrategy.SQ, int32(dparam));
                            else
                                pynodes{i}.setRouting(pyclasses{k}, L.RoutingStrategy.SQ);
                            end
                            if ~useLinkMethod
                                for j = find(connections(i, :))
                                    if m2p(j) >= 0
                                        pynet.addLink(pynodes{i}, pynodes{j});
                                    end
                                end
                            end
                        case RoutingStrategy.SDR
                            % Krzesinski state-dependent routing: the native
                            % side takes the whole subnetwork in one call, so
                            % the declared branch topology is translated from
                            % MATLAB node objects to the native nodes. Branch
                            % index 1 is the complement M-V and crosses as an
                            % empty list, keeping the paper's own numbering.
                            % See _kb/16-state-dependent-routing.md
                            if useLinkMethod
                                line_error(mfilename, 'State-dependent routing cannot be used together with the link() command.');
                            end
                            for j = find(connections(i, :))
                                if m2p(j) >= 0
                                    pynet.addLink(pynodes{i}, pynodes{j});
                                end
                            end
                            decl = output_strat{3}{1};
                            B = numel(decl.branch);
                            brCell = cell(1, B);
                            brCell{1} = py.list(cell(1, 0)); % the complement M-V
                            for b = 2:B
                                nb = numel(decl.branch{b});
                                centers = cell(1, nb);
                                for q = 1:nb
                                    cidx = decl.branch{b}{q}.index;
                                    % A center absent from the native model would
                                    % otherwise cross as an empty cell and be read
                                    % as a DIFFERENT branch, silently.
                                    if m2p(cidx) < 0
                                        line_error(mfilename, sprintf('Node ''%s'' takes part in state-dependent routing but is not present in the native model.', decl.branch{b}{q}.getName()));
                                    end
                                    centers{q} = pynodes{cidx};
                                end
                                brCell{b} = py.list(centers);
                            end
                            if m2p(decl.departure.index) < 0
                                line_error(mfilename, sprintf('The departure center ''%s'' of the state-dependent routing is not present in the native model.', decl.departure.getName()));
                            end
                            % level carries one entry per branch INDEX, the
                            % unused index 1 included, so it is the same length
                            % as branches on both sides.
                            pynodes{i}.set_state_dep_routing(pyclasses{k}, ...
                                pynodes{decl.departure.index}, ...
                                py.list(brCell), ...
                                py.list(num2cell(int32(decl.level(:)'))), ...
                                py.list(num2cell(double(decl.C(:)'))), ...
                                PYLINE.from_line_matrix(decl.d));
                        otherwise
                            line_error(mfilename, sprintf('PYLINE (lang=python) does not support the ''%s'' routing strategy.', output_strat{2}));
                    end
                end
            end

            if useLinkMethod
                pynet.link(rm);
            end

        end

        %% ---- Marshalling helpers ----

        function pyarr = from_line_matrix(matrix)
            % FROM_LINE_MATRIX  MATLAB double matrix -> numpy ndarray.
            np = py.importlib.import_module('numpy');
            if isempty(matrix)
                pyarr = np.zeros(py.tuple({int32(0), int32(0)}));
                return
            end
            [rows, cols] = size(matrix);
            % Build a nested Python list to avoid ambiguous auto-conversion,
            % then let numpy assemble the 2-D array.
            rowsCell = cell(1, rows);
            for r = 1:rows
                rowsCell{r} = py.list(num2cell(double(matrix(r, :))));
            end
            pyarr = np.array(py.list(rowsCell));
            % Preserve column-vector / row-vector shape when a singleton dim.
            if rows == 1 || cols == 1
                % ndarray METHOD, not np.reshape: MATLAB resolves a py.module
                % function whose name it also owns (reshape, zeros, ...) to its
                % own builtin, which then rejects the py.tuple shape.
                pyarr = pyarr.reshape(py.tuple({int32(rows), int32(cols)}));
            end
        end

        function matrix = from_pyline_matrix(pyarr)
            % FROM_PYLINE_MATRIX  numpy ndarray / scalar -> MATLAB double.
            if isa(pyarr, 'py.NoneType')
                matrix = [];
                return
            end
            np = py.importlib.import_module('numpy');
            matrix = double(np.asarray(pyarr, pyargs('dtype', 'float64')));
        end

        function pyopts = parseSolverOptions(options, solverName)
            % PARSESOLVEROPTIONS  MATLAB options struct -> native solver kwargs.
            % The native SolverXOptions constructors have per-solver parameter
            % lists, so build a candidate native-name->value map and forward
            % only the keys the target options class accepts.
            accepted = PYLINE.acceptedKwargs(solverName);

            cand = {}; % {nativeName, value} pairs
            % THE METHOD NAMES THE ALGORITHM, so it can never be defaulted
            % silently. Every solver reached through PYLINE.Solver already
            % carries it POSITIONALLY, so it is offered as a kwarg only to the
            % constructors that take no positional method -- SolverLN, whose
            % options were reaching the native side with 'method' stripped, so
            % lqn_moment3's `lnoptions2.method='moment3'` selected the DEFAULT
            % layered engine. Listing 'method' in acceptedKwargs for a solver
            % that also takes it positionally would pass it twice.
            if isfield(options, 'method') && ~isempty(options.method)
                cand = [cand, {'method', char(options.method)}];
            end
            if isfield(options, 'tol') && ~isempty(options.tol)
                cand = [cand, {'tol', options.tol}];
            end
            if isfield(options, 'iter_tol') && ~isempty(options.iter_tol)
                cand = [cand, {'iter_tol', options.iter_tol}];
            end
            if isfield(options, 'iter_max') && ~isempty(options.iter_max)
                % Different options classes name this max_iter vs iter_max.
                cand = [cand, {'max_iter', int32(options.iter_max), ...
                    'iter_max', int32(options.iter_max)}];
            end
            if isfield(options, 'seed') && ~isempty(options.seed)
                cand = [cand, {'seed', int32(options.seed)}];
            end
            % THE AMVA WARM START IS PART OF THE ANSWER, not a speed knob.
            % SolverLN reseeds every layer's MVA with the previous iteration's
            % chain queue lengths (SolverLN/analyze); dropping it made each
            % bridged layer restart cold, and the two fixed points parted
            % company from the third LN iteration (tut11_lqn_basics:
            % DBProcessor Util 0.994 native vs 0.830 through this bridge).
            if isfield(options, 'init_sol') && ~isempty(options.init_sol) ...
                    && isnumeric(options.init_sol)
                cand = [cand, {'init_sol', PYLINE.from_line_matrix(options.init_sol)}];
            end
            % A cutoff may be per station-class, not scalar: options.cutoff is a
            % (nstations x nclasses) matrix on models that truncate each queue
            % separately, and the native side reads either shape.
            if isfield(options, 'cutoff') && ~isempty(options.cutoff) ...
                    && all(isfinite(options.cutoff(:)))
                if isscalar(options.cutoff)
                    cand = [cand, {'cutoff', int32(options.cutoff)}];
                else
                    cand = [cand, {'cutoff', PYLINE.from_line_matrix(options.cutoff)}];
                end
            end
            if isfield(options, 'samples') && ~isempty(options.samples)
                cand = [cand, {'samples', int32(options.samples)}];
            end
            % THE TRANSIENT HORIZON IS PART OF THE ANSWER for every solver that
            % integrates one: a stage of a random environment is read over its
            % sojourn, so a stage solver given the default horizon instead of
            % the caller's answers a different question (see getEnvAvg).
            %
            % BUT AN INFINITE START IS A SENTINEL, NOT A HORIZON. MATLAB spells
            % "solve for the steady state" as timespan(1) = Inf -- it is what
            % every runAnalyzer tests with isinf(options.timespan(1)) -- and
            % Solver.defaultOptions carries [Inf,Inf]. Forwarded literally, the
            % native side reads it as a transient over a degenerate interval and
            % integrates nothing: init_state_ps, whose FLD is built from the
            % GENERIC defaults, came back with the uniform split QLen
            % [1.5 0.5 1.5 0.5] and Tput 4.5 against ArvR 1.8615 -- values that
            % do not even balance -- instead of QLen(1) = 0.33559. Reproduced in
            % native python with no bridge: absent, [0,Inf] and [0,1000] all give
            % 0.33559 and only [Inf,Inf] degenerates. A FINITE start is a real
            % transient request and still crosses, which is what the random
            % environment's per-stage [0,1e3] needs.
            if isfield(options, 'timespan') && numel(options.timespan) == 2 ...
                    && isfinite(options.timespan(1))
                cand = [cand, {'timespan', py.list({double(options.timespan(1)), ...
                    double(options.timespan(2))})}];
            end
            % Subprocess controls of the simulator wrappers. keep governs
            % whether the working directory survives the run, timeout the
            % child's deadline; Inf crosses as float('inf'), which is the
            % native default and means "wait".
            % The ODE integrator's own settings. `stiff` selects the stiff
            % solver and is TRUE by default, so a caller turning it off
            % (renv_threestages_repairmen) was being given the stiff one
            % anyway; `timestep` pins a fixed transient grid, and a transient
            % read on a different grid is a different set of points.
            if isfield(options, 'stiff') && ~isempty(options.stiff) ...
                    && (islogical(options.stiff) || isnumeric(options.stiff))
                cand = [cand, {'stiff', logical(options.stiff)}];
            end
            if isfield(options, 'timestep') && ~isempty(options.timestep) ...
                    && isnumeric(options.timestep) && isscalar(options.timestep)
                cand = [cand, {'timestep', double(options.timestep)}];
            end
            if isfield(options, 'keep') && ~isempty(options.keep) ...
                    && (islogical(options.keep) || isnumeric(options.keep))
                cand = [cand, {'keep', logical(options.keep)}];
            end
            if isfield(options, 'timeout') && ~isempty(options.timeout) ...
                    && isnumeric(options.timeout) && isscalar(options.timeout)
                cand = [cand, {'timeout', double(options.timeout)}];
            end
            % options.config CROSSES WHOLE, because the two codebases ported one
            % design and read the same key names ('multiserver', 'interlocking',
            % 'map_env_method', ...). Dropping it substituted an APPROXIMATION
            % silently: MVA(model, config.multiserver='softmin') on
            % cqn_repairmen_multi answers QLen 2.4671 at the Delay in MATLAB and
            % 1.7794 through this bridge, which is the DEFAULT rule's answer, so
            % lang='python' returned a different approximation under the same
            % name. Only the entries a py.dict can hold travel: a function
            % handle, struct or array has no dict form and is left behind rather
            % than mistranslated.
            cfgCrossed = {};
            if isfield(options, 'config') && isstruct(options.config) ...
                    && isscalar(options.config)
                cfgPairs = {};
                cfgFields = fieldnames(options.config);
                for ci = 1:numel(cfgFields)
                    cval = options.config.(cfgFields{ci});
                    if ischar(cval) || isstring(cval)
                        cfgPairs = [cfgPairs, {cfgFields{ci}, char(cval)}]; %#ok<AGROW>
                    elseif (isnumeric(cval) || islogical(cval)) && isscalar(cval)
                        cfgPairs = [cfgPairs, {cfgFields{ci}, cval}]; %#ok<AGROW>
                    end
                end
                if ~isempty(cfgPairs)
                    cand = [cand, {'config', py.dict(pyargs(cfgPairs{:}))}];
                    cfgCrossed = cfgPairs(1:2:end);
                end
            end
            % SOME CONFIG ENTRIES ARE TOP-LEVEL FIELDS NATIVELY. SolverSSAOptions
            % and SolverMAMOptions carry no config= parameter at all, and native
            % FLD reads options.pstar rather than config['pstar'], so an entry
            % MATLAB spells under options.config has to be LIFTED to a kwarg or
            % it does not cross even when 'config' itself is accepted. Each pair
            % below is a name the native options dataclass declares.
            liftedCfg = {'pstar','hide_immediate','aoi_preemption', ...
                'warmupfrac','space_max','softmin_alpha','odemaxstep','record_events'};
            % Names the caller NAMED rather than inherited from the defaults.
            % Only these are worth refusing over: a default that the native
            % class has no field for was never a request.
            explicit = {};
            if isfield(options, 'config') && isstruct(options.config) && isscalar(options.config)
                for ci = 1:numel(liftedCfg)
                    nm = liftedCfg{ci};
                    if isfield(options.config, nm) && ~isempty(options.config.(nm))
                        cval = options.config.(nm);
                        if islogical(cval)
                            cand = [cand, {nm, logical(cval)}]; %#ok<AGROW>
                            explicit{end+1} = nm; %#ok<AGROW>
                        elseif isnumeric(cval) && isscalar(cval)
                            cand = [cand, {nm, double(cval)}]; %#ok<AGROW>
                            explicit{end+1} = nm; %#ok<AGROW>
                        end
                    end
                end
            end
            % The simulator's event budget and confidence level. options.events
            % is NaN when unset, which means "use samples" and must not cross as
            % a NaN budget; options.confint is false, true (95%) or the level.
            if isfield(options, 'events') && ~isempty(options.events) ...
                    && isnumeric(options.events) && isscalar(options.events) ...
                    && isfinite(options.events)
                cand = [cand, {'events', int32(options.events)}];
                explicit{end+1} = 'events';
            end
            if isfield(options, 'confint') && ~isempty(options.confint)
                if islogical(options.confint)
                    if options.confint
                        cand = [cand, {'confidence_level', 0.95}];
                        explicit{end+1} = 'confidence_level';
                    end
                elseif isnumeric(options.confint) && isscalar(options.confint) ...
                        && options.confint > 0 && options.confint < 1
                    cand = [cand, {'confidence_level', double(options.confint)}];
                    explicit{end+1} = 'confidence_level';
                end
            end
            % VERBOSITY IS THE CALLER'S, not a constant. It was pinned to false
            % here, so a caller raising it saw nothing and a caller who had
            % silenced a noisy solver could not tell the bridge from the native
            % run. VerboseLevel is an integer level (SILENT=0), so anything
            % above SILENT is verbose.
            cand = [cand, {'verbose', PYLINE.verboseFlag(options)}];

            % AN ENTRY THAT CROSSED INSIDE `config` HAS CROSSED. The lift above
            % exists for the native classes that read a config key as a
            % top-level field, and it is the only route when the class declares
            % no `config` at all. Where the class DOES declare one, the whole
            % dict already carries the key and the failed lift is a second route
            % that was never needed: refusing there rejected every M2P row of a
            % solver whose MATLAB defaults set config.hide_immediate (54 rows of
            % the 2026-08-14 parity-static run), an AMBIENT default the comment
            % below already says must not be refused over.
            crossedInConfig = any(strcmp('config', accepted));
            args = {};
            unroutable = {};
            for i = 1:2:numel(cand)
                if any(strcmp(cand{i}, accepted))
                    args = [args, {cand{i}, cand{i+1}}]; %#ok<AGROW>
                elseif any(strcmp(cand{i}, explicit))
                    if crossedInConfig && any(strcmp(cand{i}, cfgCrossed))
                        continue;
                    end
                    unroutable{end+1} = cand{i}; %#ok<AGROW>
                end
            end
            if ~isempty(unroutable)
                % A DROPPED OPTION ANSWERS A DIFFERENT QUESTION UNDER THE SAME
                % NAME. Saying so beats running the default silently. Only a
                % field the caller NAMED reaches here; the ambient defaults
                % MATLAB fills in for every solver are still dropped quietly,
                % since they were never a request.
                line_error(mfilename, sprintf(['PYLINE (lang=''python'') cannot forward %s ' ...
                    'to %s: the native options class declares no such field. Remove the ' ...
                    'option, or add it to the native %sOptions dataclass.'], ...
                    strjoin(unique(unroutable), ', '), char(solverName), char(solverName)));
            end
            pyopts = pyargs(args{:});
        end

        function flag = verboseFlag(options)
            % VERBOSEFLAG  options.verbose (VerboseLevel or logical) -> logical.
            flag = false;
            if ~isfield(options, 'verbose') || isempty(options.verbose)
                return
            end
            v = options.verbose;
            if islogical(v)
                flag = any(v(:));
            elseif isnumeric(v)
                flag = any(v(:) > VerboseLevel.SILENT);
            end
        end

        function name = canonicalSolverName(solverName)
            % CANONICALSOLVERNAME  Short alias class -> canonical Solver* name.
            % `class(FLD(model))` is 'FLD', not 'SolverFluid': the short classes
            % SUBCLASS the canonical ones, so anything keyed on class() must
            % canonicalize first or a perfectly valid solver reads as unported.
            name = char(solverName);
            switch name
                case {'FLD','Fluid','SolverFLD'}
                    name = 'SolverFluid';
                case 'CTMC'
                    name = 'SolverCTMC';
                case 'MVA'
                    name = 'SolverMVA';
                case 'NC'
                    name = 'SolverNC';
                case 'MAM'
                    name = 'SolverMAM';
                case {'AG','SolverAG'}
                    name = 'SolverAG';
                case 'SSA'
                    name = 'SolverSSA';
                case 'JMT'
                    name = 'SolverJMT';
                case 'BA'
                    name = 'SolverBA';
                case {'AUTO','SolverAuto'}
                    name = 'SolverAUTO';
            end
        end

        function names = acceptedKwargs(solverName)
            % ACCEPTEDKWARGS  Native SolverXOptions constructor parameter names.
            %
            % READ OFF THE NATIVE DATACLASS, not transcribed. A hand-kept list
            % drifts the moment a native option is added or renamed, and the
            % failure is silent: the field is dropped and the default answers
            % under the caller's name. dataclasses.fields is the same source
            % the constructor itself uses, so the two cannot disagree.
            canon = PYLINE.canonicalSolverName(solverName);
            [modName, clsName] = PYLINE.nativeOptionsClass(canon);
            if isempty(clsName)
                % NOT a silent {'verbose'}. Falling through here means a
                % solver was bridged without stating which of its options
                % survive the crossing, and the default answer to that must
                % be "say so", not "drop them all" -- that default is what
                % let SolverLN reach the native side carrying nothing but
                % verbose. Mirrors CPPLINE's carry-or-refuse policy.
                line_error(mfilename, sprintf(['PYLINE (lang=''python'') has no option ' ...
                    'class for %s, so its options cannot be forwarded. Add a case to ' ...
                    'PYLINE.nativeOptionsClass naming the native %sOptions class.'], ...
                    char(solverName), char(solverName)));
            end
            PYLINE.assertPythonReady();
            dc = py.importlib.import_module('dataclasses');
            mod = py.importlib.import_module(modName);
            cls = py.getattr(mod, clsName);
            flds = cell(dc.fields(cls));
            names = cell(1, numel(flds));
            for i = 1:numel(flds)
                names{i} = char(flds{i}.name);
            end
            % Names that must never be forwarded as kwargs:
            %  - 'method' is POSITIONAL for every solver PYLINE.Solver builds;
            %    listing it here would pass it twice. SolverLN is the one
            %    constructor that takes no positional method, so it keeps it.
            %  - 'lang' and 'arith' select the backend of the NATIVE solver.
            %    Forwarding MATLAB's lang='python' into them would ask the
            %    native side to re-dispatch, which is the bridge's own job.
            drop = {'lang','arith'};
            if ~strcmp(canon, 'SolverLN')
                drop = [drop, {'method'}];
            end
            names = setdiff(names, drop, 'stable');
        end

        function [modName, clsName] = nativeOptionsClass(canon)
            % NATIVEOPTIONSCLASS  Canonical solver name -> native options class.
            modName = '';
            clsName = '';
            switch canon
                case 'SolverMVA'
                    modName = 'line_solver.solvers.solver_mva'; clsName = 'SolverMVAOptions';
                case 'SolverBA'
                    % Native SolverBA subclasses SolverMVA, so it takes the MVA
                    % options; 'level' rides on the method name, not on kwargs.
                    modName = 'line_solver.solvers.solver_mva'; clsName = 'SolverMVAOptions';
                case 'SolverNC'
                    modName = 'line_solver'; clsName = 'SolverNCOptions';
                case 'SolverCTMC'
                    modName = 'line_solver'; clsName = 'SolverCTMCOptions';
                case 'SolverMAM'
                    modName = 'line_solver'; clsName = 'SolverMAMOptions';
                case 'SolverAG'
                    % The agent-based solver's options are its own class: they
                    % carry maxStates and the execution backend, which the
                    % generic container has no field for.
                    modName = 'line_solver.solvers.solver_ag'; clsName = 'SolverAGOptions';
                case {'SolverFluid','SolverFLD'}
                    modName = 'line_solver'; clsName = 'SolverFLDOptions';
                case 'SolverSSA'
                    modName = 'line_solver'; clsName = 'SolverSSAOptions';
                case 'SolverJMT'
                    modName = 'line_solver'; clsName = 'SolverJMTOptions';
                case 'SolverLN'
                    modName = 'line_solver'; clsName = 'SolverLNOptions';
                case 'SolverAUTO'
                    % SolverAUTOOptions carries the selection knobs only; the
                    % chosen solver's own options come from the model it picks.
                    modName = 'line_solver'; clsName = 'SolverAUTOOptions';
            end
        end

        %% ---- Solver constructors ----

        function pysolver = Solver(name, pynet, options, L)
            % SOLVER  Dispatch to the native Python solver constructor by name.
            method = 'default';
            if isfield(options, 'method') && ~isempty(options.method)
                method = char(options.method);
            end
            name = PYLINE.canonicalSolverName(name);
            pyopts = PYLINE.parseSolverOptions(options, name);
            switch name
                case 'SolverMVA'
                    pysolver = L.SolverMVA(pynet, method, pyopts);
                case 'SolverNC'
                    pysolver = L.SolverNC(pynet, method, pyopts);
                case 'SolverCTMC'
                    pysolver = L.SolverCTMC(pynet, method, pyopts);
                case 'SolverMAM'
                    pysolver = L.SolverMAM(pynet, method, pyopts);
                case 'SolverAG'
                    % The agent-based (RCAT) solver. nativeOptionsClass already
                    % named its SolverAGOptions class, so the options crossed;
                    % without this case the dispatch fell through to the
                    % catch-all and every AG solve under lang='python' died with
                    % "does not support SolverAG yet" -- which is not one of the
                    % harness's refusal phrases, so the five ag_* [M2P] parity
                    % rows failed rather than reporting a named gap.
                    pysolver = L.SolverAG(pynet, method, pyopts);
                case {'SolverFluid', 'SolverFLD'}
                    pysolver = L.SolverFluid(pynet, method, pyopts);
                case 'SolverSSA'
                    pysolver = L.SolverSSA(pynet, method, pyopts);
                case {'SolverAUTO', 'SolverAuto'}
                    pysolver = L.SolverAUTO(pynet, method, pyopts);
                case 'SolverBA'
                    pysolver = L.SolverBA(pynet, method, pyopts);
                case 'SolverJMT'
                    % Native SolverJMT reads 'default' as the alias of 'jsim',
                    % so the MATLAB method name crosses unmapped.
                    pysolver = L.SolverJMT(pynet, method, pyopts);
                otherwise
                    line_error(mfilename, sprintf('PYLINE (lang=python) does not support %s yet.', name));
            end
        end

        function [QNt, UNt, TNt] = tranAvg(solverName, model, options)
            % [QNT, UNT, TNT] = TRANAVG(SOLVERNAME, MODEL, OPTIONS)
            % The transient means from the native getTranAvg(), in the contract
            % of self.result.Tran.Avg: (nstations x nclasses) cells whose entries
            % are [value, t] two-column matrices.
            %
            % THAT IS THE RESULT TABLE, NOT THE GETTER'S ANSWER, exactly as in
            % CPPLINE.tranAvg: getTranAvg returns a metricVal struct per cell and
            % builds it in @NetworkSolver/getTranAvg from these tables, so the
            % callers store what comes back and delegate there. One wrapping
            % implementation, not a native one and a bridge one that can differ
            % on the disabled case.
            %
            % A DISABLED PAIR STAYS EMPTY. The native getter returns None for a
            % (station, class) it has no series for, and filling it with zeros
            % here would report an idle station where there is no measurement.
            PYLINE.assertPythonReady();
            L = py.importlib.import_module('line_solver');
            pynet = PYLINE.line_to_pyline(model);
            pysolver = PYLINE.Solver(solverName, pynet, options, L);
            res = cell(pysolver.getTranAvg());
            QNt = PYLINE.tranGrid(res{1});
            UNt = PYLINE.tranGrid(res{2});
            TNt = PYLINE.tranGrid(res{3});
        end

        function v = pySeries(x)
            % V = PYSERIES(X)
            % One native 1-D series as a MATLAB row vector, WITHOUT crossing
            % the buffer protocol.
            %
            % `double(numpy_array)` binds MATLAB's own libmwbuffer, and inside
            % the embedded interpreter that import can fail outright --
            % "ImportError: PyCapsule_Import could not import module
            % libmwbuffer", raised where a transient should have been. It is
            % also dtype-sensitive: this interpreter warns that its numpy has
            % "broken support" for longdouble, and a series carrying one
            % converts for `t` and dies for `metric` in the same call. Going
            % through `tolist()` hands MATLAB native Python floats and no
            % buffer at all, so neither failure is reachable. A caller that
            % already has a plain list (the native getters return either) is
            % served by the same path.
            if isa(x, 'py.NoneType')
                v = [];
                return
            end
            if any(strcmp(methods(x), 'tolist'))
                x = x.tolist();
            end
            c = cell(py.list(x));
            v = zeros(1, numel(c));
            for k = 1:numel(c)
                v(k) = double(c{k});
            end
        end

        function G = tranGrid(pygrid)
            % G = TRANGRID(PYGRID)
            % One [M][K] nested list of native TranResult(t, metric) objects as
            % an (M x K) cell of [value, t] matrices.
            rows = cell(py.list(pygrid));
            G = cell(numel(rows), 0);
            for i = 1:numel(rows)
                cols = cell(py.list(rows{i}));
                if i == 1
                    G = cell(numel(rows), numel(cols));
                end
                for r = 1:numel(cols)
                    cell_ir = cols{r};
                    if isa(cell_ir, 'py.NoneType')
                        continue
                    end
                    t = PYLINE.pySeries(cell_ir.t);
                    v = PYLINE.pySeries(cell_ir.metric);
                    n = min(numel(t), numel(v));
                    G{i, r} = [v(1:n).', t(1:n).'];
                end
            end
        end

        function [QN, UN, RN, TN, AN, WN, runtime] = getAvg(solverName, model, options)
            % GETAVG  Build the native model+solver, run getAvg(), marshal back.
            % Guarded import: line_to_pyline below calls assertPythonReady, but
            % this line runs FIRST, so an unconfigured pyenv escaped as a raw
            % MATLAB:Python:PyException instead of the classified message. That
            % is not cosmetic -- the parity harness reads the message to tell a
            % missing environment (a named SKIP) from a wrong answer (a FAIL),
            % so every [M2P] row on a host without the native package was scored
            % as a failure of the solver under test.
            PYLINE.assertPythonReady();
            L = py.importlib.import_module('line_solver');
            Tstart = tic;
            pynet = PYLINE.line_to_pyline(model);
            pysolver = PYLINE.Solver(solverName, pynet, options, L);
            res = cell(pysolver.getAvg());
            % Native getAvg() returns (Q,U,R,T,A,W) as (M x R) numpy arrays.
            QN = PYLINE.from_pyline_matrix(res{1});
            UN = PYLINE.from_pyline_matrix(res{2});
            RN = PYLINE.from_pyline_matrix(res{3});
            TN = PYLINE.from_pyline_matrix(res{4});
            AN = PYLINE.from_pyline_matrix(res{5});
            WN = PYLINE.from_pyline_matrix(res{6});
            PYLINE.restoreCacheResults(model, pynet);
            % THE TWO CODEBASES DISAGREE ON WHAT THE 5th OUTPUT IS. MATLAB's AN
            % is the OFFERED arrival rate; native python's getAvg returns the
            % CARRIED rate and sets A = T outright (solvers/base.py getAvg),
            % recomputing the offered rate only when it builds the table. Passed
            % through unchanged, every station here would report ArvR = Tput --
            % on a closed 2-station cycle that swaps the two stations' arrival
            % rates, which is wrong in MATLAB's contract even though python's own
            % table is right. Translate on MATLAB's side, with MATLAB's own
            % routing-based derivation, so lang='python' answers the question
            % lang='matlab' was asked.
            % TN stands in for the throughput handles: sn_get_arvr_from_tput
            % only tests them for emptiness before deriving AN from the routing
            % matrix, and this path holds the throughputs themselves, not the
            % solver's handle objects.
            %
            % A SIMULATOR IS THE EXCEPTION: JMT MEASURES the arrival rate. Its
            % JSIM document declares an 'Arrival Rate' measure per station-class
            % and both codebases parse it (MATLAB getResults.m, native
            % _parse_jsim_results), so lang='matlab' reports the SIMULATED rate
            % while the routing derivation above reports a flow-balance estimate
            % of it. The two agree only up to simulation noise, which broke the
            % cross-codebase JMT invariant at fixed seed and sample count (1e-4
            % to 1e-3 on closed and open multiclass models, against 1e-13 on
            % every other metric). Native getAvg() drops the measure by setting
            % A = T like every other solver, so read it off the solver instead.
            % NaN means "JMT reported no such measure", which MATLAB's own
            % zero-initialized result store reports as 0.
            if strcmp(solverName, 'SolverJMT')
                AJMT = PYLINE.from_pyline_matrix(pysolver.getAvgArvR());
                if ~isempty(AJMT)
                    AJMT(isnan(AJMT)) = 0;
                    AN = AJMT;
                end
                % A FINITE-CAPACITY REGION IS NOT A STATION, so getAvg() does not
                % return it: its rows live past the last station, which is where
                % getAvgNode reads them from to fill the FCR pseudo-node. Without
                % them getAvgNodeTable finds every FCR metric zero and filters the
                % node out (fcr_mm1waitq[M2P], 'row FCR1 missing'). getAvgFcr
                % stacks Q,U,R,W,A,T, each (nregions x nclasses).
                FCR = PYLINE.from_pyline_matrix(pysolver.getAvgFcr());
                snm = model.getStruct();
                F = 0;
                if isfield(snm, 'nregions') && ~isempty(snm.nregions)
                    F = double(snm.nregions);
                end
                if F > 0 && size(FCR, 1) == 6 * F
                    QN = [QN; FCR(1:F, :)];
                    UN = [UN; FCR(F+1:2*F, :)];
                    RN = [RN; FCR(2*F+1:3*F, :)];
                    WN = [WN; FCR(3*F+1:4*F, :)];
                    AN = [AN; FCR(4*F+1:5*F, :)];
                    TN = [TN; FCR(5*F+1:6*F, :)];
                end
            elseif ~isempty(TN)
                AN = sn_get_arvr_from_tput(model.getStruct(), TN, TN);
            end
            runtime = toc(Tstart);
        end

        function restoreCacheResults(model, pynet)
            % RESTORECACHERESULTS  Copy each Cache node's solved hit/miss split
            % back from the native model onto the LINE node.
            %
            % A cache's hit and miss probabilities are a SOLVER RESULT, not model
            % state: the analyzer computes them and writes them onto the node,
            % and getAvgNode then rebuilds the hit-class and miss-class node
            % throughputs from them (sn_get_node_tput_from_tput). Solving in the
            % native engine leaves the LINE node's copy empty, so that
            % reconstruction silently falls back to the nodevisits split, which
            % is the 0.5/0.5 guess link() laid down before any cache was
            % analysed: on cache_replc_fifo (5 items, capacity 2) the hit and
            % miss throughputs came back 0.5/0.5 instead of 0.4/0.6. getHitRatio
            % on the LINE node returned nothing for the same reason.
            %
            % The struct is then refreshed HARD, because the visits that carry
            % the hit and miss classes are derived from the split: without it a
            % Router downstream of the cache still saw the 0.5/0.5 routing
            % (cache_replc_routing Router ArvR 1/1 instead of 0.8/1.2, and the
            % residence times scaled by the same ratio). Every MATLAB solver
            % that writes a hit probability does the same, see
            % SolverNC/runAnalyzer.m (refreshChains / refreshStruct(true)).
            if ~isa(model, 'Network')
                return
            end
            nodes = model.getNodes;
            restored = false;
            pysn = pynet.getStruct();
            pynodes = cell(py.list(pynet.getNodes()));
            byName = configureDictionary('string', 'cell');
            for k = 1:numel(pynodes)
                byName{char(pynodes{k}.getName())} = pynodes{k};
            end
            for ind = 1:numel(nodes)
                if ~isa(nodes{ind}, 'Cache')
                    continue
                end
                key = nodes{ind}.getName;
                if ~isKey(byName, key)
                    continue
                end
                pynode = byName{key};
                % THE SPLIT IS READ FROM THE STRUCT FIRST. Every native solver
                % writes sn.nodeparam[i].actualhitprob, but only some also write
                % it onto the node object (SolverSSA writes the node in its table
                % builder, not in runAnalyzer, so after a plain getAvg the node
                % is still empty while the struct already holds the simulated
                % 0.394). NaN entries are kept: they mean "this class does not
                % read the cache", and sn_get_node_tput_from_tput tests for them.
                [hit, miss, dhit, residt] = PYLINE.cacheSplitFromStruct(pysn, ind - 1);
                if isempty(hit)
                    hit = PYLINE.from_pyline_matrix(pynode.getHitRatio());
                    miss = PYLINE.from_pyline_matrix(pynode.getMissRatio());
                end
                % THE DELAYED-HIT SHARE AND THE RETRIEVAL LATENCY COME FROM THE
                % STRUCT FOR THE SAME REASON as the hit/miss split, and reading
                % them from the node alone silently dropped both: the native SSA
                % writes actualdelayedhitprob/actualresidt on sn.nodeparam and
                % never on the node, so getAvgCacheTable computed
                % ArvR = arvr*(missprob+delayedprob) with delayedprob = 0 and
                % reported retrieval_simple SSA ArvR 0.47165 (the bare miss
                % share) where the same engine run natively reports 0.56387.
                if isempty(dhit)
                    dhit = PYLINE.from_pyline_matrix(pynode.get_delayed_hit_ratio());
                end
                if isempty(residt)
                    residt = PYLINE.from_pyline_matrix(pynode.get_residt());
                end
                % ABSENT MUST CLEAR, NOT KEEP. Every one of the six results is
                % written on every solve, empty included, exactly as the
                % lang='java' path does (SolverMVA/runAnalyzer.m:68-75). Writing
                % only the non-empty ones leaves the PREVIOUS solver's value on
                % the node, and an example that runs LDES then MVA then reports
                % one table built from both: on retrieval_simple the MVA row
                % carried its own hit/miss beside LDES's delayed-hit 0.091368 and
                % LDES's latency 0.49958, so the three probabilities summed to
                % 1.09 and ArvR = arvr*(missprob+delayedprob) matched no solver.
                if isempty(hit)
                    hit = [];
                    miss = [];
                elseif isempty(miss)
                    miss = 1 - hit;
                end
                nodes{ind}.setResultHitProb(hit(:)');
                nodes{ind}.setResultMissProb(miss(:)');
                nodes{ind}.setResultDelayedHitProb(PYLINE.asRow(dhit));
                nodes{ind}.setResultHitProbList( ...
                    PYLINE.from_pyline_matrix(pynode.get_hit_ratio_by_list()));
                nodes{ind}.setResultItemProb( ...
                    PYLINE.from_pyline_matrix(pynode.get_item_prob()));
                nodes{ind}.setResultResidT(PYLINE.asRow(residt));
                restored = restored || ~isempty(hit);
            end
            if restored
                model.refreshStruct(true);
            end
        end

        function v = asRow(v)
            % ASROW  A per-class result as a row vector, [] when absent.
            % setResult* stores what it is given, so an empty here is what
            % CLEARS a stale value from an earlier solve on the same model.
            if ~isempty(v)
                v = v(:)';
            end
        end

        function [hit, miss, dhit, residt] = cacheSplitFromStruct(pysn, nodeIndex0)
            % CACHESPLITFROMSTRUCT  sn.nodeparam[i].actualhitprob /
            % actualmissprob / actualdelayedhitprob / actualresidt of the native
            % struct, each [] when this solver wrote it nowhere.
            [hit, miss, dhit, residt] = deal([]);
            entry = py.dict(pysn.nodeparam).get(py.int(nodeIndex0));
            if isa(entry, 'py.NoneType')
                return
            end
            hit = PYLINE.from_pyline_matrix(py.getattr(entry, 'actualhitprob', py.None));
            miss = PYLINE.from_pyline_matrix(py.getattr(entry, 'actualmissprob', py.None));
            dhit = PYLINE.from_pyline_matrix(py.getattr(entry, 'actualdelayedhitprob', py.None));
            residt = PYLINE.from_pyline_matrix(py.getattr(entry, 'actualresidt', py.None));
        end

        %% ---- JSON model bridge (advanced Network features, Environment) ----

        function pynet = from_line_via_json(model)
            % FROM_LINE_VIA_JSON  Bridge a model through the canonical line-model
            % JSON schema: MATLAB linemodel_save -> temp .json -> native
            % load_model. Used for models whose element-by-element construction
            % is impractical over py.* (Cache, SPN Place/Transition, Environment).
            PYLINE.assertPythonReady();
            L = py.importlib.import_module('line_solver');
            jsonfile = [tempname, '.json'];
            linemodel_save(model, jsonfile);
            cleanup = onCleanup(@() PYLINE.tryDelete(jsonfile)); %#ok<NASGU>
            pynet = L.load_model(jsonfile);
        end

        %% ---- LayeredNetwork (LQN) ----

        function pynet = from_line_layered_network(model)
            % FROM_LINE_LAYERED_NETWORK  LINE LayeredNetwork -> native LQN.
            % Bridged through the canonical LQN XML interchange format: the
            % element-by-element LQN graph (processors/tasks/entries/activities/
            % calls/precedences) is large and error-prone to marshal over py.*,
            % whereas writeXML/parseXML is a faithful, well-established round-trip.
            PYLINE.assertPythonReady();
            L = py.importlib.import_module('line_solver');
            xmlfile = [tempname, '.lqnx'];
            model.writeXML(xmlfile);
            cleanup = onCleanup(@() PYLINE.tryDelete(xmlfile));
            pynet = L.LayeredNetwork.parseXML(xmlfile);
            PYLINE.restoreNonRefThinkTimes(model, pynet);
        end

        function restoreNonRefThinkTimes(model, pynet)
            % RESTORENONREFTHINKTIMES  Reapply the one piece of LINE state that
            % the .lqnx interchange cannot carry. lqns rejects think-time on a
            % non-reference task, so writeXML omits it there; LINE nonetheless
            % gives such a think time to the task's callers as a delay, and
            % without this the bridged model would be a DIFFERENT model (the
            % layer of the called task loses its delay, which moves that layer's
            % throughput and every quantity derived from it).
            pytasks = cell(py.list(pynet.tasks));
            byName = configureDictionary('string', 'cell');
            for k = 1:numel(pytasks)
                byName{char(pytasks{k}.name)} = pytasks{k};
            end
            for t = 1:numel(model.tasks)
                task = model.tasks{t};
                if SchedStrategy.fromText(task.scheduling) == SchedStrategy.REF
                    continue
                end
                if isempty(task.thinkTimeMean) || task.thinkTimeMean <= 0
                    continue
                end
                if isKey(byName, task.name)
                    pytask = byName{task.name};
                    pytask.set_think_time(task.thinkTimeMean);
                end
            end
        end

        function tryDelete(f)
            % TRYDELETE  Best-effort temp-file removal.
            if exist(f, 'file')
                delete(f);
            end
        end

        function name = lnLayerSolverName(solver)
            % NAME = LNLAYERSOLVERNAME(SOLVER)
            % Class of the layer solver a MATLAB SolverLN was built with, or ''
            % when it kept the default. Under lang='python' no layer is
            % constructed, so the factory is probed on a throwaway network
            % rather than read off self.solvers -- the same resolution
            % CPPLINE.lnLayerSolver performs for --layer-solver. The layer
            % solver is what the fixed point is a fixed point OF, so leaving it
            % unnamed hands the native side a different question to answer.
            % CANONICALIZE: class(MVA(model)) is 'MVA', not 'SolverMVA', because
            % the short solver classes SUBCLASS the canonical ones. Both spellings
            % happen to resolve today (line_solver exports the short aliases too),
            % so this is fragility rather than a live fault -- but every other
            % consumer of a solver name on this bridge goes through
            % canonicalSolverName, and one that does not is the shape that made
            % PYLINE.Solver refuse a valid FLD as "not supported yet".
            name = '';
            if ~isempty(solver.solvers) && ~isempty(solver.solvers{1})
                name = PYLINE.canonicalSolverName(class(solver.solvers{1}));
                return
            end
            if isprop(solver, 'solverFactory') && ~isempty(solver.solverFactory) && ...
                    isa(solver.solverFactory, 'function_handle')
                name = PYLINE.canonicalSolverName(class(solver.solverFactory(CPPLINE.probeNetwork())));
            end
        end

        function pysolver = SolverLN(pynet, options, L, layerSolverName)
            % SOLVERLN  Native SolverLN over an LQN model.
            %
            % LAYERSOLVERNAME is the MATLAB class of the LAYER solvers the
            % caller built (SolverLN holds one instance per layer). Without it
            % the native side falls back to its own default layer solver, so
            % LN(model, @(l)NC(l,...)) under lang='python' silently answered
            % with AMVA layers: on lqn_twotasks that is E1 RespT 390.1 against
            % the exact 390. The native constructor takes the factory as its
            % second positional argument, the same shape as MATLAB's.
            opts = PYLINE.parseSolverOptions(options, 'SolverLN');
            if nargin < 4 || isempty(layerSolverName)
                pysolver = L.SolverLN(pynet, opts);
                return
            end
            facsrc = sprintf('lambda m: getattr(__import__("line_solver"), "%s")(m)', char(layerSolverName));
            pysolver = L.SolverLN(pynet, py.eval(facsrc, py.dict()), opts);
        end

        function [QN, UN, RN, TN, AN, WN, runtime] = getEnsembleAvg(model, options, elemNames, layerSolverName)
            % GETENSEMBLEAVG  Build the native LQN + SolverLN, run
            % get_ensemble_avg(), marshal the per-element metric vectors back.
            % Native get_ensemble_avg() returns (Q,U,R,T,A,W) as column vectors
            % positionally aligned with the NATIVE LQN element index, and the
            % native struct keeps a leading index-0 placeholder (empty name).
            %
            % THE TWO ELEMENT ORDERS ARE NOT THE SAME ORDER. MATLAB indexes the
            % activities in CREATION order, while the .lqnx interchange this
            % bridge goes through writes them grouped per task, so parseXML
            % rebuilds them in DOCUMENT order. On lcq_threehosts, which creates
            % A2 before A1, MATLAB's names run [... A2 A1 Ac Ac_hit Ac_miss] and
            % the native ones [... A1 Ac Ac_hit Ac_miss A2]: a positional copy
            % rotated the whole activity block by one and reported Ac_hit's
            % answer under Ac. Marshal by NAME, so the two orders never have to
            % agree. elemNames is self.lqn.names; a numeric argument is still
            % accepted and means "positional, nElem elements".
            PYLINE.assertPythonReady();
            L = py.importlib.import_module('line_solver');
            Tstart = tic;
            if nargin < 4
                layerSolverName = '';
            end
            pynet = PYLINE.from_line_layered_network(model);
            pysolver = PYLINE.SolverLN(pynet, options, L, layerSolverName);
            res = cell(pysolver.get_ensemble_avg());
            if isnumeric(elemNames)
                pyNames = {};
                mlNames = {};
                nElem = elemNames;
            else
                pyNames = PYLINE.pyNameList(pynet.getStruct().names);
                mlNames = cellstr(elemNames(:));
                nElem = numel(mlNames);
            end
            out = cell(1, 6);
            for k = 1:6
                out{k} = PYLINE.alignLNVector(PYLINE.from_pyline_matrix(res{k}), ...
                    nElem, pyNames, mlNames);
            end
            [QN, UN, RN, TN, AN, WN] = deal(out{:});
            runtime = toc(Tstart);
        end

        function [SensTable, sens] = getLNSensitivityTable(model, options, varargin)
            % GETLNSENSITIVITYTABLE  Build the native LQN + SolverLN and marshal
            % get_sensitivity_table() back into the MATLAB layer-wise table.
            % The native call runs the fixed-point loop itself when the layers
            % have not been solved, so the layer derivatives are taken at the
            % converged parameterization (see solver_ln.py getSensitivityTable).
            %
            % The name-value options are those of
            % @NetworkSolver/getSensitivityTable ('method', 'step', 'scheme');
            % an empty step forwards as None, which selects the native default.
            PYLINE.assertPythonReady();
            L = py.importlib.import_module('line_solver');
            pynet = PYLINE.from_line_layered_network(model);
            pysolver = PYLINE.SolverLN(pynet, options, L);

            method = 'auto'; step = []; scheme = 'forward';
            for a = 1:2:numel(varargin)
                switch lower(varargin{a})
                    case 'method', method = lower(varargin{a+1});
                    case 'step',   step = varargin{a+1};
                    case 'scheme', scheme = lower(varargin{a+1});
                    otherwise
                        line_error(mfilename, sprintf('Unknown option ''%s''.', varargin{a}));
                end
            end
            if isempty(step)
                pystep = py.None;
            else
                pystep = step;
            end
            T = pysolver.getSensitivityTable(pyargs('method', method, ...
                'step', pystep, 'scheme', scheme));

            Layer    = PYLINE.pyStringColumn(T, 'Layer');
            Station  = PYLINE.pyStringColumn(T, 'Station');
            JobClass = PYLINE.pyStringColumn(T, 'JobClass');
            dTput  = PYLINE.pyNumericColumn(T, 'dTput_dRate');
            dRespT = PYLINE.pyNumericColumn(T, 'dRespT_dRate');
            dQLen  = PYLINE.pyNumericColumn(T, 'dQLen_dRate');
            dUtil  = PYLINE.pyNumericColumn(T, 'dUtil_dRate');

            SensTable = table(Layer, Station, JobClass, dTput, dRespT, dQLen, dUtil, ...
                'VariableNames', {'Layer', 'Station', 'JobClass', 'dTput_dRate', ...
                'dRespT_dRate', 'dQLen_dRate', 'dUtil_dRate'});

            attrs = T.attrs;
            summary = char(attrs.get('method'));
            pyMethods = cell(py.list(attrs.get('layer_methods')));
            methods = cell(1, numel(pyMethods));
            for e = 1:numel(pyMethods)
                if isa(pyMethods{e}, 'py.NoneType')
                    methods{e} = '';
                else
                    methods{e} = char(pyMethods{e});
                end
            end
            SensTable.Properties.UserData = struct('method', summary, ...
                'layerMethods', {methods});

            pySens = cell(py.list(attrs.get('sens')));
            sens = cell(1, numel(pySens));
            for e = 1:numel(pySens)
                sens{e} = PYLINE.pySensStruct(pySens{e});
            end
        end

        function c = pyStringColumn(T, name)
            % PYSTRINGCOLUMN  pandas column of str -> column cell of char.
            vals = cell(py.list(T.get(name).tolist()));
            c = cell(numel(vals), 1);
            for i = 1:numel(vals)
                c{i} = char(vals{i});
            end
        end

        function v = pyNumericColumn(T, name)
            % PYNUMERICCOLUMN  pandas column of float -> column double vector.
            v = PYLINE.from_pyline_matrix(T.get(name).to_numpy());
            v = v(:);
        end

        function s = pySensStruct(pyobj)
            % PYSENSSTRUCT  native pfqn_sens result -> the MATLAB pfqn_sens
            % struct returned as the second output of getSensitivityTable.
            % Empty for a layer that took the finite-difference branch, which
            % carries no analytic Jacobian.
            if isempty(pyobj) || isa(pyobj, 'py.NoneType')
                s = [];
                return
            end
            s = struct();
            fields = {'X','Q','U','R','dX','dQ','dU','dR','QCov','QVar', ...
                'QTotVar','QCovAsym'};
            for k = 1:numel(fields)
                f = fields{k};
                if ~py.hasattr(pyobj, f)
                    continue
                end
                val = py.getattr(pyobj, f);
                if isa(val, 'py.NoneType')
                    continue
                end
                s.(f) = PYLINE.from_pyline_matrix(val);
            end
            if isfield(s, 'QCovAsym')
                s.QCovAsym = double(s.QCovAsym);
            end
            % Native params is a list of dicts with 0-based station/jobclass
            % (station -1 for a think-time parameter); pfqn_sens.m returns a
            % 1 x P struct array with 1-based .station (0 for Z) and .class.
            if py.hasattr(pyobj, 'params')
                pyParams = cell(py.list(py.getattr(pyobj, 'params')));
                P = numel(pyParams);
                params = struct('type', cell(1, P), 'station', cell(1, P), ...
                    'class', cell(1, P));
                for p = 1:P
                    d = pyParams{p};
                    params(p).type = char(d.get('type'));
                    st = double(d.get('station'));
                    if st < 0
                        params(p).station = 0;
                    else
                        params(p).station = st + 1;
                    end
                    params(p).class = double(d.get('jobclass')) + 1;
                end
                s.params = params;
            end
        end

        function [QN, UN, TN, runtime] = getEnvAvg(model, options, stageSolverName, stageOptions)
            % GETENVAVG  Build a native Environment (via JSON) + SolverENV and
            % marshal the environment-weighted (Q,U,T) metrics back.
            %
            % STAGESOLVERNAME is the MATLAB class of the stage solvers the
            % caller built (SolverENV holds one instance per stage) and
            % STAGEOPTIONS is that solver's own options struct. Both must be
            % forwarded: the stage solver is what decides the cost and the
            % accuracy of the whole ensemble, and substituting one here makes
            % lang='python' answer a different question from lang='matlab'.
            % Hardcoding SolverCTMC once turned the MVA random-environment
            % image of a two-class open FCFS model (mapEnvApprox, gallery
            % mmap1_multiclass) into a truncated CTMC at the default cutoff,
            % which exhausted the host. Dropping the stage options is just as
            % destructive on the other side: a stage FLD built here without the
            % caller's timespan integrates an UNSTABLE down-stage to the default
            % horizon instead of to 1e3, and renv_node_breakdown came back with
            % QLen 475.78 against MATLAB's 0.462.
            PYLINE.assertPythonReady();
            L = py.importlib.import_module('line_solver');
            Tstart = tic;
            penv = PYLINE.from_line_via_json(model);
            if nargin < 3 || isempty(stageSolverName)
                stageSolverName = 'SolverCTMC';
            end
            if nargin < 4
                stageOptions = struct();
            end
            stageSolverName = PYLINE.canonicalSolverName(stageSolverName);
            isCTMCStage = strcmp(stageSolverName, 'SolverCTMC');
            if isfield(options, 'method') && ~isempty(options.method)
                method = char(options.method);
            elseif isCTMCStage
                method = 'statevec'; % the documented bit-identical env analyzer
            else
                method = 'default';
            end
            % A CTMC stage needs a FINITE transient horizon; when the caller
            % pinned none, the ensemble's own timespan supplies it.
            if isCTMCStage && ~(isfield(stageOptions, 'timespan') ...
                    && numel(stageOptions.timespan) == 2 && isfinite(stageOptions.timespan(2)))
                Tend = 100;
                if isfield(options, 'timespan') && numel(options.timespan) == 2 && isfinite(options.timespan(2))
                    Tend = options.timespan(2);
                end
                stageOptions.timespan = [0, Tend];
            end
            % One native stage solver per stage model, built through the same
            % option marshalling every other bridged solver uses, so the method
            % and the timespan the caller pinned reach the stage.
            pstages = cell(py.list(penv.ensemble));
            psolvers = cell(1, numel(pstages));
            for e = 1:numel(pstages)
                psolvers{e} = PYLINE.Solver(stageSolverName, pstages{e}, stageOptions, L);
            end
            % The ENV fixed point stops where the CALLER's tolerance says it
            % does: the native default is iter_tol 1e-4, and an example pinning
            % 0.05 or 0.01 stops at a visibly different iterate.
            envArgs = {'method', method};
            if isfield(options, 'iter_max') && ~isempty(options.iter_max)
                envArgs = [envArgs, {'iter_max', int32(options.iter_max)}];
            end
            if isfield(options, 'iter_tol') && ~isempty(options.iter_tol)
                envArgs = [envArgs, {'iter_tol', double(options.iter_tol)}];
            end
            envArgs = [envArgs, {'verbose', false}];
            psolver = L.SolverENV(penv, py.list(psolvers), py.dict(pyargs(envArgs{:})));
            % getAvg returns (Q,U,R,T,A,W): throughput is the FOURTH output and
            % R is NaN throughout, exactly as in MATLAB SolverENV.getAvg. Taking
            % res{3} handed the caller NaN in place of every throughput.
            res = cell(psolver.getAvg());
            QN = PYLINE.from_pyline_matrix(res{1});
            UN = PYLINE.from_pyline_matrix(res{2});
            TN = PYLINE.from_pyline_matrix(res{4});
            runtime = toc(Tstart);
        end

        function v = trimLNVector(v, nElem)
            % TRIMLNVECTOR  Align a native LQN metric vector to nElem MATLAB
            % elements: drop the leading index-0 placeholder when present.
            v = v(:);
            if numel(v) == nElem + 1
                v = v(2:end);
            elseif numel(v) ~= nElem
                line_error(mfilename, sprintf('PYLINE (lang=python) LQN result length (%d) does not match the expected element count (%d or %d).', numel(v), nElem, nElem+1));
            end
        end

        function names = pyNameList(pyval)
            % PYNAMELIST  A native LQN struct name vector -> column cellstr.
            % The native struct stores it as a numpy array of str, whose first
            % entry is the empty index-0 placeholder.
            raw = cell(py.list(py.list(pyval)));
            names = cell(numel(raw), 1);
            for k = 1:numel(raw)
                names{k} = char(py.str(raw{k}));
            end
        end

        function out = alignLNVector(v, nElem, pyNames, mlNames)
            % ALIGNLNVECTOR  One native LQN metric vector in MATLAB's element
            % order. With names on both sides the values are looked up BY NAME;
            % with none (a numeric nElem was passed) it falls back to dropping
            % the index-0 placeholder positionally.
            v = v(:);
            if isempty(pyNames) || isempty(mlNames)
                out = PYLINE.trimLNVector(v, nElem);
                return
            end
            if numel(pyNames) ~= numel(v)
                line_error(mfilename, sprintf('PYLINE (lang=python) LQN result length (%d) does not match the native element-name count (%d).', numel(v), numel(pyNames)));
            end
            pos = configureDictionary('string', 'double');
            for k = 1:numel(pyNames)
                if ~isempty(pyNames{k})
                    pos(pyNames{k}) = k;
                end
            end
            out = zeros(nElem, 1);
            for e = 1:nElem
                if ~isKey(pos, mlNames{e})
                    line_error(mfilename, sprintf('PYLINE (lang=python) LQN element ''%s'' is absent from the native model.', mlNames{e}));
                end
                out(e) = v(pos(mlNames{e}));
            end
        end

        %% ---- Refusals ----

        function pyUnsupported(solverName, method, reason)
            % PYUNSUPPORTED(SOLVERNAME, METHOD, REASON)
            % Refuse a getter this bridge cannot serve, naming the getter and
            % the reason, exactly as CPPLINE.cppUnsupported does for lang='cpp'.
            % A getter that is not bridged must say so: the delegation in
            % runAnalyzerPreamble sets the AVERAGE results only, so a transient
            % getter that fell through found no self.result.Tran and died on a
            % dot-index rather than on a statement about the port.
            name = char(solverName);
            if strncmpi(name, 'Solver', 6)
                name = name(7:end);
            end
            line_error(mfilename, sprintf(['lang=''python'' cannot serve Solver%s.%s: %s. ' ...
                'Use lang=''matlab'' or lang=''java'' for it.'], name, method, reason));
        end

        %% ---- Environment checks ----

        function assertPythonReady()
            % ASSERTPYTHONREADY  Verify the native line_solver is importable and
            % pin the embedded interpreter's BLAS/OpenMP thread pools to a single
            % thread. Running numpy/scipy in-process (py.*) alongside MATLAB's own
            % worker threads otherwise causes BLAS oversubscription that spins or
            % deadlocks the session on heavy linear-algebra algorithms (e.g. CTMC).
            persistent ready tpHandle %#ok<PUSE>
            if ~isempty(ready) && ready
                return
            end
            % Set before the first py.* call, so the interpreter the probe in
            % bindSupportedInterpreter loads picks single-threaded BLAS from the
            % start rather than being clamped after numpy is already up.
            setenv('OMP_NUM_THREADS', '1');
            setenv('OPENBLAS_NUM_THREADS', '1');
            setenv('MKL_NUM_THREADS', '1');
            setenv('NUMEXPR_NUM_THREADS', '1');
            PYLINE.bindSupportedInterpreter();
            % THE BRIDGE BINDS TO THIS CHECKOUT'S python/, NOT TO WHATEVER
            % line_solver THE EMBEDDED INTERPRETER HAPPENS TO CARRY. Nothing in
            % the repository configured the embedded sys.path, so the import
            % resolved against the interpreter's site-packages and had two
            % failure modes, both of which the [M2P] parity row hit:
            %   * NO line_solver AT ALL -- a bare ModuleNotFoundError, i.e. the
            %     bridge is simply unusable on that host;
            %   * A FOREIGN line_solver SHADOWING THE TREE -- e.g. the retired
            %     JPype wrapper distribution (line-solver-wrapper) left in
            %     ~/.local/lib/python3.x/site-packages. That one imports, so the
            %     row RUNS, and answers from a different engine and a different
            %     vintage than the PYTHON row it is supposed to mirror.
            % The python-hosted parity rows never had either problem because
            % their runner puts line-dev.git/python on PYTHONPATH; the MATLAB
            % rows get no such environment, so the path is asserted here instead.
            pytree = fullfile(fileparts(lineRootFolder()), 'python');
            if exist(fullfile(pytree, 'line_solver'), 'dir') == 7
                syspath = py.sys.path();
                if ~any(cellfun(@(q) strcmp(char(q), pytree), cell(syspath)))
                    syspath.insert(int32(0), pytree);
                end
            end
            try
                mod = py.importlib.import_module('line_solver');
            catch ME
                line_error(mfilename, sprintf('lang=python could not import the native line_solver package via pyenv (%s). Set pyenv to a CPython where python/ is on sys.path.', ME.message));
            end
            % IMPORTED IS NOT THE SAME AS IMPORTED FROM HERE. A line_solver
            % already resident in the interpreter (imported before this call, or
            % pinned by a .pth) keeps its own origin no matter what is prepended
            % to sys.path, and a bridge answering from a foreign package is a
            % silent wrong answer rather than a missing one -- so say which file
            % was loaded and stop.
            if exist(fullfile(pytree, 'line_solver'), 'dir') == 7
                origin = '';
                try
                    % py.getattr rather than dot access: MATLAB's dot syntax
                    % does not reach a dunder attribute on a py.module.
                    origin = char(py.getattr(mod, '__file__'));
                end
                if ~isempty(origin) && ~strncmp(origin, pytree, numel(pytree))
                    line_error(mfilename, sprintf(['lang=python imported line_solver from %s, not from ' ...
                        'this checkout''s %s. A foreign line_solver (e.g. the retired JPype ' ...
                        'line-solver-wrapper distribution) is shadowing the native package; remove it ' ...
                        'from the pyenv interpreter''s site-packages, or point pyenv at an ' ...
                        'interpreter that does not carry it.'], origin, fullfile(pytree, 'line_solver')));
                end
            end
            % Runtime fallback: if numpy was already loaded with multithreaded
            % BLAS in this session, clamp the live thread pools too. The handle
            % is kept persistent so the limit is not reverted on GC.
            try
                tp = py.importlib.import_module('threadpoolctl');
                tpHandle = tp.threadpool_limits(pyargs('limits', int32(1)));
            catch
                % threadpoolctl unavailable; the env vars above cover fresh
                % interpreters. On an already-loaded multithreaded numpy without
                % threadpoolctl, restart MATLAB after setting the env vars.
            end
            ready = true;
        end

        function bindSupportedInterpreter()
            % BINDSUPPORTEDINTERPRETER  Ensure pyenv names a CPython that this
            % MATLAB release can actually load, and bind one that it can when it
            % does not.
            %
            % TWO THINGS WENT WRONG HERE AND EITHER ALONE IS ENOUGH TO LOSE THE
            % DIAGNOSIS. First, with nothing configured MATLAB auto-discovers
            % `python3` from PATH, and that is routinely NEWER than the release
            % supports -- every host of this cluster answers `python3 -V` with
            % 3.14 while R2026a loads at most 3.13, so the bridge had no usable
            % interpreter at all. Second, the guard meant to say so was DEAD
            % CODE: `pyenv` returns Executable as a STRING SCALAR, and
            % isempty("") is false for a 1x1 string, so `isempty(pe.Executable)`
            % never fired however unconfigured the environment was. Hence the
            % clear message was never reached and the first py.* call died with
            % MATLAB's own "Python commands require a supported version of
            % CPython" -- an error naming neither the interpreter nor the bridge,
            % which reached the [M2P] parity rows as "solver missing from
            % output". Emptiness is therefore tested on char(), never on the
            % string property.
            %
            % So the interpreter is asserted the way sys.path is asserted in
            % assertPythonReady: probe it, and when it cannot serve, bind one
            % that can. LINE_PYTHON pins an executable and is tried first.
            % pyenv PERSISTS its choice across MATLAB sessions, so a switch is
            % announced rather than made silently.
            if PYLINE.interpreterUsable()
                return
            end
            pe = pyenv;
            defaultExe = char(pe.Executable);
            if isempty(defaultExe)
                defaultExe = '<unset>';
            end
            % An interpreter already loaded into this MATLAB cannot be swapped:
            % pyenv('Version',...) errors while Status is Loaded/OutOfProcess.
            if ~strcmp(pe.Status, 'NotLoaded')
                line_error(mfilename, sprintf(['lang=python: the loaded Python (%s, version %s) cannot run ' ...
                    'py.* commands, and MATLAB cannot replace an interpreter once it is loaded. Restart ' ...
                    'MATLAB and set pyenv(''Version'', <supported CPython>), or set LINE_PYTHON, before solving.'], ...
                    defaultExe, char(pe.Version)));
            end
            tried = {};
            cands = PYLINE.candidateInterpreters();
            for k = 1:numel(cands)
                exe = cands{k};
                tried{end+1} = exe; %#ok<AGROW>
                try
                    pyenv('Version', exe);
                catch
                    continue % this MATLAB refuses that version outright
                end
                if PYLINE.interpreterUsable()
                    line_warning(mfilename, sprintf(['lang=python: the pyenv default (%s) is not a CPython ' ...
                        'this MATLAB release can load; bound %s instead. pyenv persists this choice; set ' ...
                        'LINE_PYTHON to pin a different interpreter.'], defaultExe, exe));
                    return
                end
            end
            if isempty(tried)
                triedMsg = 'no candidate interpreter was found on PATH';
            else
                triedMsg = sprintf('tried: %s', strjoin(tried, ', '));
            end
            line_error(mfilename, sprintf(['lang=python requires a MATLAB Python environment this release ' ...
                'can load; the pyenv default (%s) cannot run py.* commands and no alternative worked (%s). ' ...
                'Install a supported CPython, or set LINE_PYTHON / pyenv(''Version'',...) to one.'], ...
                defaultExe, triedMsg));
        end

        function tf = interpreterUsable()
            % INTERPRETERUSABLE  True when py.* commands actually execute.
            % Loading is the only reliable test: pyenv reports the CONFIGURED
            % executable, not whether MATLAB can bind it.
            %
            % THE PROBE MUST BE A PLAIN FUNCTION CALL, and `py.sys.version_info`
            % is not one. `sys.version_info` is a structseq TYPE, so MATLAB is
            % free to read `py.sys.version_info` as a CONSTRUCTOR and call it,
            % which raises `TypeError: cannot create 'sys.version_info'
            % instances` against a perfectly healthy interpreter. Which reading
            % MATLAB takes is not stable across a session: measured 2026-09-08,
            % every probe below answered before `dispatch_closed` ran and only
            % this one failed after it, with `py.str`, `py.list` and
            % `importlib.import_module` all still fine.
            %
            % The cost of getting it wrong is the whole rest of the session:
            % bindSupportedInterpreter reports "cannot run py.* commands", and
            % MATLAB cannot swap a loaded interpreter, so every later row sees
            % the same refusal. That is what skipped all 174 [M2P] parity rows
            % of run 20260908_104240 from `dispatch_closed` onwards, against an
            % interpreter that was answering the whole time.
            tf = false;
            pe = pyenv;
            if isempty(char(pe.Executable))
                return
            end
            try
                rlim = py.sys.getrecursionlimit(); %#ok<NASGU> the call itself forces the load
                tf = true;
            catch
                tf = false;
            end
        end

        function exes = candidateInterpreters()
            % CANDIDATEINTERPRETERS  Executables to try, newest first, when the
            % pyenv default cannot serve. No MATLAB-release-to-CPython table is
            % kept: a version this release refuses simply fails its probe in
            % bindSupportedInterpreter and the next candidate is tried, so the
            % list stays correct as both sides move.
            exes = {};
            pinned = getenv('LINE_PYTHON');
            if ~isempty(pinned)
                exes{end+1} = pinned;
            end
            names = {'python3.13','python3.12','python3.11','python3.10','python3.9'};
            if ispc
                lookup = 'where ';
            else
                lookup = 'command -v ';
            end
            for i = 1:numel(names)
                [rc, out] = system([lookup names{i}]);
                out = strtrim(out);
                if rc == 0 && ~isempty(out)
                    % `where` can return several matches, one per line.
                    lines = strsplit(out, newline);
                    exes{end+1} = strtrim(lines{1}); %#ok<AGROW>
                end
            end
            exes = unique(exes, 'stable');
        end

    end
end
