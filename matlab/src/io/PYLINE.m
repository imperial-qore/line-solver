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

            % see CLAUDE.md (lang='python' backend: Model bridge) for rationale
            if any(cellfun(@(nd) isa(nd, 'Cache') || isa(nd, 'Place') || isa(nd, 'Transition'), line_nodes))
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

            % 4) Routing
            PYLINE.from_line_links(model, pynet, pynodes, pyclasses, L);
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
                if isinf(nservers)
                    pynode.setNumberOfServers(int32(intmax('int32')));
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
            if isa(line_class, 'OpenClass')
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
                    pynode.setService(pyclasses{r}, PYLINE.from_line_distribution(matlab_dist, L));
                end
            end
        end

        function pydist = from_line_distribution(line_dist, L)
            % FROM_LINE_DISTRIBUTION  LINE distribution -> py.line_solver process.
            if isa(line_dist, 'Exp')
                pydist = L.Exp(line_dist.getParam(1).paramValue);
            elseif isa(line_dist, 'Erlang')
                pydist = L.Erlang(line_dist.getParam(1).paramValue, int32(line_dist.getParam(2).paramValue));
            elseif isa(line_dist, 'HyperExp')
                pydist = L.HyperExp(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue, line_dist.getParam(3).paramValue);
            elseif isa(line_dist, 'APH') || isa(line_dist, 'PH')
                alpha = line_dist.getParam(1).paramValue;
                T = line_dist.getParam(2).paramValue;
                pydist = L.PH(PYLINE.from_line_matrix(alpha(:)'), PYLINE.from_line_matrix(T));
            elseif isa(line_dist, 'Coxian') % includes Cox2
                [alpha, T] = PYLINE.coxian_to_ph(line_dist);
                pydist = L.PH(PYLINE.from_line_matrix(alpha(:)'), PYLINE.from_line_matrix(T));
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

        function [alpha, T] = coxian_to_ph(line_dist)
            % COXIAN_TO_PH  Bidiagonal (alpha,T) PH representation of a Coxian.
            mu = line_dist.getParam(1).paramValue;
            phi = line_dist.getParam(2).paramValue;
            mu = mu(:)';
            phi = phi(:)';
            k = length(mu);
            alpha = zeros(1, k);
            alpha(1) = 1;
            T = zeros(k, k);
            for i = 1:k
                T(i, i) = -mu(i);
                if i < k
                    T(i, i+1) = mu(i) * (1 - phi(i));
                end
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
                        case RoutingStrategy.RL
                            % The native line_solver exposes no API to carry an RL
                            % value function / action map, so it cannot reproduce a
                            % learned RL policy. Fail clearly rather than silently
                            % solving a different (default) policy.
                            line_error(mfilename, 'PYLINE (lang=python) cannot bridge RL routing: the native line_solver has no value-function API to reproduce the learned policy. Use lang=''matlab'' for RL routing.');
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
                pyarr = np.reshape(pyarr, py.tuple({int32(rows), int32(cols)}));
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
            if isfield(options, 'cutoff') && ~isempty(options.cutoff) && isfinite(options.cutoff)
                cand = [cand, {'cutoff', int32(options.cutoff)}];
            end
            if isfield(options, 'samples') && ~isempty(options.samples)
                cand = [cand, {'samples', int32(options.samples)}];
            end
            cand = [cand, {'verbose', false}];

            args = {};
            for i = 1:2:numel(cand)
                if any(strcmp(cand{i}, accepted))
                    args = [args, {cand{i}, cand{i+1}}]; %#ok<AGROW>
                end
            end
            pyopts = pyargs(args{:});
        end

        function names = acceptedKwargs(solverName)
            % ACCEPTEDKWARGS  Native SolverXOptions constructor parameter names.
            switch solverName
                case 'SolverMVA'
                    names = {'max_iter','tol','verbose','seed','cutoff','samples'};
                case 'SolverNC'
                    names = {'tol','iter_max','iter_tol','verbose','seed','cutoff','samples'};
                case 'SolverCTMC'
                    names = {'tol','cutoff','seed','samples','verbose'};
                case 'SolverMAM'
                    names = {'tol','max_iter','verbose'};
                case {'SolverFluid','SolverFLD'}
                    names = {'tol','iter_max','iter_tol','verbose','seed','cutoff','samples'};
                case 'SolverSSA'
                    names = {'tol','samples','seed','cutoff','verbose'};
                otherwise
                    names = {'verbose'};
            end
        end

        %% ---- Solver constructors ----

        function pysolver = Solver(name, pynet, options, L)
            % SOLVER  Dispatch to the native Python solver constructor by name.
            method = 'default';
            if isfield(options, 'method') && ~isempty(options.method)
                method = char(options.method);
            end
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
                case {'SolverFluid', 'SolverFLD'}
                    pysolver = L.SolverFluid(pynet, method, pyopts);
                case 'SolverSSA'
                    pysolver = L.SolverSSA(pynet, method, pyopts);
                case 'SolverAuto'
                    pysolver = L.SolverAuto(pynet, method, pyopts);
                otherwise
                    line_error(mfilename, sprintf('PYLINE (lang=python) does not support %s yet.', name));
            end
        end

        function [QN, UN, RN, TN, AN, WN, runtime] = getAvg(solverName, model, options)
            % GETAVG  Build the native model+solver, run getAvg(), marshal back.
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
            runtime = toc(Tstart);
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
            byName = containers.Map('KeyType', 'char', 'ValueType', 'any');
            for k = 1:numel(pytasks)
                byName(char(pytasks{k}.name)) = pytasks{k};
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
                    pytask = byName(task.name);
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

        function pysolver = SolverLN(pynet, options, L)
            % SOLVERLN  Native SolverLN over an LQN model.
            pysolver = L.SolverLN(pynet, PYLINE.parseSolverOptions(options, 'SolverLN'));
        end

        function [QN, UN, RN, TN, AN, WN, runtime] = getEnsembleAvg(model, options, nElem)
            % GETENSEMBLEAVG  Build the native LQN + SolverLN, run
            % get_ensemble_avg(), marshal the per-element metric vectors back.
            % Native get_ensemble_avg() returns (Q,U,R,T,A,W) as column vectors
            % positionally aligned with the LQN element index. The native LQN
            % struct keeps a leading index-0 placeholder (empty name), so the
            % vectors are length nidx+1; MATLAB self.lqn.names has nidx entries.
            % nElem = numel(self.lqn.names) is used to drop that placeholder so
            % the vectors line up 1:1 with the MATLAB element order.
            L = py.importlib.import_module('line_solver');
            Tstart = tic;
            pynet = PYLINE.from_line_layered_network(model);
            pysolver = PYLINE.SolverLN(pynet, options, L);
            res = cell(pysolver.get_ensemble_avg());
            QN = PYLINE.trimLNVector(PYLINE.from_pyline_matrix(res{1}), nElem);
            UN = PYLINE.trimLNVector(PYLINE.from_pyline_matrix(res{2}), nElem);
            RN = PYLINE.trimLNVector(PYLINE.from_pyline_matrix(res{3}), nElem);
            TN = PYLINE.trimLNVector(PYLINE.from_pyline_matrix(res{4}), nElem);
            AN = PYLINE.trimLNVector(PYLINE.from_pyline_matrix(res{5}), nElem);
            WN = PYLINE.trimLNVector(PYLINE.from_pyline_matrix(res{6}), nElem);
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

        function [QN, UN, TN, runtime] = getEnvAvg(model, options)
            % GETENVAVG  Build a native Environment (via JSON) + SolverENV and
            % marshal the environment-weighted (Q,U,T) metrics back. Uses the
            % statevec method (the documented bit-identical env analyzer); its
            % inner per-stage solver is a CTMC with a finite transient timespan.
            L = py.importlib.import_module('line_solver');
            Tstart = tic;
            penv = PYLINE.from_line_via_json(model);
            method = 'statevec';
            if isfield(options, 'method') && ~isempty(options.method)
                method = char(options.method);
            end
            Tend = 100;
            if isfield(options, 'timespan') && numel(options.timespan) == 2 && isfinite(options.timespan(2))
                Tend = options.timespan(2);
            end
            facsrc = sprintf('lambda m: __import__("line_solver").SolverCTMC(m, timespan=[0.0,%.15g], verbose=False)', Tend);
            fac = py.eval(facsrc, py.dict());
            psolver = L.SolverENV(penv, fac, py.dict(pyargs('method', method)));
            res = cell(psolver.getAvg());
            QN = PYLINE.from_pyline_matrix(res{1});
            UN = PYLINE.from_pyline_matrix(res{2});
            TN = PYLINE.from_pyline_matrix(res{3});
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
            pe = pyenv;
            if strcmp(pe.Status, 'NotLoaded') && isempty(pe.Executable)
                line_error(mfilename, 'lang=python requires a configured MATLAB Python environment (see pyenv).');
            end
            % Set before the first numpy import so a fresh interpreter picks
            % single-threaded BLAS from the start.
            setenv('OMP_NUM_THREADS', '1');
            setenv('OPENBLAS_NUM_THREADS', '1');
            setenv('MKL_NUM_THREADS', '1');
            setenv('NUMEXPR_NUM_THREADS', '1');
            try
                py.importlib.import_module('line_solver');
            catch ME
                line_error(mfilename, sprintf('lang=python could not import the native line_solver package via pyenv (%s). Set pyenv to a CPython where python/ is on sys.path.', ME.message));
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

    end
end
