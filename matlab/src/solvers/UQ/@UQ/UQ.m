classdef UQ < EnsembleSolver
    % UQ Solver wrapper for models with Prior distributions
    %
    % UQ detects Prior distributions in a model, expands the model
    % into a family of concrete networks (one per Prior alternative), solves
    % each using the specified solver, and aggregates results using prior
    % probabilities as weights.
    %
    % @brief Solver wrapper for Bayesian-style uncertainty analysis
    %
    % Key characteristics:
    % - Detects and expands Prior distributions in models
    % - Orchestrates multiple solver runs
    % - Aggregates results with prior-weighted expectations
    % - Provides posterior distribution access
    % - Reports the support-only (interval) range through getInterval, which
    %   ignores the weights and keeps only the endpoints of each Prior
    %
    % Example:
    % @code
    % model = Network('UncertainService');
    % source = Source(model, 'Source');
    % queue = Queue(model, 'Queue', SchedStrategy.FCFS);
    % sink = Sink(model, 'Sink');
    %
    % class = OpenClass(model, 'Jobs');
    % source.setArrival(class, Exp(1.0));
    % queue.setService(class, Prior({Exp(1), Exp(2)}, [0.5, 0.5]));
    %
    % model.link(model.serialRouting(source, queue, sink));
    %
    % post = UQ(model, @SolverMVA);
    % avgTable = post.getAvgTable();            % Prior-weighted expectations
    % postTable = post.getPosteriorTable();     % Per-alternative results
    % postDist = post.getPosteriorDist('R', queue, class); % Response time distribution
    % ival = post.getIntervalTable();           % Support-only range, no weights
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Constant)
        % Cap on the tensor-product design size. A design point is one full
        % solver run, so this bounds the cost of a quadrature design over
        % several Priors; beyond it the Monte Carlo design is the right tool.
        MaxDesignPoints = 4096;
    end

    properties
        solverFactory;      % Function handle to create solvers: @(model) SolverXXX(model)
        priorInfo;          % Struct array with Prior detection info, one entry per Prior
        design;             % Struct array of design points: .weight and .dists (one per Prior)
        aggregatedResult;   % Prior-weighted aggregate metrics
        originalModel;      % Reference to original model with Prior
    end

    methods
        function [bool, reason] = supportsModelMethod(self, method)
            % [BOOL, REASON] = SUPPORTSMODELMETHOD(METHOD)
            %
            % UQ's refusal IN ITS OWN WORDS, which the base gate cannot supply
            % here. That one falls back to SUPPORTS(MODEL), which answers with
            % a bare logical, and a caller left to reconstruct a reason from
            % the feature set gets a list of every feature the model uses: this
            % class's set declares the Prior UQ ITSELF consumes, every other
            % feature being the inner solver's to accept (see SUPPORTS), so
            % comparing it against a model answers a question nobody asked.
            %
            % The rule is the one SUPPORTS states: a model with no uncertain
            % parameter is not a UQ model. On top of it, the tensor-product
            % design of the discrete/quadrature methods must fit
            % MaxDesignPoints, a rule buildDesign used to state alone on the
            % run path; 'montecarlo' has no such cap.
            bool = UQ.modelHasPrior(self.model);
            if bool
                reason = '';
                if ~strcmp(UQ.resolveUQMethod(method), 'montecarlo')
                    [bool, reason] = self.tensorDesignRefusal();
                end
            else
                reason = ['SolverUQ analyses a model with an UNCERTAIN parameter: it expands ' ...
                    'the Prior into one model per design point and runs an inner solver on ' ...
                    'each. No service or arrival process of this model is a Prior, so there ' ...
                    'is one design point and its posterior is the point estimate the inner ' ...
                    'solver already returns.'];
            end
        end

        function self = UQ(model, solverFactory, varargin)
            % UQ Create a UQ solver wrapper
            %
            % @brief Creates a UQ wrapper for uncertainty analysis
            % @param model Network model (may contain Prior distributions)
            % @param solverFactory Function handle: @(m) SolverXXX(m) or solver class name
            % @param varargin Optional solver options
            % @return self UQ instance

            self@EnsembleSolver(model, mfilename);

            % Handle solver factory - accept class name or function handle
            if isa(solverFactory, 'function_handle')
                self.solverFactory = solverFactory;
            elseif ischar(solverFactory) || isstring(solverFactory)
                % Convert solver class name to factory
                className = char(solverFactory);
                self.solverFactory = str2func(['@(m) ', className, '(m)']);
            else
                line_error(mfilename, 'solverFactory must be a function handle or solver class name');
            end

            if ~isempty(varargin)
                self.setOptions(Solver.parseOptions(varargin, UQ.defaultOptions));
            else
                self.setOptions(UQ.defaultOptions);
            end

            self.priorInfo = [];
            self.design = [];
            self.ensemble = {};
            self.solvers = {};
            self.results = {};
            self.aggregatedResult = [];
            self.originalModel = model;

            % Detect and validate Prior usage
            self.detectPriors();
        end

        function detectPriors(self)
            % DETECTPRIORS Find Prior distributions in the model
            %
            % Scans all nodes for Prior distributions and stores location info.
            % Currently supports only a single Prior in the model.

            priors = {};
            model = self.model;
            nodes = model.getNodes();
            classes = model.getClasses();

            for i = 1:length(nodes)
                node = nodes{i};

                % Check service distributions (Queue, Delay)
                if isa(node, 'Queue') || isa(node, 'Delay')
                    for c = 1:length(classes)
                        try
                            dist = node.getService(classes{c});
                            if ~isempty(dist) && isa(dist, 'Prior')
                                priors{end+1} = struct(...
                                    'type', 'service', ...
                                    'node', node, ...
                                    'nodeName', node.name, ...
                                    'nodeIdx', i, ...
                                    'class', classes{c}, ...
                                    'className', classes{c}.name, ...
                                    'classIdx', c, ...
                                    'prior', dist);
                            end
                        catch
                            % Service not set for this class
                        end
                    end
                end

                % Check arrival distributions (Source)
                if isa(node, 'Source')
                    for c = 1:length(classes)
                        if isa(classes{c}, 'OpenClass')
                            try
                                if ~isempty(node.arrivalProcess) && ...
                                   length(node.arrivalProcess) >= 1 && ...
                                   size(node.arrivalProcess, 2) >= c && ...
                                   ~isempty(node.arrivalProcess{1, c})
                                    dist = node.arrivalProcess{1, c};
                                    if isa(dist, 'Prior')
                                        priors{end+1} = struct(...
                                            'type', 'arrival', ...
                                            'node', node, ...
                                            'nodeName', node.name, ...
                                            'nodeIdx', i, ...
                                            'class', classes{c}, ...
                                            'className', classes{c}.name, ...
                                            'classIdx', c, ...
                                            'prior', dist);
                                    end
                                end
                            catch
                                % Arrival not set
                            end
                        end
                    end
                end
            end

            if isempty(priors)
                self.priorInfo = [];
            else
                self.priorInfo = [priors{:}];
            end
        end

        function buildDesign(self)
            % BUILDDESIGN Reduce the detected Priors to weighted design points
            %
            % Each design point assigns one concrete Distribution to every
            % Prior in the model and carries the weight of that joint
            % assignment. Discrete and quadrature designs take the tensor
            % product of the per-Prior alternatives, so their weights are the
            % products of the marginal weights: this is the product-density
            % case of the joint f(theta_1,...,theta_l) in Trivedi and Bobbio
            % (2017), Eq. (3.67), and assumes the Priors are independent.
            % Monte Carlo instead draws all Priors jointly, so its cost is
            % independent of the number of Priors.

            if isempty(self.priorInfo)
                self.design = struct('weight', 1, 'dists', {{}});
                return;
            end

            L = length(self.priorInfo);
            method = self.getUQMethod();
            n = self.getUQNodes();

            if strcmp(method, 'montecarlo')
                self.design = repmat(struct('weight', 1/n, 'dists', {{}}), 1, n);
                for i = 1:n
                    dists = cell(1, L);
                    for l = 1:L
                        [d, ~] = self.priorInfo(l).prior.discretize(1, 'montecarlo');
                        dists{l} = d{1};
                    end
                    self.design(i).weight = 1/n;
                    self.design(i).dists = dists;
                end
                return;
            end

            % Tensor product of the per-Prior alternatives
            margDists = cell(1, L);
            margWeights = cell(1, L);
            for l = 1:L
                [margDists{l}, margWeights{l}] = self.priorInfo(l).prior.discretize(n, 'quadrature');
            end

            counts = cellfun(@length, margDists);
            total = prod(counts);
            % the predicate supportsModelMethod asks, so the run raises the
            % sentence the report showed
            [okd, whyd] = self.tensorDesignRefusal(counts);
            if ~okd
                line_error(mfilename, whyd);
            end

            self.design = repmat(struct('weight', 1, 'dists', {{}}), 1, total);
            for i = 1:total
                idx = UQ.unrankIndex(i, counts);
                dists = cell(1, L);
                w = 1;
                for l = 1:L
                    dists{l} = margDists{l}{idx(l)};
                    w = w * margWeights{l}(idx(l));
                end
                self.design(i).weight = w;
                self.design(i).dists = dists;
            end
        end

        function method = getUQMethod(self)
            % METHOD = GETUQMETHOD()
            % Resolve the discretization method from the solver options.
            %
            % 'default' keeps the historical behaviour: discrete Priors are
            % expanded as given, continuous Priors are discretized by
            % quadrature. The name table is RESOLVEUQMETHOD, shared with the
            % gate so the two cannot disagree on a method name.
            method = 'quadrature';
            if isfield(self.options, 'method') && ~isempty(self.options.method)
                method = UQ.resolveUQMethod(self.options.method);
            end
        end

        function n = getUQNodes(self)
            % N = GETUQNODES()
            % Number of nodes per continuous Prior, from options.samples.
            n = 11;
            if isfield(self.options, 'samples') && ~isempty(self.options.samples)
                n = round(self.options.samples);
            end
            if n < 1
                line_error(mfilename, 'options.samples must be at least 1');
            end
        end

        function [ok, reason] = tensorDesignRefusal(self, counts)
            % [OK, REASON] = TENSORDESIGNREFUSAL(COUNTS)
            % Whether the tensor-product design of the discrete/quadrature
            % methods fits MaxDesignPoints. COUNTS, the per-Prior alternative
            % counts, is derived from the Priors when omitted. Asked by
            % supportsModelMethod, so model.help refuses 'quadrature' on a
            % model whose Priors multiply out past the cap and still offers
            % 'montecarlo', and by buildDesign, so the run raises the same
            % sentence.
            ok = true;
            reason = '';
            if isempty(self.priorInfo)
                return
            end
            if nargin < 2 || isempty(counts)
                n = self.getUQNodes();
                L = length(self.priorInfo);
                counts = zeros(1, L);
                for l = 1:L
                    counts(l) = length(self.priorInfo(l).prior.discretize(n, 'quadrature'));
                end
            end
            total = prod(counts);
            if total > UQ.MaxDesignPoints
                ok = false;
                reason = sprintf(['Tensor-product design has %d points, above the limit of %d. ', ...
                    'Use options.method = ''montecarlo'' or reduce options.samples.'], ...
                    total, UQ.MaxDesignPoints);
            end
        end

        function hasPrior = hasPriorDistribution(self)
            % HASPRIOR = HASPRIORDISTRIBUTION()
            % Return true if model contains a Prior distribution
            hasPrior = ~isempty(self.priorInfo);
        end

        function n = getNumAlternatives(self)
            % N = GETNUMALTERNATIVES()
            % Return number of design points (1 if no Prior)
            if isempty(self.design)
                self.buildDesign();
            end
            n = length(self.design);
        end

        function probs = getProbabilities(self)
            % PROBS = GETPROBABILITIES()
            % Return vector of design-point weights
            if isempty(self.design)
                self.buildDesign();
            end
            probs = [self.design.weight];
        end

        function E = getNumberOfModels(self)
            % E = GETNUMBEROFMODELS()
            % Return number of ensemble models (design points)
            %
            % Overrides EnsembleSolver to return count based on the design,
            % since ensemble is not populated until init().
            if ~isempty(self.ensemble)
                E = length(self.ensemble);
            else
                E = self.getNumAlternatives();
            end
        end

        %% EnsembleSolver abstract method implementations

        function init(self)
            % INIT Initialize the UQ solver
            %
            % Expands the model into a family of concrete models,
            % one for each Prior alternative.

            if isempty(self.design)
                self.buildDesign();
            end
            line_debug('UQ init: expanding model into %d alternatives', length(self.design));
            if isempty(self.priorInfo)
                % No Prior - just use the original model
                self.ensemble = {self.model};
                self.solvers{1} = self.solverFactory(self.model);
                return;
            end

            n = length(self.design);
            L = length(self.priorInfo);
            self.ensemble = cell(1, n);
            self.solvers = cell(1, n);

            for i = 1:n
                % Deep copy the model
                modelCopy = self.originalModel.copy();
                % setService does NOT invalidate a cached sn (Queue.setService:
                % deliberate, for SolverLN's iteration cost), and copy() carries
                % the cache, so without this every design point would be solved
                % with the Prior's MIXTURE moments -- see _kb/06-solver-catalog.md
                modelCopy.resetStruct();

                nodes = modelCopy.getNodes();
                classes = modelCopy.getClasses();

                % Replace every Prior with its concrete alternative at this
                % design point
                for l = 1:L
                    node = nodes{self.priorInfo(l).nodeIdx};
                    class = classes{self.priorInfo(l).classIdx};
                    concreteDist = self.design(i).dists{l};

                    if strcmp(self.priorInfo(l).type, 'service')
                        node.setService(class, concreteDist);
                    elseif strcmp(self.priorInfo(l).type, 'arrival')
                        node.setArrival(class, concreteDist);
                    end
                end

                % Rename model to indicate alternative
                modelCopy.name = sprintf('%s_alt%d', self.originalModel.name, i);

                self.ensemble{i} = modelCopy;
                self.solvers{i} = self.solverFactory(modelCopy);
            end
        end

        function pre(self, it)
            % PRE Pre-iteration operations (no-op for UQ)
            % UQ only needs a single iteration
        end

        function varargout = iterate(self, varargin)
            % ITERATE Run the ensemble, narrating the run on the console
            %
            % Solver console: UQ drives an ensemble of alternative models and
            % does not pass through runAnalyzerChecks. Every entry point
            % (runAnalyzer, getAvg, getAvgTable) reaches the analysis through
            % iterate, so the run is opened here. The guard must live until
            % this function returns.
            consoleGuard = LineConsole.beginRun(self, self.options); %#ok<NASGU>
            LineConsole.step('uncertainty quantification over %d alternative models', ...
                self.getNumAlternatives());
            [varargout{1:nargout}] = iterate@EnsembleSolver(self, varargin{:});
        end

        function [result, runtime] = analyze(self, it, e)
            % ANALYZE Run solver for ensemble model e
            %
            % @param it Iteration number
            % @param e Ensemble model index
            % @return result Solver result structure
            % @return runtime Solver execution time

            T0 = tic;
            solver = self.solvers{e};
            solver.runAnalyzer();
            result = solver.result;
            runtime = toc(T0);
        end

        function post(self, it)
            % POST Post-iteration operations
            %
            % Aggregates results from all ensemble models using prior weights.

            self.aggregateResults();
        end

        function finish(self)
            % FINISH Finalization (no-op for UQ)
        end

        function bool = converged(self, it)
            % CONVERGED Check convergence
            %
            % UQ converges after the first iteration (it >= 1).
            bool = (it >= 1);
        end

        function [QN, UN, RN, TN, AN, WN] = getEnsembleAvg(self)
            % GETENSEMBLEAVG Get per-model average metrics
            %
            % Returns cell arrays with metrics from each ensemble model.

            n = self.getNumberOfModels();
            QN = cell(1, n);
            UN = cell(1, n);
            RN = cell(1, n);
            TN = cell(1, n);
            AN = cell(1, n);
            WN = cell(1, n);

            for e = 1:n
                if ~isempty(self.results) && size(self.results, 2) >= e && ~isempty(self.results{1, e})
                    res = self.results{1, e};
                    if isfield(res, 'Avg')
                        if isfield(res.Avg, 'Q'), QN{e} = res.Avg.Q; end
                        if isfield(res.Avg, 'U'), UN{e} = res.Avg.U; end
                        if isfield(res.Avg, 'R'), RN{e} = res.Avg.R; end
                        if isfield(res.Avg, 'T'), TN{e} = res.Avg.T; end
                        if isfield(res.Avg, 'A'), AN{e} = res.Avg.A; end
                        if isfield(res.Avg, 'W'), WN{e} = res.Avg.W; end
                    end
                end
            end
        end

        %% Result aggregation methods

        function aggregateResults(self)
            % AGGREGATERESULTS Compute prior-weighted aggregate metrics

            if isempty(self.results)
                return;
            end

            n = self.getNumberOfModels();
            probs = self.getProbabilities();

            % Initialize aggregated result
            self.aggregatedResult = struct();
            self.aggregatedResult.solver = 'UQ';
            self.aggregatedResult.Avg = struct();

            % Get dimensions from first result
            firstResult = [];
            for e = 1:n
                if ~isempty(self.results{1, e})
                    firstResult = self.results{1, e};
                    break;
                end
            end

            if isempty(firstResult) || ~isfield(firstResult, 'Avg')
                return;
            end

            % Aggregate each metric field
            fields = {'Q', 'U', 'R', 'T', 'A', 'W', 'C', 'X'};
            for f = 1:length(fields)
                fname = fields{f};
                if isfield(firstResult.Avg, fname) && ~isempty(firstResult.Avg.(fname))
                    self.aggregatedResult.Avg.(fname) = zeros(size(firstResult.Avg.(fname)));
                    for e = 1:n
                        if ~isempty(self.results{1, e}) && ...
                           isfield(self.results{1, e}, 'Avg') && ...
                           isfield(self.results{1, e}.Avg, fname) && ...
                           ~isempty(self.results{1, e}.Avg.(fname))
                            self.aggregatedResult.Avg.(fname) = self.aggregatedResult.Avg.(fname) + ...
                                probs(e) * self.results{1, e}.Avg.(fname);
                        end
                    end
                end
            end
        end

        %% Override NetworkSolver methods for aggregated results

        function [QN, UN, RN, TN, AN, WN] = getAvg(self, varargin)
            % GETAVG Return prior-weighted average metrics
            %
            % Returns aggregated metrics weighted by prior probabilities.

            % Run solver if not done
            if isempty(self.results)
                self.iterate();
            end

            QN = []; UN = []; RN = []; TN = []; AN = []; WN = [];
            if ~isempty(self.aggregatedResult) && isfield(self.aggregatedResult, 'Avg')
                if isfield(self.aggregatedResult.Avg, 'Q'), QN = self.aggregatedResult.Avg.Q; end
                if isfield(self.aggregatedResult.Avg, 'U'), UN = self.aggregatedResult.Avg.U; end
                if isfield(self.aggregatedResult.Avg, 'R'), RN = self.aggregatedResult.Avg.R; end
                if isfield(self.aggregatedResult.Avg, 'T'), TN = self.aggregatedResult.Avg.T; end
                if isfield(self.aggregatedResult.Avg, 'A'), AN = self.aggregatedResult.Avg.A; end
                if isfield(self.aggregatedResult.Avg, 'W'), WN = self.aggregatedResult.Avg.W; end
            end
        end

        function varargout = getAvgTable(self, varargin)
            % GETAVGTABLE Return prior-weighted average table
            %
            % Returns a table of aggregated metrics weighted by prior probabilities.
            % The result recorder captures the returned table together with the solver
            % that produced it -- see LineResultRecorder. Recording an ensemble here
            % rather than in the member solver it delegates to is what keeps an
            % AUTO/LN/ENV/UQ answer from being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getAvgTable_impl(varargin{:});
            LineResultRecorder.capture(scope, self, 'avg', varargout{1});
        end

        function AvgTable = getAvgTable_impl(self, varargin)
            % GETAVGTABLE_IMPL Implementation of GETAVGTABLE; see the wrapper above.

            % Run solver if not done
            if isempty(self.results)
                self.iterate();
            end

            [QN, UN, RN, TN, AN, WN] = self.getAvg();

            % Get model structure
            sn = self.originalModel.getStruct(false);
            M = sn.nstations;
            K = sn.nclasses;

            % Build table data
            Station = {};
            JobClass = {};
            QLen = [];
            Util = [];
            RespT = [];
            ResidT = [];
            ArvR = [];
            Tput = [];

            for ist = 1:M
                for k = 1:K
                    Q_val = 0; U_val = 0; R_val = 0; T_val = 0; A_val = 0; W_val = 0;

                    if ~isempty(QN), Q_val = QN(ist, k); end
                    if ~isempty(UN), U_val = UN(ist, k); end
                    if ~isempty(RN), R_val = RN(ist, k); end
                    if ~isempty(TN), T_val = TN(ist, k); end
                    if ~isempty(AN), A_val = AN(ist, k); end
                    if ~isempty(WN), W_val = WN(ist, k); end

                    % Only include rows with non-zero values
                    if Q_val > 0 || U_val > 0 || T_val > 0
                        Station{end+1, 1} = sn.nodenames{sn.stationToNode(ist)};
                        JobClass{end+1, 1} = sn.classnames{k};
                        QLen(end+1, 1) = Q_val;
                        Util(end+1, 1) = U_val;
                        RespT(end+1, 1) = R_val;
                        ResidT(end+1, 1) = W_val;
                        ArvR(end+1, 1) = A_val;
                        Tput(end+1, 1) = T_val;
                    end
                end
            end

            if isempty(Station)
                AvgTable = table();
                return;
            end

            Station = categorical(Station);
            JobClass = categorical(JobClass);

            AvgTable = table(Station, JobClass, QLen, Util, RespT, ResidT, ArvR, Tput);
        end

        %% UQ-specific result methods

        function ptable = getPosteriorTable(self)
            % GETPOSTERIORTABLE Return table with per-alternative results
            %
            % Returns a table showing metrics for each Prior alternative
            % along with its probability.

            % Run solver if not done
            if isempty(self.results)
                self.iterate();
            end

            probs = self.getProbabilities();
            n = self.getNumberOfModels();
            sn = self.originalModel.getStruct(false);
            M = sn.nstations;
            K = sn.nclasses;

            % Build table data
            Alternative = [];
            Probability = [];
            Station = {};
            JobClass = {};
            QLen = [];
            Util = [];
            RespT = [];
            Tput = [];

            for alt = 1:n
                if isempty(self.results{1, alt})
                    continue;
                end
                res = self.results{1, alt};
                if ~isfield(res, 'Avg')
                    continue;
                end

                for ist = 1:M
                    for k = 1:K
                        Q_val = 0; U_val = 0; R_val = 0; T_val = 0;

                        if isfield(res.Avg, 'Q') && ~isempty(res.Avg.Q)
                            Q_val = res.Avg.Q(ist, k);
                        end
                        if isfield(res.Avg, 'U') && ~isempty(res.Avg.U)
                            U_val = res.Avg.U(ist, k);
                        end
                        if isfield(res.Avg, 'R') && ~isempty(res.Avg.R)
                            R_val = res.Avg.R(ist, k);
                        end
                        if isfield(res.Avg, 'T') && ~isempty(res.Avg.T)
                            T_val = res.Avg.T(ist, k);
                        end

                        % Only include rows with non-zero values
                        if Q_val > 0 || U_val > 0 || T_val > 0
                            Alternative(end+1, 1) = alt;
                            Probability(end+1, 1) = probs(alt);
                            Station{end+1, 1} = sn.nodenames{sn.stationToNode(ist)};
                            JobClass{end+1, 1} = sn.classnames{k};
                            QLen(end+1, 1) = Q_val;
                            Util(end+1, 1) = U_val;
                            RespT(end+1, 1) = R_val;
                            Tput(end+1, 1) = T_val;
                        end
                    end
                end
            end

            if isempty(Alternative)
                ptable = table();
                return;
            end

            Station = categorical(Station);
            JobClass = categorical(JobClass);

            ptable = table(Alternative, Probability, Station, JobClass, QLen, Util, RespT, Tput);
        end

        function empDist = getPosteriorDist(self, metric, station, class)
            % GETPOSTERIORDIST Return empirical distribution of a metric
            %
            % Returns an EmpiricalCDF object representing the posterior
            % distribution of the specified metric across Prior alternatives.
            %
            % @param metric Metric name: 'Q', 'U', 'R', 'T', 'A', 'W'
            % @param station Station node or station index
            % @param class JobClass object or class index
            % @return empDist EmpiricalCDF object

            % Run solver if not done
            if isempty(self.results)
                self.iterate();
            end

            probs = self.getProbabilities();
            n = self.getNumberOfModels();

            % Get station index
            if isnumeric(station)
                ist = station;
            elseif isa(station, 'Station')
                ist = station.stationIndex;
            else
                line_error(mfilename, 'station must be a Station object or numeric index');
            end

            % Get class index
            if isnumeric(class)
                k = class;
            elseif isa(class, 'JobClass')
                k = class.index;
            else
                line_error(mfilename, 'class must be a JobClass object or numeric index');
            end

            % Extract metric values for each alternative
            values = zeros(n, 1);
            for e = 1:n
                if ~isempty(self.results{1, e}) && isfield(self.results{1, e}, 'Avg')
                    res = self.results{1, e}.Avg;
                    if isfield(res, metric) && ~isempty(res.(metric))
                        values(e) = res.(metric)(ist, k);
                    else
                        values(e) = NaN;
                    end
                else
                    values(e) = NaN;
                end
            end

            % Build empirical CDF
            % Sort values and corresponding probabilities
            [sortedVals, sortIdx] = sort(values);
            sortedProbs = probs(sortIdx);
            cdfVals = cumsum(sortedProbs);

            % Create data matrix for EmpiricalCDF: [CDF_value, X_value]
            cdfData = [cdfVals(:), sortedVals(:)];

            % Create EmpiricalCDF object
            empDist = EmpiricalCDF(cdfData);
        end

        %% Interval (support-only) uncertainty

        function ival = getInterval(self)
            % IVAL = GETINTERVAL()
            % Range of every metric over the support of the Priors.
            %
            % This drops the weights and keeps only the endpoints, which is
            % the epistemic case in which the modeller can bound a parameter
            % but not distribute it. Two regimes, distinguished by ival.exact:
            %
            %  exact = true   The model is a single-class closed product-form
            %                 network with load-independent single-server
            %                 queues and delays, so pfqn_mva_interval returns
            %                 the exact hull of MVA over the whole demand box
            %                 by the monotonicity of Luthi and Haring (1998).
            %                 No ensemble run is needed and the interval is
            %                 attained, not sampled.
            %  exact = false  Fallback: the range across the design points
            %                 that were actually solved. For a discrete Prior
            %                 this is again exact, because the design visits
            %                 the whole support; for a continuous one it is an
            %                 INNER approximation of the true range, since a
            %                 quadrature node is not an endpoint of the
            %                 support. It is therefore not an enclosure.
            %
            % The interval is conditional on the true parameters lying inside
            % the Prior supports. It is not a bound on the exact solution of
            % the network and must not be composed with SolverBA brackets.
            %
            % @return ival Struct with fields Q, U, R, T, W of size
            %              nstations x nclasses x 2, the trailing index
            %              selecting the lower and the upper endpoint, plus X
            %              and Rtot (1 x 2, exact path only), exact (logical)
            %              and method (char).

            [ok, why] = self.qualifiesForIntervalMVA();
            if ok
                ival = self.intervalByMVA();
            else
                line_warning(mfilename, ...
                    ['Exact interval MVA does not apply (%s); reporting the range over the design ', ...
                    'points, which spans the sampled support only and is not an enclosure.'], why);
                ival = self.intervalBySampling();
            end
        end

        function itable = getIntervalTable(self)
            % ITABLE = GETINTERVALTABLE()
            % Tabular form of getInterval, two columns per metric.

            ival = self.getInterval();
            sn = self.originalModel.getStruct(false);
            M = sn.nstations;
            K = sn.nclasses;

            Station = {}; JobClass = {};
            QLen_lo = []; QLen_up = []; Util_lo = []; Util_up = [];
            RespT_lo = []; RespT_up = []; Tput_lo = []; Tput_up = [];

            for ist = 1:M
                for k = 1:K
                    if ival.Q(ist,k,2) <= 0 && ival.U(ist,k,2) <= 0 && ival.T(ist,k,2) <= 0
                        continue
                    end
                    Station{end+1, 1} = sn.nodenames{sn.stationToNode(ist)}; %#ok<AGROW>
                    JobClass{end+1, 1} = sn.classnames{k}; %#ok<AGROW>
                    QLen_lo(end+1, 1) = ival.Q(ist,k,1); %#ok<AGROW>
                    QLen_up(end+1, 1) = ival.Q(ist,k,2); %#ok<AGROW>
                    Util_lo(end+1, 1) = ival.U(ist,k,1); %#ok<AGROW>
                    Util_up(end+1, 1) = ival.U(ist,k,2); %#ok<AGROW>
                    RespT_lo(end+1, 1) = ival.R(ist,k,1); %#ok<AGROW>
                    RespT_up(end+1, 1) = ival.R(ist,k,2); %#ok<AGROW>
                    Tput_lo(end+1, 1) = ival.T(ist,k,1); %#ok<AGROW>
                    Tput_up(end+1, 1) = ival.T(ist,k,2); %#ok<AGROW>
                end
            end

            if isempty(Station)
                itable = table();
                return;
            end

            Station = categorical(Station);
            JobClass = categorical(JobClass);
            itable = table(Station, JobClass, QLen_lo, QLen_up, Util_lo, Util_up, ...
                RespT_lo, RespT_up, Tput_lo, Tput_up);
        end

        function [ok, why] = qualifiesForIntervalMVA(self)
            % [OK, WHY] = QUALIFIESFORINTERVALMVA()
            % Whether the monotonicity theorems behind pfqn_mva_interval hold
            % for this model. WHY names the first violated condition.

            ok = false;
            sn = self.originalModel.getStruct(false);
            if ~isempty(self.priorInfo) && ~all(strcmp({self.priorInfo.type}, 'service'))
                why = 'a Prior sits on an arrival process, so the model is open'; return
            end
            if sn.nclasses ~= 1
                why = 'the theorems are proved for a single class only'; return
            end
            if sn.nclosedjobs <= 0
                why = 'the class is not closed'; return
            end
            if sn.nnodes ~= sn.nstations
                why = 'the model has nodes that are not stations'; return
            end
            isinf_ = (sn.sched == SchedStrategy.INF);
            if any(sn.nservers(~isinf_) > 1)
                why = 'a queueing station has more than one server'; return
            end
            if ~all(isinf_ | sn.sched == SchedStrategy.PS | sn.sched == SchedStrategy.FCFS)
                why = 'a station is neither delay, PS nor FCFS'; return
            end
            stations = arrayfun(@(p) sn.nodeToStation(p.nodeIdx), self.priorInfo);
            if numel(unique(stations)) ~= numel(stations)
                why = 'two Priors sit on the same station'; return
            end
            ok = true;
            why = '';
        end

        function ival = intervalByMVA(self)
            % IVAL = INTERVALBYMVA()
            % Exact hull through pfqn_mva_interval. The demand box is the
            % nominal demand vector with the prior-carrying stations widened
            % to the range of mean service times over the Prior support.

            sn = self.originalModel.getStruct(false);
            M = sn.nstations;
            V = sn.visits{1}(:);
            isinf_ = (sn.sched == SchedStrategy.INF);
            ST = 1 ./ sn.rates(:,1);
            ST(isnan(ST)) = 0;
            STlo = ST; STup = ST;

            for l = 1:numel(self.priorInfo)
                ist = sn.nodeToStation(self.priorInfo(l).nodeIdx);
                [lo, up] = UQ.priorMeanRange(self.priorInfo(l).prior, self.getUQNodes());
                STlo(ist) = lo; STup(ist) = up;
            end

            Dlo = V .* STlo; Dup = V .* STup;
            Zint = [sum(Dlo(isinf_)), sum(Dup(isinf_))];
            qidx = find(~isinf_);
            [X, Qq, Uq, Rq, Rtot] = pfqn_mva_interval([Dlo(qidx), Dup(qidx)], sn.nclosedjobs, Zint);

            Q = zeros(M,2); U = zeros(M,2); R = zeros(M,2); W = zeros(M,2);
            Q(qidx,:) = Qq;
            U(qidx,:) = Uq;
            W(qidx,:) = Rq;
            R(qidx,:) = Rq ./ V(qidx);
            % A delay station never queues, so its residence time is its own
            % demand interval and its population is the throughput times that
            % demand, enclosed as a product of two intervals.
            W(isinf_,:) = [Dlo(isinf_), Dup(isinf_)];
            R(isinf_,:) = [STlo(isinf_), STup(isinf_)];
            Q(isinf_,:) = [X(1)*Dlo(isinf_), X(2)*Dup(isinf_)];
            U(isinf_,:) = Q(isinf_,:);

            ival = struct();
            ival.X = X;
            ival.Rtot = Rtot;
            ival.Q = reshape(Q, [M,1,2]);
            ival.U = reshape(U, [M,1,2]);
            ival.R = reshape(R, [M,1,2]);
            ival.W = reshape(W, [M,1,2]);
            ival.T = reshape(V * X, [M,1,2]);
            ival.exact = true;
            ival.method = 'mvainterval';
        end

        function ival = intervalBySampling(self)
            % IVAL = INTERVALBYSAMPLING()
            % Range of each metric across the design points that were solved.

            if isempty(self.results)
                self.iterate();
            end
            n = self.getNumberOfModels();
            fields = {'Q', 'U', 'R', 'T', 'W'};
            ival = struct();
            for f = 1:length(fields)
                fname = fields{f};
                lo = []; up = [];
                for e = 1:n
                    res = self.results{1, e};
                    if isempty(res) || ~isfield(res, 'Avg') || ~isfield(res.Avg, fname) || isempty(res.Avg.(fname))
                        continue
                    end
                    v = res.Avg.(fname);
                    if isempty(lo)
                        lo = v; up = v;
                    else
                        lo = min(lo, v); up = max(up, v);
                    end
                end
                if isempty(lo)
                    ival.(fname) = [];
                else
                    ival.(fname) = cat(3, lo, up);
                end
            end
            ival.X = [];
            ival.Rtot = [];
            ival.exact = false;
            ival.method = 'sampled';
        end

        %% Required Solver abstract methods

        function runtime = runAnalyzer(self, options)
            % RUNANALYZER Run the UQ analysis
            %
            % @param options Solver options (optional)
            % @return runtime Total runtime in seconds

            T0 = tic;
            if nargin >= 2 && ~isempty(options)
                self.setOptions(options);
            end
            line_debug('UQ solver starting: nalternatives=%d', self.getNumAlternatives());
            line_debug('Default method: using UQ ensemble analysis\n');
            self.iterate();
            runtime = toc(T0);
        end

        function sn = getStruct(self)
            % GETSTRUCT Return model structure
            %
            % Returns the structure of the original model.
            sn = self.originalModel.getStruct(false);
        end

        %% Static methods

        function [allMethods] = listValidMethods(self)
            % LISTVALIDMETHODS Return valid methods
            allMethods = {'default', 'discrete', 'quadrature', 'montecarlo'};
        end

        %% Uncertainty summaries

        function [m, v] = getMoments(self, metric, station, class)
            % [M, V] = GETMOMENTS(METRIC, STATION, CLASS)
            % Weighted mean and variance of a metric over the design.
            %
            % The mean is the unconditional expectation of Trivedi and Bobbio
            % (2017), Eq. (3.68); the variance is the second moment of the
            % same weighting, as derived for the cold-standby case in their
            % Sec. 8.5.1. Both are exact for a discrete Prior and quadrature-
            % or sample-approximate for a continuous one.
            %
            % @param metric Metric name ('Q','U','R','T','A','W')
            % @param station Station object or index
            % @param class Class object or index
            % @return m Weighted mean
            % @return v Weighted variance

            [vals, w] = self.getSamples(metric, station, class);
            m = sum(w .* vals);
            v = sum(w .* (vals - m).^2);
        end

        function ci = getCredibleInterval(self, metric, station, class, level)
            % CI = GETCREDIBLEINTERVAL(METRIC, STATION, CLASS, LEVEL)
            % Equal-tailed credible interval from the weighted empirical CDF.
            %
            % @param metric Metric name ('Q','U','R','T','A','W')
            % @param station Station object or index
            % @param class Class object or index
            % @param level Coverage level in (0,1), default 0.95
            % @return ci Two-element vector [lower, upper]

            if nargin < 5 || isempty(level)
                level = 0.95;
            end
            if level <= 0 || level >= 1
                line_error(mfilename, 'level must lie strictly between 0 and 1');
            end

            [vals, w] = self.getSamples(metric, station, class);
            [vals, ord] = sort(vals);
            w = w(ord);
            cw = cumsum(w) / sum(w);

            alpha = (1 - level) / 2;
            lo = vals(find(cw >= alpha, 1, 'first'));
            hi = vals(find(cw >= 1 - alpha, 1, 'first'));
            if isempty(lo), lo = vals(1); end
            if isempty(hi), hi = vals(end); end
            ci = [lo, hi];
        end

        function [vals, w] = getSamples(self, metric, station, class)
            % [VALS, W] = GETSAMPLES(METRIC, STATION, CLASS)
            % Per-design-point metric values and their weights.

            if isempty(self.results)
                self.iterate();
            end

            ist = self.resolveStationIndex(station);
            k = self.resolveClassIndex(class);

            n = self.getNumberOfModels();
            w = self.getProbabilities();
            vals = zeros(1, n);
            for e = 1:n
                res = self.results{1, e};
                if isempty(res) || ~isfield(res, 'Avg') || ~isfield(res.Avg, metric) || isempty(res.Avg.(metric))
                    line_error(mfilename, sprintf('Metric %s unavailable for design point %d', metric, e));
                end
                vals(e) = res.Avg.(metric)(ist, k);
            end
        end

        function ist = resolveStationIndex(self, station)
            % IST = RESOLVESTATIONINDEX(STATION)
            % Accept a Station object or a numeric index.
            if isnumeric(station)
                ist = station;
            elseif isa(station, 'Station')
                ist = station.stationIndex;
            else
                line_error(mfilename, 'station must be a Station object or numeric index');
            end
        end

        function k = resolveClassIndex(self, class)
            % K = RESOLVECLASSINDEX(CLASS)
            % Accept a JobClass object or a numeric index.
            if isnumeric(class)
                k = class;
            elseif isa(class, 'JobClass')
                k = class.index;
            else
                line_error(mfilename, 'class must be a JobClass object or numeric index');
            end
        end
    end

    methods (Static)
        function [lo, up] = priorMeanRange(prior, n)
            % [LO, UP] = PRIORMEANRANGE(PRIOR, N)
            % Range of the mean of a Prior over its alternatives.
            %
            % Exact for a discrete Prior, whose alternatives are the support.
            % A continuous Prior is first discretized on N quadrature nodes,
            % so the range is that of the discretized support: an unbounded
            % parameter density is never reached at its tails.

            dists = prior.discretize(n, 'quadrature');
            m = zeros(1, length(dists));
            for i = 1:length(dists)
                m(i) = dists{i}.getMean();
            end
            lo = min(m); up = max(m);
        end

        function idx = unrankIndex(i, counts)
            % IDX = UNRANKINDEX(I, COUNTS)
            % Map the linear index I in 1..prod(COUNTS) to a subscript vector
            % over a mixed-radix grid, first coordinate varying fastest.
            L = length(counts);
            idx = zeros(1, L);
            rem = i - 1;
            for l = 1:L
                idx(l) = mod(rem, counts(l)) + 1;
                rem = floor(rem / counts(l));
            end
        end

        function method = resolveUQMethod(name)
            % METHOD = RESOLVEUQMETHOD(NAME)
            % The design a method name selects: 'default', 'discrete' and
            % 'quadrature' are the tensor-product (quadrature) design, and
            % 'montecarlo' the sampled one. One table for getUQMethod and the
            % gate, so the two cannot disagree on a name.
            method = 'quadrature';
            if isempty(name)
                return
            end
            switch char(name)
                case {'default', 'discrete', 'quadrature'}
                    method = 'quadrature';
                case 'montecarlo'
                    method = 'montecarlo';
                otherwise
                    line_error(mfilename, sprintf('Unknown UQ method: %s', char(name)));
            end
        end

        function featSupported = getFeatureSet()
            % GETFEATURESET Return supported features
            featSupported = SolverFeatureSet;
            featSupported.setTrue('Prior');
        end

        function [bool, featSupported] = supports(model)
            % SUPPORTS Check if model is supported
            %
            % A MODEL WITH NO UNCERTAIN PARAMETER IS NOT A UQ MODEL. This
            % returned true unconditionally, which made UQ claim every model in
            % the language: SolverAUTO.listValidMethods offered every 'uq.*'
            % method name on an ordinary network, whose posterior is a single design
            % point equal to the point estimate the inner solver already
            % returns. The JAR (SolverUQ.supports -> detectPrior() != null) and
            % native python (hasPriorDistribution) have always gated it this
            % way, so this closes an m/j/p divergence rather than tightening a
            % rule the other two share.
            %
            % The feature set stays the answer to a DIFFERENT question and is
            % returned unchanged: what UQ itself consumes is the Prior, while
            % every other feature is the INNER solver's to accept or refuse, on
            % a model from which the Prior has already been removed.
            featSupported = UQ.getFeatureSet();
            bool = UQ.modelHasPrior(model);
        end

        function tf = modelHasPrior(model)
            % TF = MODELHASPRIOR(MODEL)
            % True when the model carries an uncertain parameter, i.e. when a
            % service or arrival process is a Prior and the model therefore
            % expands into one instance per design point. SolverAUTO asks the
            % same question through this method, so the constructor's routing
            % and this gate cannot disagree.
            tf = false;
            if ~isa(model,'Network')
                return
            end
            nodes = model.getNodes();
            classes = model.getClasses();
            for i = 1:length(nodes)
                node = nodes{i};
                if isa(node,'Queue') || isa(node,'Delay')
                    for c = 1:length(classes)
                        try
                            dist = node.getService(classes{c});
                        catch
                            dist = [];
                        end
                        if ~isempty(dist) && isa(dist,'Prior')
                            tf = true;
                            return
                        end
                    end
                elseif isa(node,'Source')
                    for c = 1:length(classes)
                        if size(node.arrivalProcess,2) >= c && ~isempty(node.arrivalProcess{1,c}) ...
                                && isa(node.arrivalProcess{1,c},'Prior')
                            tf = true;
                            return
                        end
                    end
                end
            end
        end

        function options = defaultOptions()
            % DEFAULTOPTIONS Return default options
            options = Solver.defaultOptions();
            options.iter_max = 1;  % Single iteration for UQ
            % Nodes per continuous Prior. Each node is a full solver run, so
            % the simulation-oriented default of Solver.defaultOptions is far
            % too large here.
            options.samples = 11;
        end
    end
end
