classdef Posterior < EnsembleSolver
    % Posterior Solver wrapper for models with Prior distributions
    %
    % Posterior detects Prior distributions in a model, expands the model
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
    % post = Posterior(model, @SolverMVA);
    % avgTable = post.getAvgTable();            % Prior-weighted expectations
    % postTable = post.getPosteriorTable();     % Per-alternative results
    % postDist = post.getPosteriorDist('R', queue, class); % Response time distribution
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        solverFactory;      % Function handle to create solvers: @(model) SolverXXX(model)
        priorInfo;          % Structure with Prior detection info
        aggregatedResult;   % Prior-weighted aggregate metrics
        originalModel;      % Reference to original model with Prior
    end

    methods
        function self = Posterior(model, solverFactory, varargin)
            % POSTERIOR Create a Posterior solver wrapper
            %
            % @brief Creates a Posterior wrapper for uncertainty analysis
            % @param model Network model (may contain Prior distributions)
            % @param solverFactory Function handle: @(m) SolverXXX(m) or solver class name
            % @param varargin Optional solver options
            % @return self Posterior instance

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
                self.setOptions(Solver.parseOptions(varargin, Posterior.defaultOptions));
            else
                self.setOptions(Posterior.defaultOptions);
            end

            self.priorInfo = [];
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

            % Validate: currently only single Prior supported
            if length(priors) > 1
                line_error(mfilename, 'Currently only a single Prior distribution per model is supported');
            end

            if isempty(priors)
                self.priorInfo = [];
            else
                self.priorInfo = priors{1};
            end
        end

        function hasPrior = hasPriorDistribution(self)
            % HASPRIOR = HASPRIORDISTRIBUTION()
            % Return true if model contains a Prior distribution
            hasPrior = ~isempty(self.priorInfo);
        end

        function n = getNumAlternatives(self)
            % N = GETNUMALTERNATIVES()
            % Return number of Prior alternatives (1 if no Prior)
            if isempty(self.priorInfo)
                n = 1;
            else
                n = self.priorInfo.prior.getNumAlternatives();
            end
        end

        function probs = getProbabilities(self)
            % PROBS = GETPROBABILITIES()
            % Return vector of prior probabilities
            if isempty(self.priorInfo)
                probs = 1;
            else
                probs = self.priorInfo.prior.probabilities;
            end
        end

        function E = getNumberOfModels(self)
            % E = GETNUMBEROFMODELS()
            % Return number of ensemble models (alternatives)
            %
            % Overrides EnsembleSolver to return count based on Prior,
            % since ensemble is not populated until init().
            if ~isempty(self.ensemble)
                E = length(self.ensemble);
            else
                E = self.getNumAlternatives();
            end
        end

        %% EnsembleSolver abstract method implementations

        function init(self)
            % INIT Initialize the Posterior solver
            %
            % Expands the model into a family of concrete models,
            % one for each Prior alternative.

            line_debug('Posterior init: expanding model into %d alternatives', self.getNumAlternatives());
            if isempty(self.priorInfo)
                % No Prior - just use the original model
                self.ensemble = {self.model};
                self.solvers{1} = self.solverFactory(self.model);
                return;
            end

            prior = self.priorInfo.prior;
            n = prior.getNumAlternatives();
            self.ensemble = cell(1, n);
            self.solvers = cell(1, n);

            for i = 1:n
                % Deep copy the model
                modelCopy = self.originalModel.copy();

                % Get the corresponding node and class in the copy
                nodes = modelCopy.getNodes();
                classes = modelCopy.getClasses();
                node = nodes{self.priorInfo.nodeIdx};
                class = classes{self.priorInfo.classIdx};

                % Replace Prior with concrete alternative
                concreteDist = prior.getAlternative(i);

                if strcmp(self.priorInfo.type, 'service')
                    node.setService(class, concreteDist);
                elseif strcmp(self.priorInfo.type, 'arrival')
                    node.setArrival(class, concreteDist);
                end

                % Rename model to indicate alternative
                modelCopy.name = sprintf('%s_alt%d', self.originalModel.name, i);

                self.ensemble{i} = modelCopy;
                self.solvers{i} = self.solverFactory(modelCopy);
            end
        end

        function pre(self, it)
            % PRE Pre-iteration operations (no-op for Posterior)
            % Posterior only needs a single iteration
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
            % FINISH Finalization (no-op for Posterior)
        end

        function bool = converged(self, it)
            % CONVERGED Check convergence
            %
            % Posterior converges after the first iteration (it >= 1).
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
            self.aggregatedResult.solver = 'Posterior';
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

        function AvgTable = getAvgTable(self, varargin)
            % GETAVGTABLE Return prior-weighted average table
            %
            % Returns a table of aggregated metrics weighted by prior probabilities.

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

        %% Posterior-specific result methods

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

        %% Required Solver abstract methods

        function runtime = runAnalyzer(self, options)
            % RUNANALYZER Run the Posterior analysis
            %
            % @param options Solver options (optional)
            % @return runtime Total runtime in seconds

            T0 = tic;
            if nargin >= 2 && ~isempty(options)
                self.setOptions(options);
            end
            line_debug('Posterior solver starting: nalternatives=%d', self.getNumAlternatives());
            line_debug('Default method: using Posterior ensemble analysis\n');
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
            allMethods = {'default'};
        end
    end

    methods (Static)
        function featSupported = getFeatureSet()
            % GETFEATURESET Return supported features
            featSupported = SolverFeatureSet;
            featSupported.setTrue('Prior');
        end

        function [bool, featSupported] = supports(model)
            % SUPPORTS Check if model is supported
            bool = true;
            featSupported = Posterior.getFeatureSet();
        end

        function options = defaultOptions()
            % DEFAULTOPTIONS Return default options
            options = Solver.defaultOptions();
            options.iter_max = 1;  % Single iteration for Posterior
        end
    end
end
