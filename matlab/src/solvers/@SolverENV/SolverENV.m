classdef SolverENV < EnsembleSolver
    % ENV - Ensemble environment solver for models with random environment changes
    %
    % SolverENV analyzes queueing networks operating in random environments where
    % system parameters (arrival rates, service rates, routing) change according
    % to an underlying environmental process. It solves ensemble models by analyzing
    % each environmental stage and computing environment-averaged performance metrics.
    %
    % @brief Environment solver for networks in random changing environments
    %
    % Key characteristics:
    % - Random environment with multiple operational stages
    % - Environmental process governing parameter changes
    % - Ensemble model analysis across environment stages
    % - Environment-averaged performance computation
    % - Stage-dependent system behavior modeling
    %
    % Environment solver features:
    % - Multi-stage environmental modeling
    % - Stage transition matrix analysis
    % - Weighted performance metric computation
    % - Environmental ensemble solution
    % - Adaptive parameter modeling
    %
    % SolverENV is ideal for:
    % - Systems with time-varying parameters
    % - Networks subject to environmental fluctuations
    % - Multi-mode operational system analysis
    % - Performance under uncertainty modeling
    % - Adaptive system behavior analysis
    %
    % Example:
    % @code
    % env_model = Environment(stages, transitions); % Define environment
    % solver = SolverENV(env_model, @SolverMVA, options);
    % metrics = solver.getEnsembleAvg(); % Environment-averaged metrics
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.


    properties
        env; % user-supplied representation of each stage transition
        envObj;
        sn;
        resetFromMarginal;
        resetEnvRates; % function implementing the reset policy for environment rates
        stateDepMethod = ''; % state-dependent method configuration
        SMPMethod = false;  % Use DTMC-based computation for Semi-Markov Processes
        % Enhanced init properties (aligned with JAR)
        ServerNum;   % Cell array of server counts per class
        SRates;      % Cell array of service rates per class
        E0;          % Rate matrix
        Eutil;       % Infinitesimal generator
        transitionCdfs;  % Transition CDF functions
        sojournCdfs;     % Sojourn time CDF functions
        dtmcP;       % DTMC transition matrix
        holdTimeMatrix;  % Hold time matrix
        newMethod = false;  % Use DTMC-based computation
        compression = false;  % Use Courtois decomposition
        compressionResult;  % Results from Courtois decomposition
        Ecompress;  % Number of macro-states after compression
        MS;  % Macro-state partition
    end

    methods
        function self = SolverENV(renv, solverFactory, options)
            % SELF = SOLVERENV(ENV,SOLVERFACTORY,OPTIONS)
            self@EnsembleSolver(renv, mfilename);
            if nargin>=3 %exist('options','var')
                self.setOptions(options);
            else
                self.setOptions(SolverENV.defaultOptions);
            end

            % Enable SMP method if specified in options
            if isfield(self.options, 'method') && strcmpi(self.options.method, 'smp')
                self.SMPMethod = true;
                line_debug('ENV solver: SMP method enabled via options.method=''smp''');
            end

            self.envObj = renv;
            self.ensemble = renv.getEnsemble;
            env = renv.getEnv;
            self.env = env;
            for e=1:length(self.env)
                self.sn{e} = self.ensemble{e}.getStruct;
                self.setSolver(solverFactory(self.ensemble{e}),e);
            end

            for e=1:length(self.env)
                for h=1:length(self.env)
                    self.resetFromMarginal{e,h} = renv.resetFun{e,h};
                end
            end

            for e=1:length(self.env)
                for h=1:length(self.env)
                    self.resetEnvRates{e,h} = renv.resetEnvRatesFun{e,h};
                end
            end

            for e=1:length(self.env)
                for h=1:length(self.env)
                    if isa(self.env{e,h},'Disabled')
                        self.env{e,h} = Exp(0);
                    elseif ~isa(self.env{e,h},'Markovian') && ~self.SMPMethod
                        line_error(mfilename,sprintf('The distribution of the environment transition from stage %d to %d is not supported by the %s solver. Use method=''smp'' for non-Markovian distributions.',e,h,self.getName));
                    end
                end
            end

            for e=1:length(self.ensemble)
                if ~self.solvers{e}.supports(self.ensemble{e})
                    line_error(mfilename,sprintf('Model in the environment stage %d is not supported by the %s solver.',e,self.getName));
                end
            end
        end

        function setStateDepMethod(self, method)
            % SETSTATEDEPMETHOD(METHOD) Sets the state-dependent method
            if isempty(method)
                line_error(mfilename, 'State-dependent method cannot be null or empty.');
            end
            self.stateDepMethod = method;
        end

        function setNewMethod(self, flag)
            % SETNEWMETHOD(FLAG) Enable/disable DTMC-based computation
            self.newMethod = flag;
        end

        function setCompression(self, flag)
            % SETCOMPRESSION(FLAG) Enable/disable Courtois decomposition
            self.compression = flag;
        end

        function [p, eps, epsMax, q] = ctmc_decompose(self, Q, MS)
            % CTMC_DECOMPOSE Perform CTMC decomposition using configured method
            % [p, eps, epsMax, q] = CTMC_DECOMPOSE(Q, MS)
            %
            % Uses options.config.decomp to select the decomposition algorithm:
            %   'courtois'  - Courtois decomposition (default)
            %   'kms'       - Koury-McAllister-Stewart method
            %   'takahashi' - Takahashi's method
            %   'multi'     - Multigrid method (requires MSS)
            %
            % Returns:
            %   p      - steady-state probability vector
            %   eps    - NCD index
            %   epsMax - max acceptable eps value
            %   q      - randomization coefficient

            % Get decomposition/aggregation method from options
            if isfield(self.options, 'config') && isfield(self.options.config, 'da')
                method = self.options.config.da;
            else
                method = 'courtois';
            end

            % Get numsteps for iterative methods
            if isfield(self.options, 'config') && isfield(self.options.config, 'da_iter')
                numsteps = self.options.config.da_iter;
            else
                numsteps = 10;
            end

            switch lower(method)
                case 'courtois'
                    [p, ~, ~, eps, epsMax, ~, ~, ~, q] = ctmc_courtois(Q, MS);
                case 'kms'
                    [p, ~, ~, eps, epsMax] = ctmc_kms(Q, MS, numsteps);
                    q = 1.05 * max(max(abs(Q)));
                case 'takahashi'
                    [p, ~, ~, ~, eps, epsMax] = ctmc_takahashi(Q, MS, numsteps);
                    q = 1.05 * max(max(abs(Q)));
                case 'multi'
                    % Multi requires MSS (macro-macro-states), default to singletons
                    nMacro = size(MS, 1);
                    MSS = cell(nMacro, 1);
                    for i = 1:nMacro
                        MSS{i} = i;
                    end
                    [p, ~, ~, ~, eps, epsMax] = ctmc_multi(Q, MS, MSS);
                    q = 1.05 * max(max(abs(Q)));
                otherwise
                    line_error(mfilename, sprintf('Unknown decomposition method: %s', method));
            end
        end

        function bool = converged(self, it) % convergence test at iteration it
            % BOOL = CONVERGED(IT) % CONVERGENCE TEST AT ITERATION IT
            % Computes max relative absolute difference of queue lengths between iterations
            % Aligned with JAR SolverEnv.converged() implementation

            bool = false;
            if it <= 1
                return
            end

            E = self.getNumberOfModels;
            M = size(self.results{it,1}.Tran.Avg.Q, 1);  % number of stations
            K = size(self.results{it,1}.Tran.Avg.Q, 2);  % number of classes

            % Check convergence per class (aligned with JAR structure)
            for k = 1:K
                % Build QEntry and QExit matrices (M x E) for this class
                QEntry = zeros(M, E);
                QExit = zeros(M, E);
                for e = 1:E
                    for i = 1:M
                        Qik_curr = self.results{it,e}.Tran.Avg.Q{i,k};
                        if isstruct(Qik_curr) && isfield(Qik_curr, 'metric') && ~isempty(Qik_curr.metric)
                            QExit(i,e) = Qik_curr.metric(1);
                        end
                        Qik_prev = self.results{it-1,e}.Tran.Avg.Q{i,k};
                        if isstruct(Qik_prev) && isfield(Qik_prev, 'metric') && ~isempty(Qik_prev.metric)
                            QEntry(i,e) = Qik_prev.metric(1);
                        end
                    end
                end

                % Compute max relative absolute difference using maxpe
                % maxpe computes max(abs(1 - approx./exact)) = max(abs((approx-exact)./exact))
                % This matches JAR's Matrix.maxAbsDiff() implementation
                maxDiff = maxpe(QExit(:), QEntry(:));
                if isempty(maxDiff)
                    maxDiff = 0;  % Handle case where all QEntry values are zero
                end
                if isnan(maxDiff) || isinf(maxDiff)
                    return  % Non-convergence on invalid values
                end
                if maxDiff >= self.options.iter_tol
                    return  % Not converged
                end
            end
            bool = true;
        end


        function runAnalyzer(self)
            % RUNANALYZER()
            % Run the ensemble solver iteration
            line_debug('ENV solver starting: nstages=%d, method=%s', self.getNumberOfModels, self.options.method);

            % Show library attribution if verbose and not yet shown
            if self.options.verbose ~= VerboseLevel.SILENT && ~GlobalConstants.isLibraryAttributionShown()
                libs = SolverENV.getLibrariesUsed([], self.options);
                if ~isempty(libs)
                    line_printf('The solver will leverage %s.\n', strjoin(libs, ', '));
                    GlobalConstants.setLibraryAttributionShown(true);
                end
            end

            iterate(self);
        end

        function init(self)
            % INIT()
            % Initialize the environment solver with enhanced data structures
            % aligned with JAR SolverEnv implementation
            line_debug('ENV solver init: initializing environment data structures');
            options = self.options;
            if isfield(options,'seed')
                Solver.resetRandomGeneratorSeed(options.seed);
            end
            self.envObj.init();

            % Initialize ServerNum and SRates matrices (aligned with JAR)
            E = self.getNumberOfModels;
            M = self.sn{1}.nstations;
            K = self.sn{1}.nclasses;

            self.ServerNum = cell(K, 1);
            self.SRates = cell(K, 1);
            for k = 1:K
                self.ServerNum{k} = zeros(M, E);
                self.SRates{k} = zeros(M, E);
                for e = 1:E
                    for m = 1:M
                        self.ServerNum{k}(m, e) = self.sn{e}.nservers(m);
                        self.SRates{k}(m, e) = self.sn{e}.rates(m, k);
                    end
                end
            end

            % Build rate matrix E0
            self.E0 = zeros(E, E);
            for e = 1:E
                for h = 1:E
                    if ~isa(self.envObj.env{e,h}, 'Disabled')
                        self.E0(e, h) = self.envObj.env{e,h}.getRate();
                    end
                end
            end
            self.Eutil = ctmc_makeinfgen(self.E0);

            % Initialize transition CDFs
            self.transitionCdfs = cell(E, E);
            for e = 1:E
                for h = 1:E
                    if ~isa(self.envObj.env{e,h}, 'Disabled')
                        envDist = self.envObj.env{e,h};
                        self.transitionCdfs{e,h} = @(t) envDist.evalCDF(t);
                    else
                        self.transitionCdfs{e,h} = @(t) 0;
                    end
                end
            end

            % Initialize sojourn CDFs
            self.sojournCdfs = cell(E, 1);
            for e = 1:E
                self.sojournCdfs{e} = @(t) self.computeSojournCdf(e, t);
            end

            % newMethod: Use DTMC-based computation instead of CTMC
            % Verified numerical integration for Semi-Markov Process DTMC transition probabilities
            if self.newMethod
                line_debug('ENV using DTMC-based computation (newMethod=true)');
                self.dtmcP = zeros(E, E);
                for k = 1:E
                    for e = 1:E
                        if k == e || isa(self.envObj.env{k,e}, 'Disabled')
                            self.dtmcP(k, e) = 0.0;
                        else
                            % Compute the upper limit of the sojourn time
                            epsilon = 1e-8;
                            T = 1;
                            while self.transitionCdfs{k,e}(T) < 1.0 - epsilon
                                T = T * 2;
                                if T > 1e6  % safety limit
                                    break;
                                end
                            end
                            % Adaptive number of integration intervals based on T
                            N = max(1000, round(T * 100));
                            dt = T / N;
                            sumVal = 0;
                            for i = 0:(N-1)
                                t0 = i * dt;
                                t1 = t0 + dt;
                                deltaF = self.transitionCdfs{k,e}(t1) - self.transitionCdfs{k,e}(t0);
                                survival = 1;
                                for h = 1:E
                                    if h ~= k && h ~= e && ~isa(self.envObj.env{k,h}, 'Disabled')
                                        % Use midpoint for better accuracy
                                        tmid = (t0 + t1) / 2.0;
                                        survival = survival * (1.0 - self.envObj.env{k,h}.evalCDF(tmid));
                                    end
                                end
                                sumVal = sumVal + deltaF * survival;
                            end
                            self.dtmcP(k, e) = sumVal;
                        end
                    end
                end

                % Solve DTMC for stationary distribution
                dtmcPie = dtmc_solve(self.dtmcP);

                % Calculate hold times using numerical integration
                self.holdTimeMatrix = self.computeHoldTime(E);

                % Compute steady-state probabilities
                pi = zeros(1, E);
                denomSum = 0;
                for e = 1:E
                    denomSum = denomSum + dtmcPie(e) * self.holdTimeMatrix(e);
                end
                for k = 1:E
                    pi(k) = dtmcPie(k) * self.holdTimeMatrix(k) / denomSum;
                end
                self.envObj.probEnv = pi;

                % Update embedding weights
                newEmbweight = zeros(E, E);
                for e = 1:E
                    sumVal = 0.0;
                    for h = 1:E
                        if h ~= e
                            sumVal = sumVal + pi(h) * self.E0(h, e);
                        end
                    end
                    for k = 1:E
                        if k == e
                            newEmbweight(k, e) = 0;
                        else
                            if sumVal > 0
                                newEmbweight(k, e) = pi(k) * self.E0(k, e) / sumVal;
                            end
                        end
                    end
                end
                self.envObj.probOrig = newEmbweight;
            end

            % Compression: Use Courtois decomposition for large environments
            if self.compression
                line_debug('ENV using compression (Courtois decomposition)');
                self.applyCompression(E, M, K);
            end
        end

        function applyCompression(self, E, M, K)
            % APPLYCOMPRESSION Apply Courtois decomposition to reduce environment size
            % This method finds a good partition of the environment states and
            % creates compressed macro-state networks.

            % Find best partition
            if E <= 10
                self.MS = self.findBestPartition(E);
            else
                % Beam search for large environments
                self.MS = self.beamSearchPartition(E);
            end

            if isempty(self.MS)
                % No compression possible, use singletons
                self.MS = cell(E, 1);
                for i = 1:E
                    self.MS{i} = i;
                end
                self.Ecompress = E;
                return;
            end

            self.Ecompress = length(self.MS);

            % Apply decomposition/aggregation
            [p, eps, epsMax, q] = self.ctmc_decompose(self.Eutil, self.MS);

            if eps > epsMax
                line_warning(mfilename, 'Environment cannot be effectively compressed (eps > epsMax).');
            end

            % Store compression results
            self.compressionResult.p = p;
            self.compressionResult.eps = eps;
            self.compressionResult.epsMax = epsMax;
            self.compressionResult.q = q;

            % Update probEnv with macro-state probabilities
            pMacro = zeros(1, self.Ecompress);
            for i = 1:self.Ecompress
                pMacro(i) = sum(p(self.MS{i}));
            end
            self.envObj.probEnv = pMacro;

            % Compute micro-state probabilities within each macro-state
            pmicro = zeros(E, 1);
            for i = 1:self.Ecompress
                blockProb = p(self.MS{i});
                if sum(blockProb) > 0
                    pmicro(self.MS{i}) = blockProb / sum(blockProb);
                end
            end
            self.compressionResult.pmicro = pmicro;
            self.compressionResult.pMacro = pMacro;

            % Update embedding weights for macro-states
            Ecomp = self.Ecompress;
            newEmbweight = zeros(Ecomp, Ecomp);
            for e = 1:Ecomp
                sumVal = 0.0;
                for h = 1:Ecomp
                    if h ~= e
                        sumVal = sumVal + pMacro(h) * self.computeMacroRate(h, e);
                    end
                end
                for k = 1:Ecomp
                    if k == e
                        newEmbweight(k, e) = 0;
                    elseif sumVal > 0
                        newEmbweight(k, e) = pMacro(k) * self.computeMacroRate(k, e) / sumVal;
                    end
                end
            end
            self.envObj.probOrig = newEmbweight;

            % Build macro-state networks with weighted-average rates
            macroEnsemble = cell(self.Ecompress, 1);
            macroSolvers = cell(self.Ecompress, 1);
            macroSn = cell(self.Ecompress, 1);

            for i = 1:self.Ecompress
                % Copy the first micro-state network
                firstMicro = self.MS{i}(1);
                macroEnsemble{i} = self.ensemble{firstMicro}.copy();

                % Compute weighted-average rates
                for m = 1:M
                    for k = 1:K
                        rateSum = 0;
                        for r = 1:length(self.MS{i})
                            microIdx = self.MS{i}(r);
                            w = pmicro(microIdx);
                            rateSum = rateSum + w * self.sn{microIdx}.rates(m, k);
                        end

                        % Update service rate
                        jobclass = macroEnsemble{i}.classes{k};
                        station = macroEnsemble{i}.stations{m};
                        if isa(station, 'Queue') || isa(station, 'Delay')
                            if rateSum > 0
                                station.setService(jobclass, Exp(rateSum));
                            end
                        end
                    end
                end

                macroEnsemble{i}.refreshStruct(true);
                macroSn{i} = macroEnsemble{i}.getStruct(true);

                % Create solver for macro-state
                % Use the same solver factory pattern as original
                macroSolvers{i} = SolverFluid(macroEnsemble{i}, self.solvers{firstMicro}.options);
            end

            % Replace ensemble and solvers with compressed versions
            self.ensemble = macroEnsemble;
            self.solvers = macroSolvers;
            self.sn = macroSn;
        end

        function rate = computeMacroRate(self, fromMacro, toMacro)
            % COMPUTEMACRORATE Compute transition rate between macro-states
            rate = 0;
            for i = 1:length(self.MS{fromMacro})
                mi = self.MS{fromMacro}(i);
                for j = 1:length(self.MS{toMacro})
                    mj = self.MS{toMacro}(j);
                    rate = rate + self.compressionResult.pmicro(mi) * self.E0(mi, mj);
                end
            end
        end

        function MS = findBestPartition(self, E)
            % FINDBESTPARTITION Find the best partition for small environments (E <= 10)
            % Uses exhaustive search over all possible partitions

            % Start with singletons
            bestMS = cell(E, 1);
            for i = 1:E
                bestMS{i} = i;
            end

            [~, bestEps, bestEpsMax] = self.ctmc_decompose(self.Eutil, bestMS);
            if isempty(bestEps) || isnan(bestEps)
                MS = bestMS;
                return;
            end

            % Try merging pairs
            for i = 1:E
                for j = (i+1):E
                    testMS = cell(E-1, 1);
                    idx = 1;
                    for k = 1:E
                        if k == i
                            testMS{idx} = [i, j];
                            idx = idx + 1;
                        elseif k ~= j
                            testMS{idx} = k;
                            idx = idx + 1;
                        end
                    end

                    [~, testEps, testEpsMax] = self.ctmc_decompose(self.Eutil, testMS);
                    if ~isempty(testEps) && ~isnan(testEps) && testEps < bestEps
                        bestEps = testEps;
                        bestEpsMax = testEpsMax;
                        bestMS = testMS;
                    end
                end
            end

            MS = bestMS;
            self.Ecompress = length(MS);
        end

        function MS = beamSearchPartition(self, E)
            % BEAMSEARCHPARTITION Beam search for large environments (E > 10)

            B = 3;  % Beam width
            alpha = 0.01;  % Coupling threshold

            if isfield(self.options, 'config') && isfield(self.options.config, 'env_alpha')
                alpha = self.options.config.env_alpha;
            end

            % Initialize with singletons
            singletons = cell(E, 1);
            for i = 1:E
                singletons{i} = i;
            end
            beam = {singletons};

            bestSeen = singletons;
            [~, bestEps] = self.ctmc_decompose(self.Eutil, bestSeen);

            % Iteratively merge blocks
            for depth = 1:(E-1)
                candidates = {};

                for b = 1:length(beam)
                    ms = beam{b};
                    nBlocks = length(ms);

                    % Try all pairwise merges
                    for i = 1:nBlocks
                        for j = (i+1):nBlocks
                            % Create merged partition
                            trial = cell(nBlocks - 1, 1);
                            idx = 1;
                            for k = 1:nBlocks
                                if k == i
                                    trial{idx} = [ms{i}(:); ms{j}(:)];
                                    idx = idx + 1;
                                elseif k ~= j
                                    trial{idx} = ms{k};
                                    idx = idx + 1;
                                end
                            end

                            [~, childEps, childEpsMax] = self.ctmc_decompose(self.Eutil, trial);

                            if ~isempty(childEps) && ~isnan(childEps) && childEps > 0
                                cost = childEps - childEpsMax + alpha * depth;
                                candidates{end+1} = {trial, cost};

                                if cost < bestEps
                                    bestEps = cost;
                                    bestSeen = trial;
                                end
                            end
                        end
                    end
                end

                if isempty(candidates)
                    break;
                end

                % Sort by cost and keep top B
                costs = cellfun(@(x) x{2}, candidates);
                [~, sortIdx] = sort(costs);
                beam = {};
                for i = 1:min(B, length(sortIdx))
                    beam{end+1} = candidates{sortIdx(i)}{1};
                end
            end

            MS = bestSeen;
            self.Ecompress = length(MS);
        end

        function holdTime = computeHoldTime(self, E)
            % COMPUTEHOLDTIME Compute expected holding times using numerical integration
            holdTime = zeros(1, E);
            for k = 1:E
                % Survival function: 1 - sojournCDF
                surv = @(t) 1 - self.sojournCdfs{k}(t);

                % Compute upper limit
                upperLimit = 10;
                while surv(upperLimit) > 1e-8
                    upperLimit = upperLimit * 2;
                    if upperLimit > 1e6  % safety limit
                        break;
                    end
                end

                % Simpson's rule integration
                N = 10000;
                dt = upperLimit / N;
                integral = 0;
                for i = 0:(N-1)
                    t0 = i * dt;
                    t1 = t0 + dt;
                    tmid = (t0 + t1) / 2;
                    % Simpson's rule: (f(a) + 4*f(mid) + f(b)) * h/6
                    integral = integral + (surv(t0) + 4*surv(tmid) + surv(t1)) * dt / 6;
                end
                holdTime(k) = integral;
            end
        end

        function cdf = computeSojournCdf(self, e, t)
            % COMPUTESOJOURNCDF Compute sojourn time CDF for environment stage e
            E = self.getNumberOfModels;
            surv = 1.0;
            for h = 1:E
                if h ~= e
                    surv = surv * (1 - self.transitionCdfs{e,h}(t));
                end
            end
            cdf = 1 - surv;
        end

        function pre(self, it)
            % PRE(IT)

            if it==1
                for e=list(self)
                    if isinf(self.getSolver(e).options.timespan(2))
                        [QN,~,~,~] = self.getSolver(e).getAvg();
                    else
                        [QNt,~,~] = self.getSolver(e).getTranAvg();
                        % Handle case where getTranAvg returns NaN for disabled/empty results
                        QN = zeros(size(QNt));
                        for i = 1:size(QNt,1)
                            for k = 1:size(QNt,2)
                                if isstruct(QNt{i,k}) && isfield(QNt{i,k}, 'metric')
                                    QN(i,k) = QNt{i,k}.metric(end);
                                else
                                    QN(i,k) = 0; % Default to 0 for NaN/invalid entries
                                end
                            end
                        end
                    end
                    self.ensemble{e}.initFromMarginal(QN);
                end
            end
        end

        % solves model in stage e
        function [results_e, runtime] = analyze(self, it, e)
            % [RESULTS_E, RUNTIME] = ANALYZE(IT, E)
            results_e = struct();
            results_e.('Tran') = struct();
            results_e.Tran.('Avg') = [];
            T0 = tic;
            runtime = toc(T0);
            %% initialize
            [Qt,Ut,Tt] = self.ensemble{e}.getTranHandles;
            self.solvers{e}.reset();
            [QNt,UNt,TNt] = self.solvers{e}.getTranAvg(Qt,Ut,Tt);
            results_e.Tran.Avg.Q = QNt; 
            results_e.Tran.Avg.U = UNt; 
            results_e.Tran.Avg.T = TNt; 
        end


        function post(self, it)
            % POST(IT)

            E = self.getNumberOfModels;
            for e=1:E
                for h = 1:E
                    Qexit{e,h} = zeros(size(self.results{it,e}.Tran.Avg.Q));
                    Uexit{e,h} = zeros(size(self.results{it,e}.Tran.Avg.U));
                    Texit{e,h} = zeros(size(self.results{it,e}.Tran.Avg.T));
                    for i=1:size(self.results{it,e}.Tran.Avg.Q,1)
                        for r=1:size(self.results{it,e}.Tran.Avg.Q,2)
                            Qir = self.results{it,e}.Tran.Avg.Q{i,r};
                            % Check if result is a valid struct with required fields
                            if isstruct(Qir) && isfield(Qir, 't') && isfield(Qir, 'metric') && ~isempty(Qir.t)
                                w{e,h} = [0, map_cdf(self.envObj.proc{e}{h}, Qir.t(2:end)) - map_cdf(self.envObj.proc{e}{h}, Qir.t(1:end-1))]';
                                if ~isnan(w{e,h})
                                    Qexit{e,h}(i,r) = Qir.metric'*w{e,h}/sum(w{e,h});
                                    Uir = self.results{it,e}.Tran.Avg.U{i,r};
                                    if isstruct(Uir) && isfield(Uir, 'metric') && ~isempty(Uir.metric)
                                        Uexit{e,h}(i,r) = Uir.metric'*w{e,h}/sum(w{e,h});
                                    end
                                    Tir = self.results{it,e}.Tran.Avg.T{i,r};
                                    if isstruct(Tir) && isfield(Tir, 'metric') && ~isempty(Tir.metric)
                                        Texit{e,h}(i,r) = Tir.metric'*w{e,h}/sum(w{e,h});
                                    end
                                else
                                    Qexit{e,h}(i,r) = 0;
                                    Uexit{e,h}(i,r) = 0;
                                    Texit{e,h}(i,r) = 0;
                                end
                            else
                                w{e,h} = 0;
                            end
                        end
                    end
                end
            end

            Qentry = cell(1,E); % average entry queue-length
            for e = 1:E
                Qentry{e} = zeros(size(Qexit{e}));
                for h=1:E
                    % probability of coming from h to e \times resetFun(Qexit from h to e
                    if self.envObj.probOrig(h,e) > 0
                        Qentry{e} = Qentry{e} + self.envObj.probOrig(h,e) * self.resetFromMarginal{h,e}(Qexit{h,e});
                    end
                end
                self.solvers{e}.reset();
                self.ensemble{e}.initFromMarginal(Qentry{e});
            end

            % Update transition rates between stages if state-dependent method
            if isfield(self.options, 'method') && strcmp(self.options.method, 'statedep')
                for e = 1:E
                    for h = 1:E
                        if ~isa(self.envObj.env{e,h}, 'Disabled') && ~isempty(self.resetEnvRates{e,h})
                            self.envObj.env{e,h} = self.resetEnvRates{e,h}(...
                                self.envObj.env{e,h}, Qexit{e,h}, Uexit{e,h}, Texit{e,h});
                        end
                    end
                end
                % Reinitialize environment after rate updates
                self.envObj.init();
            end
        end

        function finish(self)
            % FINISH()

            it = size(self.results,1); % use last iteration
            E = self.getNumberOfModels;
            for e=1:E
                QExit{e}=[];
                UExit{e}=[];
                TExit{e}=[];
                if it>0
                    for i=1:size(self.results{it,e}.Tran.Avg.Q,1)
                        for r=1:size(self.results{it,e}.Tran.Avg.Q,2)
                            Qir = self.results{it,e}.Tran.Avg.Q{i,r};
                            % Check if result is a valid struct with required fields
                            if isstruct(Qir) && isfield(Qir, 't') && isfield(Qir, 'metric') && ~isempty(Qir.t)
                                w{e} = [0, map_cdf(self.envObj.holdTime{e}, Qir.t(2:end)) - map_cdf(self.envObj.holdTime{e}, Qir.t(1:end-1))]';
                                QExit{e}(i,r) = Qir.metric'*w{e}/sum(w{e});
                                Uir = self.results{it,e}.Tran.Avg.U{i,r};
                                if isstruct(Uir) && isfield(Uir, 'metric') && ~isempty(Uir.metric)
                                    UExit{e}(i,r) = Uir.metric'*w{e}/sum(w{e});
                                else
                                    UExit{e}(i,r) = 0;
                                end
                                Tir = self.results{it,e}.Tran.Avg.T{i,r};
                                if isstruct(Tir) && isfield(Tir, 'metric') && ~isempty(Tir.metric)
                                    TExit{e}(i,r) = Tir.metric'*w{e}/sum(w{e});
                                else
                                    TExit{e}(i,r) = 0;
                                end
                            else
                                QExit{e}(i,r) = 0;
                                UExit{e}(i,r) = 0;
                                TExit{e}(i,r) = 0;
                            end
                        end
                    end
                end
                %                 for h = 1:E
                %                     QE{e,h} = zeros(size(self.results{it,e}.Tran.Avg.Q));
                %                     for i=1:size(self.results{it,e}.Tran.Avg.Q,1)
                %                         for r=1:size(self.results{it,e}.Tran.Avg.Q,2)
                %                             w{e,h} = [0, map_cdf(self.envObj.proc{e}{h}, self.results{it,e}.Tran.Avg.Q{i,r}(2:end,2)) - map_cdf(self.envObj.proc{e}{h}, self.results{it,e}.Tran.Avg.Q{i,r}(1:end-1,2))]';
                %                             if ~isnan(w{e,h})
                %                                 QE{e,h}(i,r) = self.results{it,e}.Tran.Avg.Q{i,r}(:,1)'*w{e,h}/sum(w{e,h});
                %                             else
                %                                 QE{e,h}(i,r) = 0;
                %                             end
                %                         end
                %                     end
                %                 end
            end

            Qval=0*QExit{e};
            Uval=0*UExit{e};
            Tval=0*TExit{e};
            for e=1:E
                Qval = Qval + self.envObj.probEnv(e) * QExit{e}; % to check
                Uval = Uval + self.envObj.probEnv(e) * UExit{e}; % to check
                Tval = Tval + self.envObj.probEnv(e) * TExit{e}; % to check
            end
            self.result.Avg.Q = Qval;
            %    self.result.Avg.R = R;
            %    self.result.Avg.X = X;
            self.result.Avg.U = Uval;
            self.result.Avg.T = Tval;
            %    self.result.Avg.C = C;
            %self.result.runtime = runtime;
            %if self.options.verbose
            %    line_printf('\n');
            %end
        end

        function name = getName(self)
            % NAME = GETNAME()

            name = mfilename;
        end

        function [renvInfGen, stageInfGen, renvEventFilt, stageEventFilt, renvEvents, stageEvents] = getGenerator(self)
            % [renvInfGen, stageInfGen, renvEventFilt, stageEventFilt, renvEvents, stageEvents] = getGenerator(self)
            %
            % Returns the infinitesimal generator matrices for the random environment model.
            %
            % Outputs:
            %   renvInfGen     - Combined infinitesimal generator for the random environment (flattened)
            %   stageInfGen    - Cell array of infinitesimal generators for each stage
            %   renvEventFilt  - Cell array (E x E) of event filtration matrices for environment transitions
            %   stageEventFilt - Cell array of event filtration matrices for each stage
            %   renvEvents     - Cell array of Event objects for environment transitions
            %   stageEvents    - Cell array of synchronization maps for each stage

            E = self.getNumberOfModels;
            stageInfGen = cell(1,E);
            stageEventFilt = cell(1,E);
            stageEvents = cell(1,E);
            for e=1:E
                if isa(self.solvers{e},'SolverCTMC')
                    [stageInfGen{e}, stageEventFilt{e}, stageEvents{e}] = self.solvers{e}.getGenerator();
                else
                    line_error(mfilename,'This method requires SolverENV to be instantiated with the CTMC solver.');
                end
            end

            % Get number of states for each stage
            nstates = cellfun(@(g) size(g, 1), stageInfGen);

            % Get number of phases for each transition distribution
            nphases = zeros(E, E);
            for i=1:E
                for j=1:E
                    if ~isempty(self.env{i,j}) && ~isa(self.env{i,j}, 'Disabled')
                        nphases(i,j) = self.env{i,j}.getNumberOfPhases();
                    else
                        nphases(i,j) = 1;
                    end
                end
            end
            % Adjust diagonal (self-transitions have one less phase in the Kronecker expansion)
            nphases = nphases - eye(E);

            % Initialize block cell structure for the random environment generator
            renvInfGen = cell(E,E);

            for e=1:E
                % Diagonal block: stage infinitesimal generator
                renvInfGen{e,e} = stageInfGen{e};
                for h=1:E
                    if h~=e
                        % Off-diagonal blocks: reset matrices (identity with appropriate dimensions)
                        minStates = min(nstates(e), nstates(h));
                        resetMatrix_eh = sparse(nstates(e), nstates(h));
                        for i=1:minStates
                            resetMatrix_eh(i,i) = 1;
                        end
                        renvInfGen{e,h} = resetMatrix_eh;
                    end
                end
            end

            % Build environment transition events and expand generator with phase structure
            renvEvents = cell(1,0);
            for e=1:E
                for h=1:E
                    if h~=e
                        % Get D0 (phase generator) for transition from e to h
                        if isempty(self.env{e,h}) || isa(self.env{e,h}, 'Disabled')
                            D0 = zeros(0);
                        else
                            proc = self.env{e,h}.getProcess();
                            D0 = proc{1};
                        end
                        % Kronecker sum with diagonal block
                        renvInfGen{e,e} = krons(renvInfGen{e,e}, D0);

                        % Get D1 (completion rate matrix) and initial probability vector pie
                        if isempty(self.env{h,e}) || isa(self.env{h,e}, 'Disabled') || any(isnan(map_pie(self.env{h,e}.getProcess())))
                            pie = ones(1, nphases(h,e));
                        else
                            pie = map_pie(self.env{h,e}.getProcess());
                        end

                        if isempty(self.env{e,h}) || isa(self.env{e,h}, 'Disabled')
                            D1 = zeros(0);
                        else
                            proc = self.env{e,h}.getProcess();
                            D1 = proc{2};
                        end

                        % Kronecker product for off-diagonal block
                        onePhase = ones(nphases(e,h), 1);
                        kronArg = D1 * onePhase * pie;
                        renvInfGen{e,h} = kron(renvInfGen{e,h}, sparse(kronArg));

                        % Create environment transition events
                        for i=1:self.ensemble{e}.getNumberOfNodes
                            renvEvents{1,end+1} = Event(EventType.STAGE, i, NaN, NaN, [e,h]);  %#ok<AGROW>
                        end

                        % Handle other stages (f != e, f != h)
                        for f=1:E
                            if f~=h && f~=e
                                if isempty(self.env{f,h}) || isa(self.env{f,h}, 'Disabled') || any(isnan(map_pie(self.env{f,h}.getProcess())))
                                    pie_fh = ones(1, nphases(f,h));
                                else
                                    pie_fh = map_pie(self.env{f,h}.getProcess());
                                end
                                oneVec = ones(nphases(e,h), 1);
                                renvInfGen{e,f} = kron(renvInfGen{e,f}, oneVec * pie_fh);
                            end
                        end
                    end
                end
            end

            % Build event filtration matrices for environment transitions
            % Each renvEventFilt{e,h} isolates transitions from stage e to stage h
            renvEventFilt = cell(E,E);
            for e=1:E
                for h=1:E
                    tmpCell = cell(E,E);
                    for e1=1:E
                        for h1=1:E
                            tmpCell{e1,h1} = renvInfGen{e1,h1};
                        end
                    end
                    % Zero out diagonal blocks (internal stage transitions)
                    for e1=1:E
                        tmpCell{e1,e1} = tmpCell{e1,e1} * 0;
                    end
                    % Zero out off-diagonal blocks that don't match (e,h) transition
                    for e1=1:E
                        for h1=1:E
                            if e~=e1 || h~=h1
                                if e1~=h1  % Only zero out off-diagonal entries
                                    tmpCell{e1,h1} = tmpCell{e1,h1} * 0;
                                end
                            end
                        end
                    end
                    renvEventFilt{e,h} = cell2mat(tmpCell);
                end
            end

            % Flatten block structure into single matrix and normalize
            renvInfGen = cell2mat(renvInfGen);
            renvInfGen = ctmc_makeinfgen(renvInfGen);
        end

        function varargout = getAvg(varargin)
            % [QNCLASS, UNCLASS, TNCLASS] = GETAVG()
            [varargout{1:nargout}] = getEnsembleAvg( varargin{:} );
        end

        function [QNclass, UNclass, RNclass, TNclass, ANclass, WNclass] = getEnsembleAvg(self)
            % [QNCLASS, UNCLASS, TNCLASS] = GETENSEMBLEAVG()
            RNclass=[];
            ANclass=[];
            WNclass=[];

            if isempty(self.result) || (isfield(self.options,'force') && self.options.force)
                iterate(self);
                if isempty(self.result)
                    QNclass=[];
                    UNclass=[];
                    TNclass=[];
                    return
                end
            end
            QNclass = self.result.Avg.Q;
            UNclass = self.result.Avg.U;
            TNclass = self.result.Avg.T;
            WNclass = QNclass ./ TNclass;
            RNclass = NaN*WNclass;
            ANclass = NaN*TNclass;
        end

        function [AvgTable,QT,UT,TT] = getAvgTable(self,keepDisabled)
            % [AVGTABLE,QT,UT,TT] = GETAVGTABLE(SELF,KEEPDISABLED)
            % Return table of average station metrics

            if nargin<2 %if ~exist('keepDisabled','var')
                keepDisabled = false;
            end

            [QN,UN,~,TN] = getAvg(self);
            M = size(QN,1);
            K = size(QN,2);
            Q = self.result.Avg.Q;
            U = self.result.Avg.U;
            T = self.result.Avg.T;
            if isempty(QN)
                AvgTable = Table();
                QT = Table();
                UT = Table();
                TT = Table();
            elseif ~keepDisabled
                Qval = []; Uval = []; Tval = [];
                JobClass = {};
                Station = {};
                for ist=1:M
                    for k=1:K
                        if any(sum([QN(ist,k),UN(ist,k),TN(ist,k)])>0)
                            JobClass{end+1,1} = self.sn{1}.classnames{k};
                            Station{end+1,1} = self.sn{1}.nodenames{self.sn{1}.stationToNode(ist)};
                            Qval(end+1) = QN(ist,k);
                            Uval(end+1) = UN(ist,k);
                            Tval(end+1) = TN(ist,k);
                        end
                    end
                end
                QLen = Qval(:); % we need to save first in a variable named like the column
                QT = Table(Station,JobClass,QLen);
                Util = Uval(:); % we need to save first in a variable named like the column
                UT = Table(Station,JobClass,Util);
                Tput = Tval(:); % we need to save first in a variable named like the column
                Station = categorical(Station);
                JobClass = categorical(JobClass);
                TT = Table(Station,JobClass,Tput);
                RespT = QLen ./ Tput;
                AvgTable = Table(Station,JobClass,QLen,Util,RespT,Tput);
            else
                Qval = zeros(M,K); Uval = zeros(M,K);
                JobClass = cell(K*M,1);
                Station = cell(K*M,1);
                for ist=1:M
                    for k=1:K
                        JobClass{(ist-1)*K+k} = Q{ist,k}.class.name;
                        Station{(ist-1)*K+k} = Q{ist,k}.station.name;
                        Qval((ist-1)*K+k) = QN(ist,k);
                        Uval((ist-1)*K+k) = UN(ist,k);
                        Tval((ist-1)*K+k) = TN(ist,k);
                    end
                end
                Station = categorical(Station);
                JobClass = categorical(JobClass);
                QLen = Qval(:); % we need to save first in a variable named like the column
                QT = Table(Station,JobClass,QLen);
                Util = Uval(:); % we need to save first in a variable named like the column
                UT = Table(Station,JobClass,Util);
                Tput = Tval(:); % we need to save first in a variable named like the column
                TT = Table(Station,JobClass,Tput);
                RespT = QLen ./ Tput;
                AvgTable = Table(Station,JobClass,QLen,Util,RespT,Tput);
            end
        end

        function envsn = getStruct(self)
            E = self.getNumberOfModels;
            envsn = cell(1,E);
            for e=1:E
                envsn{e} = self.ensemble{e}.getStruct;
            end
        end

        function [T, segmentResults] = getSamplePathTable(self, samplePath)
            % [T, SEGMENTRESULTS] = GETSAMPLEPATHTABLE(SELF, SAMPLEPATH)
            %
            % Compute transient performance metrics for a sample path through
            % environment states. The method runs transient analysis for each
            % segment and extracts initial and final metric values.
            %
            % Inputs:
            %   samplePath - Cell array where each row is {stage, duration}
            %                stage: string (name) or integer (1-based index)
            %                duration: positive scalar (time spent in stage)
            %
            % Outputs:
            %   T - Table with columns: Segment, Stage, Duration, Station, JobClass,
            %       InitQLen, InitUtil, InitTput, FinalQLen, FinalUtil, FinalTput
            %   segmentResults - Cell array with detailed transient results per segment
            %
            % Example:
            %   samplePath = {'Fast', 5.0; 'Slow', 10.0; 'Fast', 3.0};
            %   [T, results] = solver.getSamplePathTable(samplePath);

            % Validate input
            if isempty(samplePath)
                line_error(mfilename, 'Sample path cannot be empty.');
            end

            % Initialize if needed
            if isempty(self.envObj.probEnv)
                self.init();
            end

            E = self.getNumberOfModels;
            M = self.sn{1}.nstations;
            K = self.sn{1}.nclasses;
            nSegments = size(samplePath, 1);
            segmentResults = cell(nSegments, 1);

            % Initialize queue lengths (uniform distribution for closed classes)
            Q_current = zeros(M, K);
            for k = 1:K
                if self.sn{1}.njobs(k) > 0  % closed class
                    Q_current(:, k) = self.sn{1}.njobs(k) / M;
                end
            end

            % Process each segment
            for seg = 1:nSegments
                % Resolve stage
                stageSpec = samplePath{seg, 1};
                duration = samplePath{seg, 2};

                if ischar(stageSpec) || isstring(stageSpec)
                    e = self.envObj.envGraph.findnode(stageSpec);
                    if e == 0
                        line_error(mfilename, sprintf('Stage "%s" not found.', stageSpec));
                    end
                    stageName = char(stageSpec);
                else
                    e = stageSpec;
                    if e < 1 || e > E
                        line_error(mfilename, sprintf('Stage index %d out of range [1, %d].', e, E));
                    end
                    stageName = self.envObj.envGraph.Nodes.Name{e};
                end

                if duration <= 0
                    line_error(mfilename, 'Duration must be positive.');
                end

                % Initialize model from current queue lengths
                self.ensemble{e}.initFromMarginal(Q_current);

                % Set solver timespan
                self.solvers{e}.options.timespan = [0, duration];
                self.solvers{e}.reset();

                % Run transient analysis
                [Qt, Ut, Tt] = self.ensemble{e}.getTranHandles;
                [QNt, UNt, TNt] = self.solvers{e}.getTranAvg(Qt, Ut, Tt);

                % Extract initial and final metrics
                initQ = zeros(M, K);
                initU = zeros(M, K);
                initT = zeros(M, K);
                finalQ = zeros(M, K);
                finalU = zeros(M, K);
                finalT = zeros(M, K);

                for i = 1:M
                    for r = 1:K
                        Qir = QNt{i, r};
                        Uir = UNt{i, r};
                        Tir = TNt{i, r};

                        if isstruct(Qir) && isfield(Qir, 'metric') && ~isempty(Qir.metric)
                            initQ(i, r) = Qir.metric(1);
                            finalQ(i, r) = Qir.metric(end);
                        end
                        if isstruct(Uir) && isfield(Uir, 'metric') && ~isempty(Uir.metric)
                            initU(i, r) = Uir.metric(1);
                            finalU(i, r) = Uir.metric(end);
                        end
                        if isstruct(Tir) && isfield(Tir, 'metric') && ~isempty(Tir.metric)
                            initT(i, r) = Tir.metric(1);
                            finalT(i, r) = Tir.metric(end);
                        end
                    end
                end

                % Store segment results
                segmentResults{seg}.stage = e;
                segmentResults{seg}.stageName = stageName;
                segmentResults{seg}.duration = duration;
                segmentResults{seg}.QNt = QNt;
                segmentResults{seg}.UNt = UNt;
                segmentResults{seg}.TNt = TNt;
                segmentResults{seg}.initQ = initQ;
                segmentResults{seg}.initU = initU;
                segmentResults{seg}.initT = initT;
                segmentResults{seg}.finalQ = finalQ;
                segmentResults{seg}.finalU = finalU;
                segmentResults{seg}.finalT = finalT;

                % Update current queue lengths for next segment
                Q_current = finalQ;
            end

            % Build output table
            Segment = [];
            Stage = {};
            Duration = [];
            Station = {};
            JobClass = {};
            InitQLen = [];
            InitUtil = [];
            InitTput = [];
            FinalQLen = [];
            FinalUtil = [];
            FinalTput = [];

            for seg = 1:nSegments
                res = segmentResults{seg};
                for i = 1:M
                    for r = 1:K
                        % Only include rows with non-zero metrics
                        if any([res.initQ(i,r), res.initU(i,r), res.initT(i,r), ...
                                res.finalQ(i,r), res.finalU(i,r), res.finalT(i,r)] > 0)
                            Segment(end+1, 1) = seg;
                            Stage{end+1, 1} = res.stageName;
                            Duration(end+1, 1) = res.duration;
                            Station{end+1, 1} = self.sn{1}.nodenames{self.sn{1}.stationToNode(i)};
                            JobClass{end+1, 1} = self.sn{1}.classnames{r};
                            InitQLen(end+1, 1) = res.initQ(i, r);
                            InitUtil(end+1, 1) = res.initU(i, r);
                            InitTput(end+1, 1) = res.initT(i, r);
                            FinalQLen(end+1, 1) = res.finalQ(i, r);
                            FinalUtil(end+1, 1) = res.finalU(i, r);
                            FinalTput(end+1, 1) = res.finalT(i, r);
                        end
                    end
                end
            end

            Stage = categorical(Stage);
            Station = categorical(Station);
            JobClass = categorical(JobClass);
            T = Table(Segment, Stage, Duration, Station, JobClass, ...
                      InitQLen, InitUtil, InitTput, FinalQLen, FinalUtil, FinalTput);
        end
        function [allMethods] = listValidMethods(self)
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            sn = self.model.getStruct();
            allMethods = {'default'};
        end
    end

    methods (Static)
        
        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)

            featUsed = model.getUsedLangFeatures();

            featSupported = SolverFeatureSet;

            % Nodes
            featSupported.setTrue('ClassSwitch');
            featSupported.setTrue('Delay');
            featSupported.setTrue('DelayStation');
            featSupported.setTrue('Queue');
            featSupported.setTrue('Sink');
            featSupported.setTrue('Source');

            % Distributions
            featSupported.setTrue('Coxian');
            featSupported.setTrue('Cox2');
            featSupported.setTrue('Erlang');
            featSupported.setTrue('Exp');
            featSupported.setTrue('HyperExp');

            % Sections
            featSupported.setTrue('StatelessClassSwitcher'); % Section
            featSupported.setTrue('InfiniteServer'); % Section
            featSupported.setTrue('SharedServer'); % Section
            featSupported.setTrue('Buffer'); % Section
            featSupported.setTrue('Dispatcher'); % Section
            featSupported.setTrue('Server'); % Section (Non-preemptive)
            featSupported.setTrue('JobSink'); % Section
            featSupported.setTrue('RandomSource'); % Section
            featSupported.setTrue('ServiceTunnel'); % Section

            % Scheduling strategy
            featSupported.setTrue('SchedStrategy_INF');
            featSupported.setTrue('SchedStrategy_PS');
            featSupported.setTrue('SchedStrategy_FCFS');
            featSupported.setTrue('RoutingStrategy_PROB');
            featSupported.setTrue('RoutingStrategy_RAND');
            featSupported.setTrue('RoutingStrategy_RROBIN'); % with SolverJMT

            % Customer Classes
            featSupported.setTrue('ClosedClass');
            featSupported.setTrue('OpenClass');

            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end
    end

    methods (Static)
        % ensemble solver options
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = SolverOptions('ENV');
        end

        function libs = getLibrariesUsed(sn, options)
            % GETLIBRARIESUSED Get list of external libraries used by ENV solver
            % ENV uses internal algorithms, no external library attribution needed
            libs = {};
        end
    end
end

