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
        horizonWarned = []; % stages whose unspecified horizon has been announced
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
        % Analyzer strategy: which inter-stage coupling to use.
        analyzerMode = 'meanfield';                  % 'meanfield' (default) | 'statevec'
        analyzerFcn = @solver_env_meanfield_analyzer; % phase-dispatched analyzer handle
        % State-vector analyzer working data (options.method='statevec' only)
        Qgen;        % Cell{E}: per-stage infinitesimal generator
        SS;          % Cell{E}: per-stage state space (rows = states)
        SSaggr;      % Cell{E}: per-stage aggregated (marginal) state space
        piEnter;     % Cell{E}: per-stage entry distribution (row vector)
        piEnterPrev; % Cell{E}: previous-iteration entry distribution (convergence)
        piExitDest;  % Cell{E}{E}: per-stage exit distribution toward each destination
        piTimeAvg;   % Cell{E}: per-stage sojourn-end distribution used in the blend
        statevecData;% Cell{E}: struct(arvRates,depRates,sn,options) per stage
        resetStateFun; % Cell{E,E}: state-vector reset maps (default identity)
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

            % State-vector reset maps (statevec analyzer). Default to identity
            % for any (e,h) the environment left unset.
            for e=1:length(self.env)
                for h=1:length(self.env)
                    if ~isempty(renv.resetStateFun) && e<=size(renv.resetStateFun,1) ...
                            && h<=size(renv.resetStateFun,2) && ~isempty(renv.resetStateFun{e,h})
                        self.resetStateFun{e,h} = renv.resetStateFun{e,h};
                    else
                        self.resetStateFun{e,h} = @(pi) pi;
                    end
                end
            end

            % The Markovian-arc rule used to be raised HERE, which made a
            % semi-Markov environment fail to construct at all, so model.help
            % (which builds a probe) could not even list the env family. It is
            % now a rule of methodRefusal, asked by the gate and by init.
            for e=1:length(self.env)
                for h=1:length(self.env)
                    if isa(self.env{e,h},'Disabled')
                        self.env{e,h} = Exp(0);
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
                    % numel, not size(MS,1): a row cell array of macro-states
                    % has size(MS,1)==1 and would build one singleton here.
                    nMacro = numel(MS);
                    MSS = cell(nMacro, 1);
                    for i = 1:nMacro
                        MSS{i} = i;
                    end
                    [p, ~, ~, eps, epsMax] = ctmc_multi(Q, MS, MSS);
                    q = 1.05 * max(max(abs(Q)));
                otherwise
                    line_error(mfilename, sprintf('Unknown decomposition method: %s', method));
            end
        end

        function bool = converged(self, it)
            % BOOL = CONVERGED(IT)  Convergence test, delegated to the active analyzer.
            bool = self.analyzerFcn(self, 'converged', it, []);
        end


        function runAnalyzer(self)
            % RUNANALYZER()
            % Run the ensemble solver iteration
            LineConsole.step('random environment with %d stages, method ''%s''', ...
                self.getNumberOfModels, self.options.method);
            line_debug('ENV solver starting: nstages=%d, method=%s', self.getNumberOfModels, self.options.method);

            % Show library attribution if verbose and not yet shown
            if self.options.verbose ~= VerboseLevel.SILENT && ~GlobalConstants.isLibraryAttributionShown()
                libs = SolverENV.getLibrariesUsed([], self.options);
                if ~isempty(libs)
                    line_printf('The solver will leverage %s.\n', strjoin(libs, ', '));
                    GlobalConstants.setLibraryAttributionShown(true);
                end
            end

            % Closed-form fast/slow environment limits bypass the transient
            % mean-field/statevec iteration (see solveEnvLimit).
            if isfield(self.options,'method') && any(strcmpi(self.options.method,{'avg','dec'}))
                self.solveEnvLimit();
                return
            end

            iterate(self);
        end

        function init(self)
            % INIT()
            % Initialize the environment solver with enhanced data structures
            % aligned with JAR SolverEnv implementation
            line_debug('ENV solver init: initializing environment data structures');
            options = self.options;
            % The gate's own sentence before Environment.init can index a
            % non-Markovian arc as a MAP or a stage can come back all zero.
            method = 'default';
            if isfield(options,'method') && ~isempty(options.method)
                method = options.method;
            end
            [okm, whym] = self.methodRefusal(method);
            if ~okm
                line_error(mfilename, whym);
            end
            if isfield(options,'seed')
                Solver.resetRandomGeneratorSeed(options.seed);
            end
            self.envObj.init();

            % Initialize ServerNum and SRates (flat station x class tables used
            % only by the statevec analyzer; skipped for LQN stages).
            % see _kb/06-solver-catalog.md for rationale
            E = self.getNumberOfModels;
            if isfield(self.sn{1}, 'nstations')
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

            % Select the inter-stage coupling analyzer (statevec vs default
            % mean-field). see _kb/06-solver-catalog.md for rationale
            if isfield(self.options,'method') && any(strcmpi(self.options.method,{'statevec','blend'}))
                % statevec needs a single per-stage CTMC generator, which an LQN
                % stage cannot provide. see _kb/06-solver-catalog.md for rationale
                for e=1:E
                    if isa(self.ensemble{e}, 'LayeredNetwork')
                        line_error(mfilename, sprintf(['The state-vector (statevec) analyzer does not support ' ...
                            'LayeredNetwork stages (stage %d): an LQN has no single stage generator. ' ...
                            'Use the default mean-field analyzer (omit method=''statevec'').'], e));
                    end
                end
                self.analyzerMode = 'statevec';
                self.analyzerFcn = @solver_env_statevec_analyzer;
            else
                self.analyzerMode = 'meanfield';
                self.analyzerFcn = @solver_env_meanfield_analyzer;
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
            % PRE(IT)  Delegated to the active analyzer (marginal | statevec).
            self.analyzerFcn(self, 'pre', it, []);
        end

        % solves model in stage e
        function [results_e, runtime] = analyze(self, it, e)
            % [RESULTS_E, RUNTIME] = ANALYZE(IT, E)  Delegated to the active analyzer.
            [results_e, runtime] = self.analyzerFcn(self, 'analyze', it, e);
        end


        function post(self, it)
            % POST(IT)  Delegated to the active analyzer.
            self.analyzerFcn(self, 'post', it, []);
        end

        function finish(self)
            % FINISH()  Delegated to the active analyzer.
            self.analyzerFcn(self, 'finish', [], []);
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
            % [QNCLASS, UNCLASS, RNCLASS, TNCLASS, ANCLASS, WNCLASS] = GETAVG()
            % Throughput is the FOURTH output, not the third: RNclass comes
            % third and is always NaN because ENV computes no response time.
            [varargout{1:nargout}] = getEnsembleAvg( varargin{:} );
        end

        function [QNclass, UNclass, RNclass, TNclass, ANclass, WNclass] = getEnsembleAvg(self)
            % [QNCLASS, UNCLASS, RNCLASS, TNCLASS, ANCLASS, WNCLASS] = GETENSEMBLEAVG()
            % Solver console: SolverENV drives an ensemble of stage models and
            % does not pass through runAnalyzerChecks, so it opens its own run
            % here, at the entry point every caller goes through. The guard
            % must live until this function returns.
            consoleGuard = LineConsole.beginRun(self, self.options); %#ok<NASGU>
            RNclass=[];
            ANclass=[];
            WNclass=[];
            % ANNOUNCE THE SOLVER, as every NetworkSolver does at
            % NetworkSolver.setAvgResults: this ensemble is not one, so its
            % result table used to arrive with NO banner and was read as a table
            % belonging to nobody -- parity-static tags it UNKNOWN, keeps a
            % golden under that name for want of anything better, and its
            % agreement-gated generator drops it as a parsing artifact, so
            % renv_threestages_repairmen could never earn a gated golden. The
            % text is the one the JAR's SolverENV prints, so both spell the
            % solver the same way. printBanner_ is called at each exit rather
            % than through onCleanup: the cleanup object is destroyed after the
            % method's own output assignment and its call was not reaching the
            % console at all, so an explicit call is the one that prints.
            envTstart = tic;

            if isfield(self.options,'lang') && strcmp(self.options.lang,'python')
                % the stage solvers this ensemble was built with decide the
                % bridge's stage-solver class AND its options; see
                % PYLINE.getEnvAvg
                % The forwarded horizon is the RESOLVED one, so the engine
                % integrates the interval the native path does instead of
                % dropping a non-finite value and falling back to its own
                % default of 100. On renv_fourstages_repairmen that was 30
                % against 100 -- a different question asked of each codebase,
                % and most of what the example's parity slack used to measure.
                % stageOptions is a struct COPY, so this states the horizon
                % without leaving it on the stage solver.
                self.assertOneStageHorizon_('python');
                stageSolverName = '';
                stageOptions = struct();
                if ~isempty(self.solvers) && ~isempty(self.solvers{1})
                    stageSolverName = class(self.solvers{1});
                    if isprop(self.solvers{1}, 'options') && isstruct(self.solvers{1}.options)
                        stageOptions = self.solvers{1}.options;
                        stageTs = self.stageHorizon_(1);
                        if numel(stageTs) >= 2
                            stageOptions.timespan = stageTs;
                        end
                    end
                end
                [QNclass, UNclass, TNclass] = PYLINE.getEnvAvg(self.envObj, self.options, stageSolverName, stageOptions);
                WNclass = QNclass ./ TNclass;
                RNclass = NaN*WNclass;
                ANclass = NaN*TNclass;
                self.result.Avg.Q = QNclass;
                self.result.Avg.U = UNclass;
                self.result.Avg.T = TNclass;
                self.printBanner_(toc(envTstart));
                return
            end

            if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
                % THE STAGE SOLVER IS NOT IN THE MODEL, which is the root of
                % this: model.json carries the stage Networks and the transition
                % process, and the ensemble's solver choice lives on this object.
                % An engine that defaults to the fluid transient therefore
                % answered a DIFFERENT model in silence -- on
                % renv_threestages_repairmen, whose stages are SolverCTMC, that
                % returned Queue1 throughput 1.5577 against the 1.3333 the CTMC
                % stages give. line-cli now takes --stage-solver, and its
                % mean-field coupling runs the enumerated CTMC as well as the
                % fluid transient, so CPPLINE.getEnvAvg names the ensemble's own.
                % A MIXED ensemble is still refused: one coupling runs one stage
                % solver, so sending the first stage's name would answer for the
                % rest under it.
                self.assertOneStageHorizon_('cpp');
                stageClass = '';
                if ~isempty(self.solvers) && ~isempty(self.solvers{1})
                    stageClass = class(self.solvers{1});
                end
                % isa, not strcmp: the example-facing aliases subclass the solver
                % (classdef FLD < SolverFluid), so a name test refuses FLD too.
                stageTok = '';
                allSame = true;
                for si = 1:numel(self.solvers)
                    if isempty(self.solvers{si}), continue; end
                    if isa(self.solvers{si},'SolverFluid')
                        tok = 'fluid';
                    elseif isa(self.solvers{si},'SolverCTMC')
                        tok = 'ctmc';
                    else
                        tok = '';
                    end
                    if isempty(tok) || (~isempty(stageTok) && ~strcmp(tok, stageTok))
                        allSame = false;
                        break
                    end
                    stageTok = tok;
                end
                if ~allSame || isempty(stageTok)
                    line_error(mfilename, sprintf(['SolverENV does not support ' ...
                        'lang=''cpp'' with %s stages: the C++ ENV engine solves every ' ...
                        'stage by the fluid transient or by the enumerated CTMC, and one ' ...
                        'coupling runs one stage solver. Use lang=''matlab'' or ' ...
                        'lang=''python'' for an ensemble built on another stage solver.'], ...
                        stageClass));
                end
                % The horizon the STAGE transients are integrated over travels
                % with the request, RESOLVED as in the python branch above; see
                % CPPLINE.getEnvAvg.
                cppOptions = self.options;
                stageTs = self.stageHorizon_(1);
                if numel(stageTs) >= 2
                    cppOptions.stagetimespan = stageTs;
                end
                % WHICH SOLVER RUNS EACH STAGE, which model.json does not carry.
                cppOptions.stagesolver = stageTok;
                if strcmp(stageTok,'ctmc') && ~isempty(self.solvers{1}) && ...
                        isfield(self.solvers{1}.options,'cutoff') && ...
                        isfinite(self.solvers{1}.options.cutoff)
                    cppOptions.stagecutoff = self.solvers{1}.options.cutoff;
                end
                [QNclass, UNclass, TNclass] = CPPLINE.getEnvAvg(self.envObj, self.ensemble{1}, cppOptions);
                WNclass = QNclass ./ TNclass;
                RNclass = NaN*WNclass;
                ANclass = NaN*TNclass;
                self.result.Avg.Q = QNclass;
                self.result.Avg.U = UNclass;
                self.result.Avg.T = TNclass;
                self.printBanner_(toc(envTstart));
                return
            end

            % A second getEnsembleAvg (getAvgTable calls it again) returns the
            % cached result without solving, and must not re-announce: a
            % NetworkSolver banners from setAvgResults, which a cached call
            % never reaches. The bridges below are exempt because they re-solve
            % on every call.
            didSolve = isempty(self.result) || (isfield(self.options,'force') && self.options.force);
            if didSolve
                if isfield(self.options,'method') && any(strcmpi(self.options.method,{'avg','dec'}))
                    self.solveEnvLimit();
                else
                    iterate(self);
                end
                if isempty(self.result)
                    QNclass=[];
                    UNclass=[];
                    TNclass=[];
                    self.printBanner_(toc(envTstart));
                    return
                end
            end
            QNclass = self.result.Avg.Q;
            UNclass = self.result.Avg.U;
            TNclass = self.result.Avg.T;
            WNclass = QNclass ./ TNclass;
            RNclass = NaN*WNclass;
            ANclass = NaN*TNclass;
            if didSolve
                self.printBanner_(toc(envTstart));
            end
        end

        function ts = stageHorizon_(self, e)
            % TS = STAGEHORIZON_(E)  Stage e's transient horizon, resolved.
            %
            % A stage solver left at timespan(2)=Inf -- what an example writes
            % when it has no particular horizon in mind -- was resolved inside
            % NetworkSolver.getTranAvg, into a LOCAL copy of the options: the
            % solver's own timespan stayed Inf, so the value was re-derived and
            % re-announced on every iteration, the bridges never saw it, and
            % cdfGrid_ (which returns empty for a non-finite horizon) left that
            % stage out of the sojourn quadrature grid. The rule below is
            % getTranAvg's, verbatim, so the horizon is the one that path would
            % have picked.
            %
            % PURE ON PURPOSE. Writing it onto the stage solver would outlive
            % the ENV solve, and a stage solver's own getters branch on whether
            % a horizon is set: renv_fourstages_repairmen ends with
            % getEnsembleAvgTables(), whose per-stage getAvgTable then took the
            % TRANSIENT arm, which lang='python' does not bridge -- turning that
            % row from a result into a refusal. Callers set it for the call that
            % needs it and put it back (see setStageTimespan_).
            ts = [];
            if isempty(self.solvers{e}) || ~isprop(self.solvers{e},'options') ...
                    || ~isfield(self.solvers{e}.options,'timespan')
                return
            end
            ts = self.solvers{e}.options.timespan;
            if numel(ts) < 2 || ~isinf(ts(2))
                return
            end
            rates = self.sn{e}.rates;
            minrate = min(rates(isfinite(rates)));
            ts(2) = 30/minrate;
            if isinf(ts(1))
                ts(1) = 0;
            end
            if numel(self.horizonWarned) < e || ~self.horizonWarned(e)
                self.horizonWarned(e) = true;
                line_warning(mfilename, ...
                    ['End time of transient analysis unspecified for stage %d, ' ...
                    'setting its timespan option to [%g,%g]. Pass a stage solver ' ...
                    'with ''timespan'',[0,T] to customize.\n'], e, ts(1), ts(2));
            end
        end

        function setStageTimespan_(self, e, ts)
            % SETSTAGETIMESPAN_(E, TS)  Put stage e's horizon back.
            if ~isempty(self.solvers{e}) && isprop(self.solvers{e},'options')
                self.solvers{e}.options.timespan = ts;
            end
        end

        function assertOneStageHorizon_(self, lang)
            % ASSERTONESTAGEHORIZON_(LANG)  One horizon crosses the wire.
            %
            % Both bridges send a SINGLE stage timespan, taken from solvers{1},
            % because the ensemble is built from one solver factory. The
            % resolved horizon is per stage (30/minrate of THAT stage's model),
            % so stages whose slowest rate differs resolve differently and only
            % the first would cross. Say so rather than let the remote engine
            % integrate a horizon nobody chose for it.
            ends = zeros(1, numel(self.solvers));
            for e = 1:numel(self.solvers)
                ts = self.stageHorizon_(e);
                if numel(ts) < 2
                    return
                end
                ends(e) = ts(2);
            end
            if numel(unique(ends)) > 1
                line_warning(mfilename, ...
                    ['Stage horizons differ (%s) but lang=''%s'' forwards only the ' ...
                    'first: the remote engine will integrate every stage to %g. ' ...
                    'Give the stage solvers one explicit ''timespan'' to remove ' ...
                    'the ambiguity.\n'], mat2str(ends), lang, ends(1));
            end
        end

        function printBanner_(self, runtime)
            % PRINTBANNER_(RUNTIME)  The completion line for an environment solve.
            %
            % Called at every exit of getEnsembleAvg that actually solved --
            % the native fixed point, the two bridges and the early return on
            % an empty result all leave by different paths, and a banner
            % printed on only some of them is worse than none: it labels some
            % tables and leaves others anonymous.
            if ~(isfield(self.options,'verbose') && self.options.verbose)
                return
            end
            method = 'default';
            if isfield(self.options,'method') && ~isempty(self.options.method)
                method = char(self.options.method);
            end
            lang = 'matlab';
            if isfield(self.options,'lang') && ~isempty(self.options.lang)
                lang = char(self.options.lang);
            end
            iter = 0;
            if ~isempty(self.results)
                iter = size(self.results,1);
            end
            % line_printf, as NetworkSolver.setAvgResults does, so a session-wide
            % SILENT silences this banner too. That is only safe because
            % runAnalyzerChecks now scopes its write to GlobalConstants.Verbose
            % to the analysis (GlobalConstants.pushVerbose): it used to leak,
            % so a stage solver built with verbose=0 -- which
            % renv_threestages_repairmen does -- left the whole session silent
            % and would have suppressed the ensemble's own banner.
            % with the console on this line is held until after DONE
            LineConsole.deferPrint(['ENV analysis [method: %s; type: approximate, deterministic; ' ...
                'lang: %s; env: %s] completed in %fs. Iterations: %d.\n'], ...
                method, lang, version('-release'), runtime, iter);
        end

        function solveEnvLimit(self)
            % SOLVEENVLIMIT()  Closed-form fast/slow random-environment limits.
            %
            % These treat the environment (stage) process as either infinitely
            % fast or infinitely slow relative to the base-model dynamics, and
            % therefore require no inter-stage coupling iteration:
            %
            %   'avg' (fast-environment limit): the base model sees the
            %         stationary-probability-weighted average of the modulated
            %         rates. A single rate-averaged model is built and solved
            %         once. Exact as the stage-switching rate -> Inf.
            %
            %   'dec' (slow-environment / quasi-stationary decomposition): each
            %         stage is solved independently in steady state and the
            %         per-stage metrics are averaged with weights probEnv(e).
            %         Exact as the stage-switching rate -> 0.
            %
            % Both populate self.result.Avg.{Q,U,T}; for stateful Cache nodes
            % the aggregated hit/miss ratios are written back onto the stage-1
            % reference model (self.ensemble{1}), so
            % self.ensemble{1}.getNodeByName('Cache').getHitRatio() returns the
            % environment-aggregated value.
            self.init();
            E = self.getNumberOfModels;
            probEnv = self.envObj.probEnv(:);
            M = self.ensemble{1}.getNumberOfStations;
            K = self.ensemble{1}.getNumberOfClasses;
            method = lower(self.options.method);

            Qval = zeros(M,K); Uval = zeros(M,K); Tval = zeros(M,K);
            nnodes = length(self.ensemble{1}.nodes);
            cacheIdx = find(cellisa(self.ensemble{1}.nodes, 'Cache'))';
            cacheHit  = cell(1,nnodes);
            cacheMiss = cell(1,nnodes);
            cacheHitL = cell(1,nnodes);

            switch method
                case 'dec'
                    for e = 1:E
                        se = self.solvers{e};
                        se.reset();
                        [Qe,Ue,~,Te] = se.getAvg();
                        Qval = Qval + probEnv(e)*Qe;
                        Uval = Uval + probEnv(e)*Ue;
                        Tval = Tval + probEnv(e)*Te;
                        se.getAvgNodeTable(); % populate cache metrics on stage e
                        for c = cacheIdx
                            cacheNode = self.ensemble{e}.getNodeByIndex(c);
                            [cacheHit, cacheMiss, cacheHitL] = SolverENV.accumCacheMetric(...
                                cacheHit, cacheMiss, cacheHitL, c, cacheNode, probEnv(e));
                        end
                    end
                case 'avg'
                    avgModel = self.buildRateAveragedModel(probEnv);
                    innerSolver = feval(class(self.solvers{1}), avgModel, self.solvers{1}.options);
                    [Qval,Uval,~,Tval] = innerSolver.getAvg();
                    innerSolver.getAvgNodeTable(); % populate cache metrics
                    for c = cacheIdx
                        cacheNode = avgModel.getNodeByIndex(c);
                        [cacheHit, cacheMiss, cacheHitL] = SolverENV.accumCacheMetric(...
                            cacheHit, cacheMiss, cacheHitL, c, cacheNode, 1.0);
                    end
                otherwise
                    line_error(mfilename, sprintf('solveEnvLimit called with unsupported method %s.', method));
            end

            % Write aggregated cache metrics onto the stage-1 reference model
            for c = cacheIdx
                refCache = self.ensemble{1}.getNodeByIndex(c);
                if ~isempty(cacheHit{c});  refCache.setResultHitProb(cacheHit{c});   end
                if ~isempty(cacheMiss{c}); refCache.setResultMissProb(cacheMiss{c}); end
                if ~isempty(cacheHitL{c}); refCache.setResultHitProbList(cacheHitL{c}); end
            end

            self.result.Avg.Q = Qval;
            self.result.Avg.U = Uval;
            self.result.Avg.T = Tval;
        end

        function avgModel = buildRateAveragedModel(self, probEnv)
            % AVGMODEL = BUILDRATEAVERAGEDMODEL(PROBENV)  Fast-environment model.
            % Replace every environment-modulated (i.e. stage-varying) station
            % rate by its probEnv-weighted average, represented as an
            % exponential. Non-modulated parameters keep their original
            % distribution, so the base model is preserved exactly outside the
            % modulated rates.
            E = self.getNumberOfModels;
            avgModel = self.ensemble{1}.copy();
            M = avgModel.getNumberOfStations;
            K = avgModel.getNumberOfClasses;
            for i = 1:M
                node = avgModel.stations{i};
                if isa(node,'Cache') || isa(node,'Sink')
                    continue % stateful/absorbing nodes carry no service rate
                end
                for k = 1:K
                    r = zeros(1,E);
                    ok = true;
                    for e = 1:E
                        r(e) = self.sn{e}.rates(i,k);
                    end
                    if any(isnan(r)) || any(r<=0)
                        continue % disabled for some stage: leave as configured
                    end
                    if (max(r)-min(r)) <= 1e-12*max(1,max(r))
                        continue % not modulated: keep the original distribution
                    end
                    ravg = probEnv(:)'*r(:);
                    if isa(node,'Source')
                        node.setArrival(avgModel.classes{k}, Exp(ravg));
                    elseif isa(node,'Queue') || isa(node,'Delay')
                        node.setService(avgModel.classes{k}, Exp(ravg));
                    end
                end
            end
            % The copy inherited a cached NetworkStruct from the stage model;
            % force a hard rebuild so the averaged rates take effect.
            avgModel.refreshStruct(true);
        end

        function varargout = getAvgTable(self, varargin)
            % [AVGTABLE,QT,UT,TT] = GETAVGTABLE(SELF,KEEPDISABLED)
            % Return table of average station metrics
            % The result recorder captures the returned table together with the solver
            % that produced it -- see LineResultRecorder. Recording an ensemble here
            % rather than in the member solver it delegates to is what keeps an
            % AUTO/LN/ENV/UQ answer from being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getAvgTable_impl(varargin{:});
            LineResultRecorder.capture(scope, self, 'avg', varargout{1});
        end

        function [AvgTable,QT,UT,TT] = getAvgTable_impl(self,keepDisabled)
            % GETAVGTABLE_IMPL Implementation of GETAVGTABLE; see the wrapper above.

            if nargin<2 %if ~exist('keepDisabled','var')
                keepDisabled = false;
            end

            [QN,UN,~,TN] = getAvg(self);
            M = size(QN,1);
            K = size(QN,2);
            Q = self.result.Avg.Q;
            U = self.result.Avg.U;
            T = self.result.Avg.T;
            % Aggregate stage-1 station/class labels (from the layer networks
            % for LQN stages). see _kb/06-solver-catalog.md for rationale
            [stationLabels, classLabels] = self.aggregateStationClassNames(M, K);
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
                            JobClass{end+1,1} = classLabels{k};
                            Station{end+1,1} = stationLabels{ist};
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
            % List valid methods for this solver.
            %
            % Every name here selects a COUPLING, i.e. what is carried across an
            % environment switch, and each is dispatched somewhere in this class
            % or its analyzers. The list used to name only 'default', which made
            % the other six reachable only by knowing the source:
            %   'default'/'meanfield'  the mean-field coupling: the marginal
            %                          means cross a switch (analyzerMode
            %                          'meanfield', solver_env_meanfield_analyzer)
            %   'statevec'/'blend'     the state-vector coupling: the whole joint
            %                          distribution crosses (analyzerMode
            %                          'statevec', solver_env_statevec_analyzer)
            %   'smp'                  accepted for a semi-Markov environment,
            %                          but Environment.init still reads every
            %                          arc as a (D0,D1) process, so it runs the
            %                          mean-field coupling on Markovian arcs
            %                          only (see methodRefusal)
            %   'statedep'             the environment transition depends on the
            %                          state it leaves (meanfield analyzer:264)
            %   'avg'/'dec'            the closed-form fast/slow environment
            %                          limits, which replace the fixed point
            %                          rather than configure it (solveEnvLimit)
            % No getStruct here: the model is an Environment, which has none,
            % and the call made every probe of this family raise, so model.help
            % on an Environment listed no env row at all.
            allMethods = {'default','meanfield','smp','statedep','statevec','blend','avg','dec'};
        end

        function [bool, reason] = supportsModelMethod(self, method)
            % [BOOL, REASON] = SUPPORTSMODELMETHOD(METHOD)
            % The method-aware gate model.help and SolverAUTO ask; the rules
            % and their second caller, init, are in METHODREFUSAL.
            [bool, reason] = self.methodRefusal(method);
        end

        function [ok, reason] = methodRefusal(self, method)
            % [OK, REASON] = METHODREFUSAL(METHOD)
            % Whether METHOD can run on this environment with these stage
            % solvers, as a predicate. ONE PREDICATE, TWO CALLERS: the gate
            % above and INIT, which every native solve passes through, so the
            % run raises the sentence the report showed rather than an index
            % error inside Environment.init or a table of zeros.
            %
            % THE RULES. (1) Every arc is Markovian. Environment.init reads each
            % transition as a (D0,D1) process to build the holding-time MMAPs
            % and probEnv; a Det, Gamma or Uniform arc has no such form, its
            % getProcess returns bare parameters, and init indexes them as a
            % MAP (an index error for Det, garbage for the others). 'smp' used
            % to lift only the constructor's check and reach that same init, so
            % it is gated alike until a semi-Markov analyzer exists.
            % (2) The mean-field family ('default', 'meanfield', 'smp',
            % 'statedep') integrates each stage over its sojourn, so every
            % stage solver must return transient averages: analyze_ swallows a
            % stage whose getTranAvg fails and the table came back all zero.
            % This is the rule mapEnvApprox applies to itself. (3) 'statevec'/
            % 'blend' need each stage's enumerated generator: a SolverCTMC or
            % SolverMAM stage solver with a finite timespan, and no
            % LayeredNetwork stage (solver_env_statevec_analyzer). (4) 'avg'/
            % 'dec' (solveEnvLimit) read the stations of the stage Networks, so
            % a LayeredNetwork stage has no rate-averaged model there.
            ok = true;
            reason = '';
            method = lower(char(method));
            E = size(self.env, 1); % an E x E cell of arcs, as the constructor reads it
            for e = 1:E
                for h = 1:E
                    arc = self.env{e,h};
                    if isempty(arc) || isa(arc,'Disabled') || isa(arc,'Markovian')
                        continue
                    end
                    ok = false;
                    reason = sprintf(['The distribution of the environment transition from stage ' ...
                        '%d to %d is a %s, which is not Markovian: Environment.init reads every ' ...
                        'arc as a (D0,D1) process to build the holding-time laws, and no method ' ...
                        'of %s analyses a semi-Markov environment yet, ''smp'' included. Fit the ' ...
                        'arc to a phase-type law, e.g. APH.fitMeanAndSCV.'], e, h, class(arc), self.getName);
                    return
                end
            end
            nstages = numel(self.ensemble);
            if any(strcmp(method, {'statevec','blend'}))
                for e = 1:nstages
                    if isa(self.ensemble{e}, 'LayeredNetwork')
                        ok = false;
                        reason = sprintf(['The state-vector (statevec) analyzer does not support ' ...
                            'LayeredNetwork stages (stage %d): an LQN has no single stage generator. ' ...
                            'Use the default mean-field analyzer (omit method=''statevec'').'], e);
                        return
                    end
                    s = self.solvers{e};
                    if ~(isa(s,'SolverCTMC') || isa(s,'SolverMAM'))
                        ok = false;
                        reason = sprintf(['The statevec analyzer requires a SolverCTMC or SolverMAM ' ...
                            'inner solver, but environment stage %d uses %s.'], e, class(s));
                        return
                    end
                    opts_e = s.getOptions;
                    if ~isfield(opts_e,'timespan') || numel(opts_e.timespan) < 2 ...
                            || ~isfinite(opts_e.timespan(2))
                        ok = false;
                        reason = sprintf(['The statevec analyzer requires a finite inner-solver ' ...
                            'timespan for stage %d, e.g. CTMC(model,''timespan'',[0,T]).'], e);
                        return
                    end
                end
            elseif any(strcmp(method, {'avg','dec'}))
                for e = 1:nstages
                    if isa(self.ensemble{e}, 'LayeredNetwork')
                        ok = false;
                        reason = sprintf(['The ''%s'' environment limit reads the stations of the ' ...
                            'stage models (solveEnvLimit), which a LayeredNetwork stage (stage %d) ' ...
                            'does not carry. Use the mean-field coupling (method=''default'').'], method, e);
                        return
                    end
                end
            else
                for e = 1:nstages
                    s = self.solvers{e};
                    if isa(self.ensemble{e}, 'LayeredNetwork') || isa(s, 'SolverLN')
                        continue % SolverLN.getTranAvg couples the layers itself
                    end
                    if ~(ismethod(s,'supportsTransientAnalysis') && s.supportsTransientAnalysis())
                        ok = false;
                        reason = sprintf(['The mean-field environment coupling integrates each stage ' ...
                            'over its sojourn, so it needs transient averages from the stage solver, ' ...
                            'which %s (stage %d) does not produce. Use a SolverFLD, SolverCTMC, ' ...
                            'SolverLDES or SolverJMT stage solver, or method=''dec'' or ''avg''.'], class(s), e);
                        return
                    end
                end
            end
        end

        function [stationLabels, classLabels] = aggregateStationClassNames(self, M, K)
            % [STATIONLABELS, CLASSLABELS] = AGGREGATESTATIONCLASSNAMES(M,K)
            % Aggregate station/class label vectors for the stage-1 result
            % layout. For flat NetworkStruct stages these are the usual node
            % and class names; for LQN stages the layout is the block-diagonal
            % union of the layer networks, so labels are drawn from each layer
            % network (prefixed by the layer name to disambiguate).
            stationLabels = cell(M,1);
            classLabels = cell(K,1);
            if isfield(self.sn{1}, 'nstations')
                for ist=1:M
                    stationLabels{ist} = self.sn{1}.nodenames{self.sn{1}.stationToNode(ist)};
                end
                for k=1:K
                    classLabels{k} = self.sn{1}.classnames{k};
                end
            elseif isa(self.ensemble{1}, 'LayeredNetwork')
                model = self.ensemble{1};
                [Roff, Coff, Msz, Ksz] = model.layerBlocks;
                for e=1:length(model.ensemble)
                    layer = model.ensemble{e};
                    lname = layer.getName;
                    for i=1:Msz(e)
                        stationLabels{Roff(e)+i} = sprintf('%s.%s', lname, layer.stations{i}.name);
                    end
                    for r=1:Ksz(e)
                        classLabels{Coff(e)+r} = sprintf('%s.%s', lname, layer.classes{r}.name);
                    end
                end
            else
                for ist=1:M
                    stationLabels{ist} = sprintf('Station%d', ist);
                end
                for k=1:K
                    classLabels{k} = sprintf('Class%d', k);
                end
            end
        end

        function Q = roundMarginalForDiscreteSolver(self, Q, snRef)
            % ROUNDMARGINALFORDISCRETESOLVER Round fractional queue lengths
            % to integers using the largest remainder method, preserving
            % closed chain populations exactly.
            for c = 1:snRef.nchains
                chainClasses = find(snRef.chains(c,:) > 0);
                njobs_chain = sum(snRef.njobs(chainClasses));
                if isinf(njobs_chain)
                    % Open chain: simple rounding
                    for idx = 1:length(chainClasses)
                        k = chainClasses(idx);
                        for i = 1:size(Q,1)
                            Q(i,k) = round(Q(i,k));
                        end
                    end
                else
                    % Closed chain: largest remainder method
                    vals = [];
                    indices = [];
                    for i = 1:size(Q,1)
                        for idx = 1:length(chainClasses)
                            k = chainClasses(idx);
                            vals(end+1) = Q(i,k);
                            indices(end+1,:) = [i, k];
                        end
                    end
                    floored = floor(vals);
                    remainders = vals - floored;
                    deficit = njobs_chain - sum(floored);
                    deficit = round(deficit); % ensure integer
                    if deficit > 0
                        [~, sortIdx] = sort(remainders, 'descend');
                        for d = 1:min(deficit, length(sortIdx))
                            floored(sortIdx(d)) = floored(sortIdx(d)) + 1;
                        end
                    end
                    for j = 1:length(vals)
                        Q(indices(j,1), indices(j,2)) = floored(j);
                    end
                end
            end
        end
    end

    methods (Static)

        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)
            % Whether MODEL is an environment this solver can be built over: an
            % Environment with at least one stage, each a Network or a
            % LayeredNetwork. What a stage USES is the stage solver's to accept,
            % and the constructor asks every stage solver so; the method-level
            % rules (Markovian arcs, stage solvers that return transients or a
            % generator) are METHODREFUSAL's, which needs the built solver.
            %
            % This used to compare model.getUsedLangFeatures() against a flat
            % set of node and discipline names, but an Environment has no such
            % method: the call raised, findSolver swallowed the error as "no
            % claim, no gate", and every env row read as runnable. The flat set
            % also named no Cache, which cache_mmap_rr_env solves through every
            % method here, so it was not a description of the analyzers either.
            featSupported = SolverFeatureSet;
            bool = isa(model, 'Environment');
            if ~bool
                return
            end
            stages = model.getEnsemble;
            bool = ~isempty(stages);
            for e = 1:numel(stages)
                bool = bool && (isa(stages{e},'Network') || isa(stages{e},'LayeredNetwork'));
            end
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

        function [cacheHit, cacheMiss, cacheHitL] = accumCacheMetric(cacheHit, cacheMiss, cacheHitL, c, cacheNode, w)
            % Accumulate w-weighted cache hit/miss ratios of node index c
            % (used by solveEnvLimit for the 'avg'/'dec' environment limits).
            h = full(cacheNode.getHitRatio);
            if ~isempty(h)
                if isempty(cacheHit{c}); cacheHit{c} = zeros(size(h)); end
                cacheHit{c} = cacheHit{c} + w*h;
            end
            mval = full(cacheNode.getMissRatio);
            if ~isempty(mval)
                if isempty(cacheMiss{c}); cacheMiss{c} = zeros(size(mval)); end
                cacheMiss{c} = cacheMiss{c} + w*mval;
            end
            hl = cacheNode.getHitRatioByList;
            if ~isempty(hl)
                hl = full(hl);
                if isempty(cacheHitL{c}); cacheHitL{c} = zeros(size(hl)); end
                cacheHitL{c} = cacheHitL{c} + w*hl;
            end
        end
    end
end

