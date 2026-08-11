function options = SolverOptions(solverName)
% SOLVEROPTIONS Create solver configuration options structure
%
% @brief Creates a configuration structure with solver-specific default options
% @param solverName Optional solver name for specific configurations (default: 'Solver')
% @return options Struct containing solver configuration parameters
%
% SolverOptions generates a standardized options structure for configuring
% LINE solvers. It provides default values for common parameters like
% convergence tolerances, iteration limits, ODE solvers, and solver-specific
% settings. The function customizes defaults based on the solver type.
%
% Common options include:
% - Convergence parameters (tol, iter_tol, iter_max)
% - Analysis parameters (samples, seed, cutoff)
% - Language selection (MATLAB vs Java)
% - ODE solver configuration for fluid methods
% - Verbosity and caching controls
% - Solver-specific configuration overrides
%
% Solver-specific customizations are available for:
% - CTMC: State space generation and transient analysis
% - Fluid: ODE solver selection and timespan
% - JMT: Simulation parameters and confidence intervals
% - MVA: Method selection and approximation settings
% - SSA: Sampling and parallel execution options
%
% Example:
% @code
% opts = SolverOptions('MVA');       % MVA-specific defaults
% opts.method = 'exact';             % Override method
% opts.iter_tol = 1e-6;             % Tighter tolerance
% solver = SolverMVA(model, opts);   % Use custom options
% @endcode

if nargin < 1
    solverName = 'Solver'; % global options unless overridden by a solver
end

%% Solver default options
options = struct();
options.cache = true;
options.cutoff = Inf;
options.config = struct(); % solver specific options
options.config.highvar = 'default';
options.config.multiserver = 'default';
options.config.np_priority = 'default';
options.config.fork_join = 'default';
% see _kb/05-solvers-overview.md for rationale
options.config.fj_warmstart = true;
options.config.nonmkv = 'bernstein'; % Method for non-Markovian distribution conversion: 'none', 'bernstein'
options.config.nonmkvorder = 20; % Order (number of phases) for non-Markovian distribution approximation
% see _kb/05-solvers-overview.md for rationale
options.config.symbolic = 'auto';
options.config.symbolic_timeout = 300; % seconds per symbolic request
global LINEDefaultLang
if ~isempty(LINEDefaultLang)
    options.lang = LINEDefaultLang;
else
    options.lang = 'matlab';
end
options.force = false;
options.init_sol = [];
options.iter_max = 1000;
options.iter_tol = 1e-4; % convergence tolerance to stop iterations
options.tol = 1e-4; % tolerance for all other uses
options.keep = true;
options.method = 'default';
%options.remote = false;
%options.remote_endpoint = '127.0.0.1';

odesfun = struct();
odesfun.fastOdeSolver = @ode23;
%odesfun.fastOdeSolver = @lsoda_fast;
odesfun.accurateOdeSolver = @ode113;
%odesfun.accurateOdeSolver = @lsoda_accurate;
odesfun.fastStiffOdeSolver = @ode23s;
%odesfun.fastStiffOdeSolver = @lsoda_fast_stiff;
odesfun.accurateStiffOdeSolver = @ode15s;
%odesfun.accurateStiffOdeSolver = @lsoda_accurate_stiff;
options.odesolvers = odesfun;

% see _kb/05-solvers-overview.md for rationale (samples vs events semantics)
options.samples = 1e4;
% see _kb/05-solvers-overview.md for rationale (samples vs events semantics)
options.events = NaN;
options.seed = randi([1,1e6]);
options.stiff = true;
options.confint = false; % confidence interval: false, true (95%), or level (0.0-1.0)
options.timespan = [Inf,Inf];
options.timestep = [];
% see _kb/05-solvers-overview.md for rationale (wall-clock, not simulated time)
options.timeout = Inf;
options.verbose = VerboseLevel.STD;

options.config.num_cdf_pts = 200;

%% Solver-specific defaults
switch solverName
    case 'CTMC'
        options.cutoff = 10; % finite per-class state-space cutoff for open/mixed models (matches native/JAR)
        options.timespan = [Inf,Inf];
        options.timestep = []; % timestep for fixed time steps in transient analysis
        options.config.hide_immediate = true; % hide immediate transitions if possible
        %options.config.state_space_gen = 'reachable'; % still buggy
        options.config.state_space_gen = 'full';
        options.rewardIterations = 1000; % number of value iterations for reward computation
        options.config.qrf_params = [];   % blocking config struct (f, MR, BB, F, MM, MM1, ZZ, ZM)
        options.config.qrf_alpha = [];    % load-dependent rates M x N matrix
    case {'ENV','Env'}
        options.method = 'default';
        options.init_sol = [];
        options.iter_max = 100;
        options.iter_tol = 1e-4;
        options.tol = 1e-4;
        options.verbose = VerboseLevel.SILENT;
        options.config.da = 'courtois'; % CTMC decomposition/aggregation: 'courtois', 'kms', 'takahashi', 'multi'
        options.config.da_iter = 10; % Number of iterations for kms/takahashi
    case {'FLD','Fluid'}
        options = Solver.defaultOptions();
        options.config.highvar = 'default';
        options.config.hide_immediate = false; % stiff ODE solver handles immediate rates accurately
        options.iter_max = 200;
        options.stiff = true;
        options.timespan = [0,Inf];
    case 'JMT'
        % use default
    case 'LN'
        options = Solver.defaultOptions();
        options.config.interlocking = true;
        options.config.multiserver = 'default';
        % Under-relaxation options for convergence improvement
        options.config.relax = 'fixed';      % 'auto' | 'fixed' | 'adaptive' | 'none'
        options.config.relax_factor = 0.5;   % Relaxation factor (0 < omega <= 1)
        options.config.relax_min = 0.1;      % Minimum relaxation factor for adaptive mode
        options.config.relax_history = 5;    % Error history window for adaptive mode
        % Stochastic iteration options, used when layer solvers are
        % simulation-based (JMT, SSA, LDES) and thus return noisy estimates
        options.config.stochiter = 'auto';        % 'auto' | 'rm' | 'crn' | 'off'
        options.config.stochiter_alpha = 0.6;     % Robbins-Monro step decay exponent, in (0.5,1]
        options.config.stochiter_a0 = 1.0;        % Robbins-Monro initial step after burn-in
        options.config.stochiter_burnin = 5;      % Picard burn-in iterations before step decay starts
        options.config.stochiter_conseq = 3;      % consecutive sub-tolerance iterations required to stop
        options.timespan = [Inf,Inf];
        options.keep = true;
        options.verbose = VerboseLevel.STD;
        options.iter_max = 200;              % More iterations for difficult LQN models
        options.iter_tol = 5e-3;             % Convergence tolerance (looser than default for LQN models)
        options.tol = 1e-4;
        % MOL (Method of Layers) options for hierarchical iteration
    case 'LQNS'
        options = Solver.defaultOptions();
        options.timespan = [Inf,Inf];
        options.keep = true;
        options.verbose = false;
        options.config.multiserver = 'default';
        options.config.remote = false;  % Enable remote execution via REST API
        options.config.remote_url = 'http://localhost:8080';  % URL of lqns-rest server
        options.config.container = '';  % Run lqns/lqsim inside Docker: '' = native (docker only as fallback if native missing), 'auto'/true = prefer docker if available, or an explicit image name (e.g. 'imperialqore/line-lqns-rest:latest')
    case 'MAM'
        options.iter_max = 100;
        options.timespan = [Inf,Inf];
        % FJ-specific options (used when Fork-Join topology detected)
        options.config.fj_accuracy = 100;    % C parameter for FJ_codes (higher = more accurate)
        options.config.fj_tmode = 'NARE';    % T matrix computation: 'NARE' or 'Sylves'
        % num_cdf_pts uses global default of 200
    case 'MVA'
        options.iter_max = 1000;
        options.iter_tol = 1e-6;
    case 'NC'
        options.samples = 1e5;
        options.timespan = [Inf,Inf];
        options.config.highvar = 'interp';
    case 'NN'
        options.iter_max = 1000;
        options.iter_tol = 1e-6;
    case 'QNS'
        options.config.multiserver = 'default';
    case 'SSA'
        options.timespan = [0,Inf];
        options.verbose = true;
        options.config.state_space_gen = 'none';
        % see _kb/05-solvers-overview.md for rationale (SSA warmup discard)
        options.config.warmupfrac = 0.0;
        % see _kb/05-solvers-overview.md for rationale (worker-count-invariant replications)
        options.config.nreplicas = 8;
        switch options.lang
            case 'java'
                options.config.eventcache = true;
            otherwise
                options.config.eventcache = false;
        end
    case 'LDES'
        options.samples = 2e5;
        options.lang = 'java';
        % Transient detection options
        options.config.tranfilter = 'mser5';  % 'mser5', 'fixed', or 'none'
        options.config.mserbatch = 5;         % MSER batch size (default: 5)
        options.config.warmupfrac = 0.2;      % Warmup fraction for fixed filter (0.0 to 1.0)
        % Confidence interval options
        options.config.cimethod = 'obm';      % 'obm', 'bm', 'spectral', or 'none'
        options.config.obmoverlap = 0.5;      % OBM overlap fraction (0.0 to 1.0)
        options.config.ciminbatch = 10;       % Minimum batch size for CI
        options.config.ciminobs = 100;        % Minimum observations for CI
        options.config.spectrallowfreqfrac = 0.25; % Low-frequency fraction, cimethod='spectral'
        % Parallel replication options. Above 1, each replication runs with its
        % own RNG stream and the confidence intervals use the cross-replication
        % variance instead of batch means.
        options.config.replications = 1;    % Number of independent replications
        options.config.numthreads = [];     % Worker threads (empty = engine default)
        % Convergence options
        options.config.cnvgon = false;      % Enable convergence-based stopping
        options.config.cnvgtol = 0.05;      % Convergence tolerance (5% relative precision)
        options.config.cnvgbatch = 20;      % Min batches before checking convergence
        options.config.cnvgchk = 0;         % Events between checks (0 = auto)
        % see _kb/05-solvers-overview.md for rationale (discrete-time/slotted)
        options.config.slotted = false;     % Run on a discrete time scale
        options.config.slotlength = 1;      % Slot length in model time units
        % see _kb/05-solvers-overview.md for rationale (LDES REST endpoint)
        options.rest_url = '';
end
end
