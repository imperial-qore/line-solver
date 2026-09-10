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
% np_priority selects the approximation used for the priority scaling at HOL
% stations: 'cl'/'chandy-lakshmi' (default) is Chandy-Lakshmi [ChaL83] in the
% arrival-instant utilization form of Eager-Lipscomb [EagL88]; 'shadow' is
% Sevcik's shadow server [Sev77]. Note [ChaL83] derives the approximation for
% PREEMPTIVE priority; LINE applies it at non-preemptive (HOL) stations, which
% is the setting Eager and Lipscomb study.
options.config.np_priority = 'default';
options.config.fork_join = 'default';
% see _kb/05-solvers-overview.md for rationale
options.config.fj_warmstart = true;
% Random-environment fallback for MAP/MMPP/MMAP models on solvers that cannot
% consume a non-renewal process: 'auto' approximates, 'off' rejects the model.
% see _kb/05-solvers-overview.md for rationale
options.config.map_env = 'auto';
options.config.map_env_method = 'auto'; % environment coupling: 'auto', 'meanfield' (transient-capable solvers), 'dec' (slow), 'avg' (fast)
options.config.map_env_maxstages = 64;  % cap on the product of the phase orders
options.config.nonmkv = 'bernstein'; % Method for non-Markovian distribution conversion: 'none', 'bernstein'
options.config.nonmkvorder = 20; % Order (number of phases) for non-Markovian distribution approximation
% see _kb/05-solvers-overview.md for rationale
options.config.symbolic = 'auto';
options.config.symbolic_timeout = 300; % seconds per symbolic request
global LINEDefaultLang
if isempty(LINEDefaultLang)
    % A test script running `clear all`/`clear java` deletes the global; the
    % environment variable survives it and re-seeds the dispatch backend.
    LINEDefaultLang = getenv('LINE_DEFAULT_LANG');
end
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
% THE STIFF ACCURATE SLOT. This is what the fluid solver reaches by default
% (options.stiff is true), and it is LSODA in all four ports: the JAR, native
% python and C++ integrate this slot with LSODA and the vendored
% lsoda_accurate_stiff is their MATLAB counterpart (order budget 12/5, BDF half
% pinned), so the four now agree in the last digits rather than by luck.
%
% IT IS LSODA BECAUSE @ode15s DOES NOT TERMINATE on a stiff LN layer.
% line-test.git testsLQNS/11-interlock.lqnx -- three tasks, twelve elements --
% carries an InfRate 1e8 in ALL FIVE of its layers (its sibling 10-interlock
% carries one). A single window of layer 1 never returns under @ode15s: 75 min
% at 4.8 GB and still growing, with the LN loop capped at iter_max = 5 and
% interlocking off, so it is one integration and not the outer iteration. That
% one model held the whole MATLAB suite twice -- the 2026-08-27 run was killed
% after 100 minutes on it and the 2026-08-28 20:14 run after 50 -- taking every
% block after it down with it. Under @lsoda_accurate_stiff the model returns in
% 42.4 s and all seven interlock models pass against their goldens. On
% 10-interlock, where @ode15s does finish, the two agree to every digit printed
% (Util 0.226604 / 0.453208 across all twelve rows) and LSODA is faster, 17.7 s
% against 28.4 s. See BUG-99.
%
% TWO THINGS HAD TO BE TRUE FIRST, and both are now:
%
%   1. A JUMP IN THE DRIFT IS AN INTEGRATION BOUNDARY. An NHPP intensity enters
%      the drift through fluid_interpcols(rt_tgrid, rt_Mmat, t), so the
%      right-hand side is PIECEWISE CONSTANT WITH JUMPS at the schedule's
%      breakpoints. @ode15s's step control happened to cope; LSODA stepped over
%      a segment and test_solver_fld_nhpp failed its `maxrel < 0.05` on the
%      queue throughput tracking lambda(t). SOLVER_FLUID_RATEMULT now reports
%      those instants and SOLVER_FLUID_ITERATION restarts the integrator at each
%      of them, which is the right thing under any integrator and costs nothing
%      on a model that has none.
%   2. THE CONSERVATION GUARD SURVIVES THE SWAP. FLUID_CONSERVATION_GUARD halts
%      a moment-closure excursion through odeset('OutputFcn'), and
%      lsoda_odesolve used to ignore that field -- so moving this slot would
%      have silently disarmed the guard and handed back the CQN_Cox_CS_9 hang it
%      was written for. lsoda_odesolve now implements the init/step/done
%      protocol and truncates on a nonzero return.
odesfun.accurateStiffOdeSolver = @lsoda_accurate_stiff;
%odesfun.accurateStiffOdeSolver = @ode15s;
% The index-1 DAE slot, read only by the fluid 'dae' method. It is separate
% from accurateStiffOdeSolver -- which is also @ode15s today -- because the two
% have different contracts: the ODE paths set NonNegative on their odeset, and
% @rodas has no hook for holding a component at or above zero and refuses that
% field rather than ignoring it. Set this to @rodas to run the vendored
% Hairer-Wanner RODAS, the same integrator the C++, native Python and JAR ports
% of this method use, so that all four agree in the last digits; @ode15s stays
% the default so an existing MATLAB answer does not move underneath a caller.
odesfun.daeSolver = @ode15s;
%odesfun.daeSolver = @rodas;
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
% INHERITED FROM THE SESSION, not hardcoded to STD, exactly as the JAR does
% (`this.verbose = GlobalConstants.getVerbose()`). A hardcoded STD here meant a
% run's own options ALWAYS named a level, so the session level was never
% consulted and line_verbosity(VerboseLevel.DEBUG) reached no solver at all.
options.verbose = GlobalConstants.Verbose;
if isempty(options.verbose)
    options.verbose = VerboseLevel.STD; % lineStart has not run yet
end
% The solver console (LineConsole) has NO option of its own: it IS
% VerboseLevel.DEBUG, so 'verbose',VerboseLevel.DEBUG is what asks for it.

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
        % Verbosity is INHERITED from the session, like every other solver
        % (see the top of this file). Pinning STD here silently downgraded a
        % session that had asked for VerboseLevel.DEBUG, so this ensemble was
        % one of the few solvers the console could never narrate. Pass
        % 'verbose',VerboseLevel.SILENT explicitly for a quiet run.
        options.config.da = 'courtois'; % CTMC decomposition/aggregation: 'courtois', 'kms', 'takahashi', 'multi'
        options.config.da_iter = 10; % Number of iterations for kms/takahashi
    case {'FLD','Fluid'}
        options = Solver.defaultOptions();
        options.config.highvar = 'default';
        % TRUE: a coordinate whose exit rate is GlobalConstants.Immediate is
        % LINE's stand-in for infinity, not a fast rate, and every fluid route
        % now stochastic-complements it out of the event set rather than
        % integrating it (ODE_ELIMINATE_IMMEDIATE, reached through
        % FLUID_HIDE_IMMEDIATE). It read false while this slot was @ode15s,
        % which coped; under @lsoda_accurate_stiff a single LN task layer
        % carrying one InfRate phase took 70 s in MATLAB and 145 s in the JAR
        % for the answer the reduced system returns in 0.35 s.
        options.config.hide_immediate = true;
        % Stop the ODE window loop once the drift residual says the state is at
        % a fixed point, instead of always running iter_max windows.
        options.config.fluid_earlystop = true;
        options.iter_max = 200;
        options.stiff = true;
        options.timespan = [0,Inf];
    case 'JMT'
        % use default
    case 'LN'
        options = Solver.defaultOptions();
        options.config.interlocking = true;
        options.config.multiserver = 'default';
        % Layering strategy: 'srvn' places each server in its own submodel,
        % 'flat' places every processor and task in a single submodel
        options.config.layering = 'srvn';
        % Under-relaxation options for convergence improvement
        options.config.relax = 'fixed';      % 'auto' | 'fixed' | 'adaptive' | 'none'
        options.config.relax_factor = 0.5;   % Relaxation factor (0 < omega <= 1)
        options.config.relax_min = 0.1;      % Minimum relaxation factor for adaptive mode
        options.config.relax_history = 5;    % Error history window for adaptive mode
        % Warm start of a fluid layer solver from the terminal ODE state of
        % its previous outer iteration, instead of the cold placement of
        % SOLVER_FLUID_INITSOL. Set false to restart every layer solve cold.
        options.config.fluid_warmstart = true;
        % Stochastic iteration options, used when layer solvers are
        % simulation-based (JMT, SSA, LDES) and thus return noisy estimates
        options.config.stochiter = 'auto';        % 'auto' | 'rm' | 'crn' | 'off'
        options.config.stochiter_alpha = 0.6;     % Robbins-Monro step decay exponent, in (0.5,1]
        options.config.stochiter_a0 = 1.0;        % Robbins-Monro initial step after burn-in
        options.config.stochiter_burnin = 5;      % Picard burn-in iterations before step decay starts
        options.config.stochiter_conseq = 3;      % consecutive sub-tolerance iterations required to stop
        options.timespan = [Inf,Inf];
        options.keep = true;
        % verbosity is INHERITED from the session, see the top of this file
        options.iter_max = 200;              % More iterations for difficult LQN models
        options.iter_tol = 5e-3;             % Convergence tolerance (looser than default for LQN models)
        options.tol = 1e-4;
        % MOL (Method of Layers) options for hierarchical iteration
    case 'LQNS'
        options = Solver.defaultOptions();
        options.timespan = [Inf,Inf];
        options.keep = true;
        % Verbosity is INHERITED from the session, like every other solver
        % (see the top of this file); the lqns binary's own advisories stay
        % suppressed until DEBUG, see runAnalyzer
        options.config.multiserver = 'default';
        options.config.remote = false;  % Enable remote execution via REST API
        options.config.remote_url = 'http://localhost:8080';  % URL of lqns-rest server
        options.config.container = '';  % JMT only: explicit Docker image ('' = default imperialqore/jmt-rest, or LINE_JMT_IMAGE); LQNS/lqsim/qnsolver ignore it, need a local binary
    case 'AG'
        options.iter_max = 100;
        options.timespan = [Inf,Inf];
        % Truncation level of an OPEN component's queue-length dimension. A
        % closed class uses its own population instead, so this bounds only the
        % open agents; 'inapinf' ignores it and solves them on the infinite
        % state space through the matrix-geometric tail.
        options.config.maxStates = 100;
        % Execution backend of the reversed-rate fixed point. Every agent is
        % solved in isolation given the current reversed rates, so the sweep is
        % order-free and each backend produces the SAME iterates:
        %   'serial'  - one agent after another, in index order (reference)
        %   'parallel' - agents fanned out over a parfor pool, barrier per sweep
        %              (alias 'para')
        %   'cluster' - agents partitioned over remote workers, barrier per sweep
        % 'parallel' is asserted bit-identical to 'serial'; 'cluster' is
        % bit-identical only against a worker running the same implementation,
        % and agrees to a few ulp against another codebase's worker. See
        % _kb/06-solver-catalog.md.
        options.config.exec = 'serial';
        % 'parallel': size of the pool, 0 = let MATLAB choose the default pool.
        options.config.nworkers = 0;
        % 'cluster': worker endpoints as a cell array of 'host:port' strings.
        options.config.endpoints = {};
        % 'cluster': seconds to wait on a worker before solving its agents
        % locally instead. A lost worker is never fatal, because any agent can
        % be solved anywhere given the reversed rates.
        options.config.worker_timeout = 30;
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
        % verbosity is INHERITED from the session, see the top of this file
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
        % see _kb/05-solvers-overview.md for rationale (LDES and JMT REST endpoint)
        options.rest_url = '';
end
end
