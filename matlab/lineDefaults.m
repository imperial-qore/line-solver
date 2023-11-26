function options = lineDefaults(solverName)
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
options.lang = 'matlab';
options.force = false;
options.init_sol = [];
options.iter_max = 10;
options.iter_tol = 1e-4; % convergence tolerance to stop iterations
options.tol = 1e-4; % tolerance for all other uses
options.keep = false;
options.method = 'default';
%options.remote = false;
%options.remote_endpoint = '127.0.0.1';

odesfun = struct();
odesfun.fastOdeSolver = Solver.fastOdeSolver;
odesfun.accurateOdeSolver = Solver.accurateOdeSolver;
odesfun.fastStiffOdeSolver = Solver.fastStiffOdeSolver;
odesfun.accurateStiffOdeSolver = Solver.accurateStiffOdeSolver;
options.odesolvers = odesfun;

options.samples = 1e4;
options.seed = randi([1,1e6]);
options.stiff = true;
options.timespan = [Inf,Inf];
options.verbose = 1;

%% Solver-specific defaults
switch solverName
    case 'CTMC'
        options.timespan = [Inf,Inf];
        options.config.hide_immediate = true; % hide immediate transitions if possible
        %options.config.state_space_gen = 'reachable'; % still buggy
        options.config.state_space_gen = 'full';
    case 'Env'
        options.method = 'default';
        options.init_sol = [];
        options.iter_max = 100;
        options.iter_tol = 1e-4;
        options.tol = 1e-4;
        options.verbose = VerboseLevel.SILENT;
    case 'Fluid'
        options = Solver.defaultOptions();
        options.config.highvar = 'default';
        options.iter_max = 200;
        options.stiff = true;
        options.timespan = [0,Inf];
    case 'JMT'
        % use default
    case 'LN'
        options = Solver.defaultOptions();
        options.config.interlocking = true;
        options.config.multiserver = 'default';
        options.timespan = [Inf,Inf];
        options.keep = false;
        options.verbose = VerboseLevel.STD;
        options.iter_tol = 5e-2;
        options.iter_max = 100;
    case 'LQNS'
        options = Solver.defaultOptions();
        options.timespan = [Inf,Inf];
        options.keep = true;
        options.verbose = false;
        options.config.multiserver = 'default';
    case 'MAM'
        options.iter_max = 100;
        options.timespan = [Inf,Inf];
    case 'MVA'
        options.iter_max = 1000;
        options.iter_max = 10^3;
        options.iter_tol = 10^-6;
    case 'QNS'
        options.config.multiserver = 'default';
    case 'NC'
        options.samples = 1e5;
        options.timespan = [Inf,Inf];
        options.config.highvar = 'interp';
    case 'SSA'
        options.timespan = [0,Inf];
        options.verbose = true;
end
end
