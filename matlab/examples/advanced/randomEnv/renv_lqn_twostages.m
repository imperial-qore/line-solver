function renv_lqn_twostages()
% RENV_LQN_TWOSTAGES  LayeredNetwork operating in a two-stage random environment.
%
% Demonstrates SolverENV running over a LayeredNetwork (LQN) base model, using
% the uniform model/solver interface (no branching in ENV). The environment
% alternates between an UP stage (fast database) and a DOWN stage (slow
% database); the environment-averaged total throughput must lie between the two
% single-stage LQN solutions.

warning off;

%% Two structurally identical LQN stages, differing only in DB host demand.
upModel   = buildLQN('LQN_UP',   0.8);   % fast DB activity
downModel = buildLQN('LQN_DOWN', 3.0);   % slow DB activity (degraded)

%% Random environment: UP <-> DOWN
env = Environment('DBReliability');
env.addStage('UP',   'operational', upModel);
env.addStage('DOWN', 'degraded',    downModel);
env.addTransition('UP',   'DOWN', Exp(0.2));   % mean UP time = 5
env.addTransition('DOWN', 'UP',   Exp(1.0));   % mean DOWN time = 1
env.init();

%% Solve with SolverENV over SolverLN(.,@SolverFluid)
% The transient window is set on SolverLN (not the layer factory): the layered
% fixed-point iteration solves each layer in steady state, and SolverLN applies
% the timespan only to the per-layer transient getTranAvg call.
T = 50;
fldFactory = @(mm) SolverFluid(mm, 'verbose', false);
lnFactory  = @(m)  SolverLN(m, fldFactory, 'timespan', [0, T], 'verbose', false);

options = SolverENV.defaultOptions;
options.iter_max = 20;
options.iter_tol = 0.02;
options.verbose  = false;

envSolver = SolverENV(env, lnFactory, options);
[QN, UN, RN, TN] = envSolver.getAvg(); %#ok<ASGLU>
fprintf('ENV over LQN ran: aggregate size = %d x %d\n', size(QN,1), size(QN,2));
envAvgTable = envSolver.getAvgTable(); %#ok<NASGU>
disp(envAvgTable);

%% Single-stage aggregate references (same block-diagonal layout)
[upQ, upT]     = stageAggregate(upModel,   fldFactory, T);
[downQ, downT] = stageAggregate(downModel, fldFactory, T);

%% Verify
% (1) The solver runs and returns finite, population-conserving metrics.
assert(all(isfinite(QN(:))), 'ENV aggregate Q contains non-finite values');
assert(abs(sum(QN(:)) - sum(upQ(:))) < 1e-2 && abs(sum(QN(:)) - sum(downQ(:))) < 1e-2, ...
    'ENV aggregate does not conserve the closed population of the stages');

% (2) A physically monotone scalar (total throughput) must lie between the two
% single-stage solutions: a slower database slows the whole system, so the
% environment-averaged throughput brackets the UP and DOWN values.
% Disabled station/class cells carry NaN (no throughput defined); they
% contribute zero to the system total, matching the stage references.
Xup = sum(upT(:)); Xdown = sum(downT(:)); Xenv = sum(TN(isfinite(TN)));
lo = min(Xup, Xdown); hi = max(Xup, Xdown);
tol = 1e-2 * max(1, hi);
fprintf('Total throughput  UP=%.4f  DOWN=%.4f  ENV=%.4f\n', Xup, Xdown, Xenv);
fprintf('Aggregate Q (sum) UP=%.4f  DOWN=%.4f  ENV=%.4f\n', sum(upQ(:)), sum(downQ(:)), sum(QN(:)));
assert(Xenv >= lo - tol && Xenv <= hi + tol, ...
    'ENV-averaged throughput is not bracketed by the single-stage solutions');

% (3) Quantitative coupling: the env-averaged throughput increases with the
% stationary probability of the fast UP stage, P(UP)=b/(a+b) for switch rates
% a (UP->DOWN) and b (DOWN->UP). Both settings remain within the bracket.
Xlow  = env2StageTput(1.0, 0.2, T);   % P(UP)=0.167 (mostly slow DOWN)
Xhigh = env2StageTput(0.2, 1.0, T);   % P(UP)=0.833 (mostly fast UP)
fprintf('Monotonicity   P(UP)=0.167 -> %.4f   P(UP)=0.833 -> %.4f\n', Xlow, Xhigh);
assert(Xhigh > Xlow + 1e-3, 'ENV throughput is not monotone in P(UP)');
assert(Xlow >= lo - tol && Xhigh <= hi + tol, 'ENV throughputs escape the single-stage bracket');

% (4) Three-stage environment (UP/MID/DOWN) stays bracketed by the extreme
% single-stage solutions (exercises the E>2 coupling).
X3 = env3StageTput(T);
fprintf('Three-stage ENV throughput = %.4f (bracket [%.4f, %.4f])\n', X3, lo, hi);
assert(X3 >= lo - tol && X3 <= hi + tol, 'Three-stage ENV-averaged throughput is not bracketed');

fprintf('PASS: ENV-over-LQN meanfield ran; throughput bracketed, monotone in P(UP), 3-stage bracketed.\n');
end

function X = env2StageTput(a, b, T)
% Env-averaged total throughput for the two-stage UP/DOWN environment at switch
% rates a (UP->DOWN) and b (DOWN->UP).
upModel   = buildLQN('UP',   0.8);
downModel = buildLQN('DOWN', 3.0);
env = Environment('R');
env.addStage('UP',   'operational', upModel);
env.addStage('DOWN', 'degraded',    downModel);
env.addTransition('UP',   'DOWN', Exp(a));
env.addTransition('DOWN', 'UP',   Exp(b));
env.init();
X = envTputSum(env, T);
end

function X = env3StageTput(T)
% Env-averaged total throughput for a three-stage UP/MID/DOWN environment.
upModel   = buildLQN('UP',   0.8);
midModel  = buildLQN('MID',  1.6);
downModel = buildLQN('DOWN', 3.0);
env = Environment('R3');
env.addStage('UP',   'operational', upModel);
env.addStage('MID',  'degraded',    midModel);
env.addStage('DOWN', 'failed',      downModel);
env.addTransition('UP',   'MID',  Exp(0.3));
env.addTransition('MID',  'DOWN', Exp(0.3));
env.addTransition('DOWN', 'MID',  Exp(0.6));
env.addTransition('MID',  'UP',   Exp(0.6));
env.init();
X = envTputSum(env, T);
end

function X = envTputSum(env, T)
% Solve an LQN-in-ENV and return the finite aggregate throughput sum.
fldFactory = @(mm) SolverFluid(mm, 'verbose', false);
lnFactory  = @(m)  SolverLN(m, fldFactory, 'timespan', [0, T], 'verbose', false);
opt = SolverENV.defaultOptions;
opt.iter_max = 10; opt.iter_tol = 0.03; opt.verbose = false;
s = SolverENV(env, lnFactory, opt);
[~,~,~,TN] = s.getAvg();
X = sum(TN(isfinite(TN)));
end

function model = buildLQN(name, dbMean)
model = LayeredNetwork(name);
P1 = Processor(model, 'ClientProcessor', 1, SchedStrategy.PS); %#ok<NASGU>
P2 = Processor(model, 'DBProcessor', 1, SchedStrategy.PS); %#ok<NASGU>
T1 = Task(model, 'ClientTask', 5, SchedStrategy.REF).on(P1);
T1.setThinkTime(Exp.fitMean(5.0));
T2 = Task(model, 'DBTask', Inf, SchedStrategy.INF).on(P2);
E1 = Entry(model, 'ClientEntry').on(T1);
E2 = Entry(model, 'DBEntry').on(T2);
A1 = Activity(model, 'ClientActivity', Exp.fitMean(1.0)).on(T1);
A1.boundTo(E1).synchCall(E2, 2.5);
A2 = Activity(model, 'DBActivity', Exp.fitMean(dbMean)).on(T2);
A2.boundTo(E2).repliesTo(E2);
end

function [Q, T] = stageAggregate(model, fldFactory, Tspan)
% Steady transient aggregate (block-diagonal M x K) queue lengths and
% throughputs for one LQN stage, using the same SolverLN layout SolverENV
% consumes.
s = SolverLN(model, fldFactory, 'timespan', [0, Tspan], 'verbose', false);
[Qt,~,Tt] = s.getTranAvg();
M = size(Qt,1); K = size(Qt,2);
Q = zeros(M,K); T = zeros(M,K);
for i=1:M
    for k=1:K
        Q(i,k) = tailValue(Qt{i,k});
        T(i,k) = tailValue(Tt{i,k});
    end
end
end

function v = tailValue(cell_ik)
v = 0;
if isstruct(cell_ik) && isfield(cell_ik,'metric') && ~isempty(cell_ik.metric)
    v = cell_ik.metric(end);
end
end
