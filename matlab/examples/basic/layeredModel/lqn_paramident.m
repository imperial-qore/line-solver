% LQN_PARAMIDENT Identify hidden LQN parameters from measured performance.
%
% Demonstrates infer_lqn, an Extended Kalman Filter that tracks hidden LQN
% parameters (host demands, think times) from measurable performance data,
% following Zheng, Yang, Woodside, Litoiu, Iszlai, "Tracking Time-Varying
% Parameters in Software Systems with Extended Kalman Filters", CASCON 2005.
%
% Here two parameters are hidden: the reference-task think time (the paper's Z)
% and the P2 host demand of activity AS3 (the paper's service demand S_d). The
% measurable vector is [R(E1), U(P1), U(P2)] (a user response time and two
% processor utilizations). We synthesise a measurement sequence with a step
% change plus noise, then recover the parameter trajectory.

clear; lineStart;
rng(12345);

% ---- base LQN (mirrors lqn_basic.m) ----------------------------------------
model = LayeredNetwork('paramident_LQN');
P1 = Processor(model, 'P1', 2, SchedStrategy.PS);
P2 = Processor(model, 'P2', 3, SchedStrategy.PS);
T1 = Task(model, 'T1', 50, SchedStrategy.REF).on(P1).setThinkTime(Exp(1/2));
T2 = Task(model, 'T2', 50, SchedStrategy.FCFS).on(P1).setThinkTime(Exp(1/3));
T3 = Task(model, 'T3', 25, SchedStrategy.FCFS).on(P2).setThinkTime(Exp(1/4));
E1 = Entry(model, 'E1').on(T1);
E2 = Entry(model, 'E2').on(T2);
E3 = Entry(model, 'E3').on(T3);
A1 = Activity(model, 'AS1', Exp(10)).on(T1).boundTo(E1).synchCall(E2, 1);
A2 = Activity(model, 'AS2', Exp(20)).on(T2).boundTo(E2).synchCall(E3, 5).repliesTo(E2);
A3 = Activity(model, 'AS3', Exp(50)).on(T3).boundTo(E3).repliesTo(E3);

% ---- what to estimate (paramSpec) and what is observed (obsSpec) -----------
paramSpec = struct('type', {'think', 'hostdem'}, 'name', {'T1', 'AS3'});
obsSpec   = struct('metric', {'RespT', 'Util', 'Util'}, 'name', {'E1', 'P1', 'P2'});

% ---- ground-truth parameter trajectory with a step change ------------------
nsteps = 24;
aTrueSeq = zeros(2, nsteps);
aTrueSeq(1, :) = 0.5;                 % think time Z
aTrueSeq(2, :) = 1/50;               % S_d = AS3 mean host demand
aTrueSeq(1, 13:end) = 1.0;            % step: think time doubles at step 13
aTrueSeq(2, 7:18)   = 1/25;          % pulse: S_d doubles for steps 7..18

% ---- synthesise the measurement matrix Z from the true model + noise -------
solveropts = SolverLN.defaultOptions; solveropts.verbose = VerboseLevel.SILENT;
no = numel(obsSpec);
Z = zeros(no, nsteps);
for k = 1:nsteps
    infer_lqn_setparams(model, paramSpec, aTrueSeq(:, k));
    [QN, UN, RN, TN] = SolverLN(model, solveropts).getEnsembleAvg();
    lsn = model.getStruct();
    ztrue = infer_lqn_getobs(lsn.names, ...
        struct('QLen', QN, 'Util', UN, 'RespT', RN, 'Tput', TN), obsSpec);
    Z(:, k) = ztrue .* (1 + 0.02 * randn(no, 1));   % ~2% measurement noise
end

% ---- run the EKF identification --------------------------------------------
opt = struct();
opt.a0 = [0.7; 1/40];       % deliberately wrong initial guess
opt.QFac = 0.1;             % eq 9a drift-noise factor
opt.RFac = 0.2;             % eq 9b measurement-noise factor
opt.gammaT = 1;             % T = T* (measurement interval == system constant)
opt.aTrue = aTrueSeq(:, end);
[estModel, info] = infer_lqn(model, paramSpec, obsSpec, Z, opt);

% ---- report -----------------------------------------------------------------
fprintf('\nStep | Z_true  Z_hat  | Sd_true  Sd_hat  | ||e||\n');
for k = 1:nsteps
    fprintf('%4d | %6.3f  %6.3f | %7.4f  %7.4f | %.3g\n', k, ...
        aTrueSeq(1, k), info.ahat(1, k), aTrueSeq(2, k), info.ahat(2, k), ...
        norm(info.e(:, k)));
end
fprintf('\nFinal estimate:  Z = %.4f (true %.4f),  Sd = %.5f (true %.5f)\n', ...
    info.ahat(1, end), aTrueSeq(1, end), info.ahat(2, end), aTrueSeq(2, end));
fprintf('RMS parameter tracking error Ea = %.4g,  prediction error Er = %.4g\n', ...
    info.Ea, info.Er);
