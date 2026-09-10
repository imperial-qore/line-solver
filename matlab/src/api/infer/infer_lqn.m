function [model, info] = infer_lqn(model, paramSpec, obsSpec, Z, options)
% INFER_LQN Identify hidden LQN parameters from measured performance data.
%
%   [MODEL, INFO] = INFER_LQN(MODEL, PARAMSPEC, OBSSPEC, Z, OPTIONS) estimates
%   the LQN parameters named in PARAMSPEC (activity host demands and/or task
%   think times) of the LayeredNetwork MODEL from the sequence of performance
%   measurements Z, using an Extended Kalman Filter over the observation model
%   defined by OBSSPEC. It implements Zheng, Yang, Woodside, Litoiu, Iszlai,
%   "Tracking Time-Varying Parameters in Software Systems with Extended Kalman
%   Filters", CASCON 2005. A single measurement column with OPTIONS.QFac = 0
%   reduces to one-shot least-squares calibration.
%
%   Inputs:
%     PARAMSPEC : struct array of parameters to estimate; see INFER_LQN_SETPARAMS.
%     OBSSPEC   : struct array of observed metrics; see INFER_LQN_GETOBS.
%     Z         : (no x nsteps) measurements, no == numel(OBSSPEC).
%     OPTIONS   : struct with optional fields (defaults in parentheses):
%         solver        (@SolverLN)  observation-model solver constructor
%         solveropts    ([])         options struct passed to the solver
%         QFac          (0.1)        drift-noise factor, Q_ii=(QFac*a0_i*cvA)^2
%         RFac          (0.2)        meas.-noise factor, R_ii=((RFac*zbar_i)/1.96)^2/gammaT
%         cvA           (1)          parameter drift coefficient of variation
%         gammaT        ([])         T/Tstar ratio for R (else from T and Tstar, else 1)
%         T, Tstar      ([])         measurement interval and system constant
%         Q, R          ([])         explicit covariances (override the above)
%         P0            ([])         initial covariance (default diag((0.5*a0)^2))
%         a0            ([])         initial estimate (default: current model values)
%         aTrue         ([])         ground truth, enables the Ea RMS metric
%         fdStep (1e-3), fdFloor (1e-6), clampPositive (true), verbose (false)
%
%   Outputs:
%     MODEL : the LayeredNetwork with the final parameter estimate applied.
%     INFO  : EKF result struct (see INFER_LQN_EKF) augmented with ahat, a0,
%             Q, R, P0.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5, options = struct(); end

no = numel(obsSpec);
if size(Z, 1) ~= no
    line_error(mfilename, 'Row count of Z must equal numel(obsSpec).');
end

% observation-model solver and its (silent) options
solverCtor = infer_lqn_optget(options, 'solver', @SolverLN);
solveropts = infer_lqn_optget(options, 'solveropts', []);
if isempty(solveropts)
    try
        solveropts = feval([func2str(solverCtor) '.defaultOptions']);
        solveropts.verbose = VerboseLevel.SILENT;
    catch
        solveropts = [];
    end
end

% initial parameter estimate (default: current model values)
a0 = infer_lqn_optget(options, 'a0', []);
if isempty(a0)
    a0 = infer_lqn_getparams(model, paramSpec);
end
a0 = a0(:);

% drift covariance Q (eq 9a) and measurement covariance R (eq 9b)
QFac = infer_lqn_optget(options, 'QFac', 0.1);
RFac = infer_lqn_optget(options, 'RFac', 0.2);
cvA  = infer_lqn_optget(options, 'cvA', 1);
gammaT = infer_lqn_optget(options, 'gammaT', []);
if isempty(gammaT)
    T = infer_lqn_optget(options, 'T', []);
    Tstar = infer_lqn_optget(options, 'Tstar', []);
    if ~isempty(T) && ~isempty(Tstar)
        gammaT = T / Tstar;
    else
        gammaT = 1;
    end
end

Q = infer_lqn_optget(options, 'Q', []);
if isempty(Q)
    qd = (QFac * abs(a0) * cvA).^2;              % eq 9a, mean(a_i) ~ a0_i
    Q = diag(max(qd, eps));
end
R = infer_lqn_optget(options, 'R', []);
if isempty(R)
    zbar = mean(Z, 2);                            % mean of z_i across steps
    rd = ((RFac * abs(zbar)) / 1.96).^2 / gammaT; % eq 9b
    R = diag(max(rd, eps));
end
P0 = infer_lqn_optget(options, 'P0', []);
if isempty(P0)
    P0 = diag(max((0.5 * abs(a0)).^2, eps));
end

% observation model h(a): inject params, solve the LQN, read the metrics
hfun = @(a) evalObs(model, paramSpec, obsSpec, a, solverCtor, solveropts);

ekfopts = struct( ...
    'fdStep',        infer_lqn_optget(options, 'fdStep', 1e-3), ...
    'fdFloor',       infer_lqn_optget(options, 'fdFloor', 1e-6), ...
    'clampPositive', infer_lqn_optget(options, 'clampPositive', true), ...
    'aTrue',         infer_lqn_optget(options, 'aTrue', []), ...
    'verbose',       infer_lqn_optget(options, 'verbose', false));

[ahat, info] = infer_lqn_ekf(hfun, a0, P0, Z, Q, R, ekfopts);
info.ahat = ahat;
info.a0 = a0;
info.Q = Q;
info.R = R;
info.P0 = P0;

% apply the final estimate to the returned model
model = infer_lqn_setparams(model, paramSpec, ahat(:, end));
end

function z = evalObs(model, paramSpec, obsSpec, a, solverCtor, solveropts)
% Inject parameter vector a, solve the LQN, extract the observation vector.
model = infer_lqn_setparams(model, paramSpec, a);
if isempty(solveropts)
    solver = solverCtor(model);
else
    solver = solverCtor(model, solveropts);
end
[QN, UN, RN, TN] = solver.getEnsembleAvg();
lsn = model.getStruct();
metrics = struct('QLen', QN, 'Util', UN, 'RespT', RN, 'Tput', TN);
z = infer_lqn_getobs(lsn.names, metrics, obsSpec);
end

function a0 = infer_lqn_getparams(model, paramSpec)
% Read the current values of the parameters named in paramSpec.
np = numel(paramSpec);
a0 = zeros(np, 1);
for i = 1:np
    switch lower(paramSpec(i).type)
        case 'hostdem'
            act = infer_lqn_findbyname(model.activities, paramSpec(i).name);
            if isempty(act)
                line_error(mfilename, sprintf('Activity ''%s'' not found.', paramSpec(i).name));
            end
            a0(i) = act.hostDemandMean;
        case 'think'
            tsk = infer_lqn_findbyname(model.tasks, paramSpec(i).name);
            if isempty(tsk)
                line_error(mfilename, sprintf('Task ''%s'' not found.', paramSpec(i).name));
            end
            a0(i) = tsk.thinkTimeMean;
        otherwise
            line_error(mfilename, sprintf('Unknown parameter type ''%s''.', paramSpec(i).type));
    end
end
end
