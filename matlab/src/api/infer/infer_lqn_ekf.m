function [ahat, info] = infer_lqn_ekf(hfun, a0, P0, Z, Q, R, options)
% INFER_LQN_EKF Extended Kalman Filter for LQN parameter identification.
%
%   [AHAT, INFO] = INFER_LQN_EKF(HFUN, A0, P0, Z, Q, R, OPTIONS) tracks a
%   hidden parameter vector across the measurement sequence Z using an Extended
%   Kalman Filter, following Zheng, Yang, Woodside, Litoiu, Iszlai, "Tracking
%   Time-Varying Parameters in Software Systems with Extended Kalman Filters",
%   CASCON 2005 (equations 1-9). The parameter is modelled as a zero-mean
%   random walk a_k = a_{k-1} + w and the measurement as z_k = h(a_k) + v,
%   where h is the (nonlinear) LQN performance model supplied as HFUN.
%
%   Inputs:
%     HFUN    : handle mapping a parameter vector to a predicted observation
%               vector z = h(a) (evaluates the LQN model).
%     A0      : (np x 1) initial parameter estimate.
%     P0      : (np x np) initial estimation-error covariance.
%     Z       : (no x nsteps) measurement matrix, one column per step.
%     Q       : (np x np) parameter-drift (process-noise) covariance.
%     R       : (no x no) measurement-error covariance.
%     OPTIONS : struct with optional fields:
%                 fdStep (1e-3), fdFloor (1e-6)  - finite-difference steps
%                 clampPositive (true)           - clamp estimates to > fdFloor
%                 aTrue ([])                      - ground truth for Ea metric
%                 verbose (false)
%
%   Outputs:
%     AHAT : (np x nsteps) parameter estimate trajectory.
%     INFO : struct with fields P (final covariance), Phist (per-step
%            covariances), e (no x nsteps prediction errors), zpred (predicted
%            measurements), Er (prediction RMS), Ea (parameter tracking RMS vs
%            aTrue when supplied, else []).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 7, options = struct(); end
fdStep        = infer_lqn_optget(options, 'fdStep', 1e-3);
fdFloor       = infer_lqn_optget(options, 'fdFloor', 1e-6);
clampPositive = infer_lqn_optget(options, 'clampPositive', true);
aTrue         = infer_lqn_optget(options, 'aTrue', []);
verbose       = infer_lqn_optget(options, 'verbose', false);

a0 = a0(:);
np = numel(a0);
no = size(Z, 1);
nsteps = size(Z, 2);

ahat = zeros(np, nsteps);
info.e = zeros(no, nsteps);
info.zpred = zeros(no, nsteps);
info.Phist = cell(1, nsteps);

a = a0;
P = P0;
I_np = eye(np);
for k = 1:nsteps
    % (1) predict: zero-mean drift keeps a_pred = a; project covariance (eq 5)
    aPred = a;
    Ppred = P + Q;

    % (2,4) predicted measurement and sensitivity matrix H = dh/da at a_pred
    [H, zpred] = infer_lqn_jacobian(hfun, aPred, fdStep, fdFloor);

    % (3) prediction error
    zk = Z(:, k);
    e = zk - zpred;

    % (6) Kalman gain (suboptimal because h is nonlinear)
    S = H * Ppred * H' + R;
    K = (Ppred * H') / S;

    % (4-update) improved parameter estimate
    a = aPred + K * e;
    if clampPositive
        a = max(a, fdFloor);
    end

    % (7) covariance update; symmetrize for numerical stability
    P = (I_np - K * H) * Ppred;
    P = 0.5 * (P + P');

    ahat(:, k) = a;
    info.e(:, k) = e;
    info.zpred(:, k) = zpred;
    info.Phist{k} = P;

    if verbose
        line_printf('[infer_lqn_ekf] step %d/%d  ||e||=%.4g\n', k, nsteps, norm(e));
    end
end
info.P = P;

% RMS tracking (Ea) and prediction (Er) errors (paper Section 4)
info.Er = sqrt(mean(info.e(:).^2));
if ~isempty(aTrue)
    aTrue = aTrue(:);
    D = ahat - repmat(aTrue, 1, nsteps);
    info.Ea = sqrt(mean(D(:).^2));
else
    info.Ea = [];
end
end
