function [estVal, fObjFun] = estimator_ubo(self, nodes)
% UBO Utilization-based optimization
% This demand estimator is based on the method proposed in:
%
% Liu, Z., Wynter, L., Xia, C. H. and Zhang, F.
% Parameter inference of queueing models for IT systems using end-to-end measurements
% Performance Evaluation, Elsevier, 2006.
%
% Implements the single-experiment QP formulation (Eq. 9-13, Section 4.2)
% applied as a Bundle over all experiments (Single-QP, Section 4.3).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
% This code is released under the 3-Clause BSD License.

sn = self.model.getStruct;

%% Collect measurements
for n = 1:size(nodes, 2)
    node = nodes{n};
    if isfinite(node.getNumberOfServers())
        U = self.getAggrUtil(node);
        if ~isempty(U)
            avgU{n} = U.data * node.getNumberOfServers();
        end
    end

    for r = 1:sn.nclasses
        avgArvR{n, r} = self.getArvR(node, self.model.classes{r});
        if isempty(avgArvR{n, r})
            error('Arrival rate data for node %d in class %d is missing.', self.model.getNodeIndex(node), r);
        else
            avgArvR{n, r} = avgArvR{n, r}.data;
        end
        avgRespT{n, r} = self.getRespT(node, self.model.classes{r});
        if isempty(avgRespT{n, r})
            error('Response time data for node %d in class %d is missing.', self.model.getNodeIndex(node), r);
        else
            avgRespT{n, r} = avgRespT{n, r}.data;
        end
    end
end

try
    avgA = cell2mat(avgArvR);
    avgA = reshape(avgA, size(avgA,1)/size(avgArvR, 1), size(avgArvR, 1), sn.nclasses);
    avgR = cell2mat(avgRespT);
    avgR = reshape(avgR, size(avgR,1)/size(avgRespT, 1), size(avgRespT, 1), sn.nclasses);
    avgU = cell2mat(avgU);
catch me
    switch me.identifier
        case 'MATLAB:catenate:dimensionMismatch'
            error('Sampled metrics have different number of samples, use interpolate() before starting this estimation algorithm.');
    end
end

[estVal, fObjFun] = ubo_qp(avgU, avgR, avgA);
end

function [demEst, fObjFun] = ubo_qp(cpuUtil, rAvgTimes, avgArvR)
% Solves the UBO demand estimation as a quadratic program following
% Liu et al. 2006, Eq. 9-13 (single experiment) applied as Bundle QP
% (Section 4.3) over all experiments.
%
% Inputs:
%   cpuUtil:    N x M   - measured aggregate utilization per station
%   rAvgTimes:  N x M x R - measured per-station response times
%   avgArvR:    N x M x R - measured arrival rates per station per class
%
% The QP minimizes (Eq. 9):
%   sum_j w_j * delta_j^2 + sum_i epsilon_i^2
% where:
%   delta_j = E_j^e - E_j^m     (end-to-end delay error, Eq. 12)
%   epsilon_i = rho_i^e - rho_i^m  (utilization error, Eq. 11)
%   w_j = lambda_j / sum(lambda)    (class weights)
% subject to s_ji >= 0 (Eq. 13)

%% Clean input data
a = sum(isnan(cpuUtil), 2);
if sum(a) > 0
    cpuUtil = cpuUtil(a == 0, :);
    rAvgTimes = rAvgTimes(a == 0, :, :);
    avgArvR = avgArvR(a == 0, :, :);
end

a = sum(sum(avgArvR, 3), 2) == 0;
if sum(a) > 0
    cpuUtil = cpuUtil(a == 0, :);
    rAvgTimes = rAvgTimes(a == 0, :, :);
    avgArvR = avgArvR(a == 0, :, :);
end

N = size(cpuUtil, 1);   % number of experiments
M = size(cpuUtil, 2);   % number of stations (I in paper)
R = size(rAvgTimes, 3); % number of classes (J in paper)
MR = M * R;

%% Build Bundle QP by summing over experiments
% Variables: s_vec of length MR, where s_vec((r-1)*M + i) = s_{ir}
%
% For each experiment n:
%   beta_i^n = 1/(1 - rho_i^n)                         (p.45)
%   E_r^n = sum_i R_{ir}^n                              (end-to-end delay, Eq. 7 with v=1)
%   w_r^n = lambda_r^n / sum_r lambda_r^n               (Eq. 9)
%   delta_r = sum_i s_{ir} * beta_i^n - E_r^n           (Eq. 12)
%   epsilon_i = sum_r lambda_{ir}^n * s_{ir} - rho_i^n  (Eq. 11)
%
% QP: min (1/2) s^T H s + h^T s + const
%   s.t. s >= 0

H = zeros(MR, MR);
h = zeros(MR, 1);

for n = 1:N
    rho_n = cpuUtil(n, :);          % 1 x M
    beta_n = 1 ./ (1 - rho_n);     % 1 x M

    % Per-station arrival rates and response times for this experiment
    if M == 1 && R == 1
        lambda_n = avgArvR(n, 1, 1);    % scalar
        R_n = rAvgTimes(n, 1, 1);       % scalar
    elseif M == 1
        lambda_n = reshape(avgArvR(n, 1, :), 1, R);    % 1 x R
        R_n = reshape(rAvgTimes(n, 1, :), 1, R);       % 1 x R
    else
        lambda_n = reshape(avgArvR(n, :, :), M, R);    % M x R
        R_n = reshape(rAvgTimes(n, :, :), M, R);       % M x R
    end

    % Class weights w_r = lambda_r / sum(lambda) (Eq. 9)
    % lambda_r = total arrival rate for class r (sum over stations for multi-station)
    lambda_r = sum(lambda_n, 1);        % 1 x R
    w_n = lambda_r / sum(lambda_r);     % 1 x R

    % End-to-end delay: E_r = sum_i R_{ir} (Eq. 7, assuming v=1)
    E_n = sum(R_n, 1);                  % 1 x R

    % Build A_delta: R x MR (Eq. 7)
    % Row r has beta_i at column (r-1)*M + i
    A_delta = zeros(R, MR);
    for r = 1:R
        A_delta(r, (r-1)*M + (1:M)) = beta_n;
    end

    % Build A_epsilon: M x MR (Eq. 8)
    % Row i has lambda_{ir} at column (r-1)*M + i
    A_eps = zeros(M, MR);
    for i = 1:M
        for r = 1:R
            A_eps(i, (r-1)*M + i) = lambda_n(i, r);
        end
    end

    W_n = diag(w_n);

    % Accumulate QP matrices (Bundle, Section 4.3)
    % H^n = 2(A_delta^T W A_delta + A_eps^T A_eps)
    H = H + 2 * (A_delta' * W_n * A_delta + A_eps' * A_eps);
    % h^n = -2(A_delta^T W E^m + A_eps^T rho^m)
    h = h - 2 * (A_delta' * W_n * E_n' + A_eps' * rho_n');
end

%% Solve QP: min (1/2) s^T H s + h^T s, s.t. s >= 0 (Eq. 13-17)
lb = zeros(MR, 1);
opts = optimoptions('quadprog', 'Display', 'off');
[s_vec, fObjFun] = quadprog(H, h, [], [], [], [], lb, [], [], opts);

%% Reshape to M x R
demEst = reshape(s_vec, M, R);
end
