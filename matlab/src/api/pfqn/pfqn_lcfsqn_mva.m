%{
%{
 % @file pfqn_lcfsqn_mva.m
 % @brief Exact MVA for 2-station LCFS queueing networks.
%}
%}

%{
%{
 % @brief Exact MVA for 2-station LCFS queueing networks.
 % @fn pfqn_lcfsqn_mva(alpha, beta, N)
 % @param alpha Service rates at LCFS station (1xR vector).
 % @param beta Service rates at LCFS-PR station (1xR vector).
 % @param N Population vector (default: ones(1,R)).
 % @return T Throughput vector.
 % @return Q Mean queue lengths (2xR matrix).
 % @return U Utilization (2xR matrix).
 % @return B Back probability matrix (2xR).
%}
%}
function [T,Q,U,B] = pfqn_lcfsqn_mva(alpha, beta, N)
% [T,Q,U,B] = PFQN_LCFSQN_MVA(ALPHA, BETA, N)
% Mean Value Analysis for multiclass LCFS queueing networks
%
% This function computes performance metrics for a 2-station closed
% queueing network with:
%   - Station 1: LCFS (Last-Come-First-Served, non-preemptive)
%   - Station 2: LCFS-PR (LCFS with Preemption-Resume)
%
% Parameters:
%   alpha - vector of inverse service rates at station 1 (LCFS)
%           alpha(r) = 1/mu(1,r) for class r
%   beta  - vector of inverse service rates at station 2 (LCFS-PR)
%           beta(r) = 1/mu(2,r) for class r
%   N     - population vector, N(r) = number of jobs of class r
%           (default: ones(1,R) - one job per class)
%
% Returns:
%   T - throughput vector, T(r) = throughput of class r
%   Q - queue length matrix, Q(i,r) = mean queue length at station i, class r
%   U - utilization matrix, U(i,r) = utilization at station i, class r
%   B - back probability matrix, B(i,r) = probability class r job is at
%       back of queue at station i
%
% Note: This implementation uses log-space arithmetic to prevent numerical
% underflow. The results are mathematically exact (up to floating-point
% precision) - no approximations are made.
%
% Reference:
%   G. Casale, "A family of multiclass LCFS queueing networks with
%   order-dependent product-form solutions", QUESTA 2026.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

R = length(alpha);
if nargin < 3
    N = ones(1, R);
end
K = sum(N);

if K == 0
    Q = zeros(2, R);
    T = zeros(1, R);
    B = zeros(2, R);
    U = zeros(2, R);
    return
else
    % Precompute log values for numerical stability
    log_alpha = log(alpha);
    log_beta = log(beta);

    prods = zeros(1, R);
    for r = 1:R
        prods(r) = prod(N(1:r-1) + 1);
    end

    nstates = prod(N + 1);
    QN = cell(1, nstates);
    TN = cell(1, nstates);
    % Store B in log-scaled form to prevent underflow:
    %   B_actual = BN_scaled * exp(log_scale_BN)
    % This representation is EXACT - no approximation.
    BN_scaled = cell(1, nstates);
    log_scale_BN = cell(1, nstates);
    for ci__=1:nstates
        QN{ci__} = zeros(2, R);
        TN{ci__} = zeros(1, R);
        BN_scaled{ci__} = zeros(2, R);
        log_scale_BN{ci__} = -inf(2, R);
    end

    idx = 1; % pre-initialize for codegen
    n = pprod(N);
    while n >= 0
        idx = hashpop(n, N, R, prods);
        QN{idx} = zeros(2, R);
        BN_scaled{idx} = zeros(2, R);
        log_scale_BN{idx} = -inf(2, R);  % -inf represents 0 in log-space
        TN{idx} = zeros(1, R);

        if sum(n) == 1
            % Base case: single job - no scaling needed
            for r = 1:R
                if n(r) > 0
                    denom = alpha(r) + beta(r);
                    QN{idx}(1, r) = alpha(r) / denom;
                    QN{idx}(2, r) = beta(r) / denom;
                    % At n=1, B = Q directly (no prod(alpha.^n) factor)
                    BN_scaled{idx}(1, r) = alpha(r) / denom;
                    BN_scaled{idx}(2, r) = beta(r) / denom;
                    log_scale_BN{idx}(1, r) = 0;  % exp(0) = 1, so actual = scaled
                    log_scale_BN{idx}(2, r) = 0;
                    TN{idx}(1, r) = 1 / denom;
                end
            end
        else
            % Recursive case: multiple jobs
            % log_scale_np = log(prod(alpha.^n)) - can be very negative
            log_scale_np = sum(n .* log_alpha);

            for k = 1:R
                if n(k) > 0
                    idx_k = hashpop(oner(n, k), N, R, prods);

                    % Compute unscaled waiting time contributions
                    % These are O(1) and don't underflow
                    Wnp_unscaled = 1 + QN{idx_k}(1, k);
                    Wpr_unscaled = 1 + QN{idx_k}(2, k);

                    for r = 1:R
                        if r == k; continue; end
                        if n(r) > 0
                            idx_r = hashpop(oner(n, r), N, R, prods);

                            % Compute B ratio in log-space:
                            % ratio = B_{n-e_k}(1,r) / B_{n-e_r}(1,k)
                            % The scaling factors largely cancel in ratios
                            if BN_scaled{idx_k}(1, r) > 0 && BN_scaled{idx_r}(1, k) > 0
                                log_ratio_np = log(BN_scaled{idx_k}(1, r)) + log_scale_BN{idx_k}(1, r) ...
                                             - log(BN_scaled{idx_r}(1, k)) - log_scale_BN{idx_r}(1, k);
                                ratio_np = exp(log_ratio_np);
                                Wnp_unscaled = Wnp_unscaled + (alpha(k) / alpha(r)) * ratio_np * QN{idx_r}(1, k);
                            end

                            if BN_scaled{idx_k}(2, r) > 0 && BN_scaled{idx_r}(2, k) > 0
                                log_ratio_pr = log(BN_scaled{idx_k}(2, r)) + log_scale_BN{idx_k}(2, r) ...
                                             - log(BN_scaled{idx_r}(2, k)) - log_scale_BN{idx_r}(2, k);
                                ratio_pr = exp(log_ratio_pr);
                                Wpr_unscaled = Wpr_unscaled + (alpha(r) / alpha(k)) * ratio_pr * QN{idx_r}(2, k);
                            end
                        end
                    end

                    % log_scale_pr_k = log(alpha(k)^(sum(n)-1) * beta(k))
                    log_scale_pr_k = (sum(n) - 1) * log_alpha(k) + log_beta(k);

                    % Compute Y(k) = n(k) / (Wnp + Wpr) using log-sum-exp
                    % Wnp = exp(log_scale_np) * Wnp_unscaled
                    % Wpr = exp(log_scale_pr_k) * Wpr_unscaled
                    log_Wnp = log_scale_np + log(Wnp_unscaled);
                    log_Wpr = log_scale_pr_k + log(Wpr_unscaled);

                    % log-sum-exp for numerical stability
                    max_log = max(log_Wnp, log_Wpr);
                    log_sum_W = max_log + log(exp(log_Wnp - max_log) + exp(log_Wpr - max_log));

                    % B(1,k) = prod(alpha.^n) * n(k) / (Wnp + Wpr)
                    %        = n(k) * exp(log_scale_np - log_sum_W)
                    BN_scaled{idx}(1, k) = n(k);
                    log_scale_BN{idx}(1, k) = log_scale_np - log_sum_W;

                    % B(2,k) = alpha(k)^(sum(n)-1) * beta(k) * n(k) / (Wnp + Wpr)
                    %        = n(k) * exp(log_scale_pr_k - log_sum_W)
                    BN_scaled{idx}(2, k) = n(k);
                    log_scale_BN{idx}(2, k) = log_scale_pr_k - log_sum_W;
                end
            end

            % Obtain queue lengths: Q(i,k) = B(i,k) + sum_r B(i,r) * Q_{n-e_r}(i,k)
            for k = 1:R
                if n(k) > 0
                    % Accumulate terms using log-sum-exp
                    % Each term: B_actual(i,r) * Q_{n-e_r}(i,k)
                    %          = BN_scaled(i,r) * exp(log_scale_BN(i,r)) * Q_{n-e_r}(i,k)
                    for station = 1:2
                        terms = [];
                        log_scales = [];

                        % First term: B(station,k)
                        if BN_scaled{idx}(station, k) > 0
                            terms(end+1) = BN_scaled{idx}(station, k);
                            log_scales(end+1) = log_scale_BN{idx}(station, k);
                        end

                        % Sum over r: B(station,r) * Q_{n-e_r}(station,k)
                        for r = 1:R
                            if n(r) > 0
                                idx_r = hashpop(oner(n, r), N, R, prods);
                                if BN_scaled{idx}(station, r) > 0 && QN{idx_r}(station, k) > 0
                                    terms(end+1) = BN_scaled{idx}(station, r) * QN{idx_r}(station, k);
                                    log_scales(end+1) = log_scale_BN{idx}(station, r);
                                end
                            end
                        end

                        QN{idx}(station, k) = pfqn_lcfsqn_mva_logsumexp(terms, log_scales);
                    end
                end
            end

            % Obtain throughput
            for k = 1:R
                if n(k) > 0
                    idx_k = hashpop(oner(n, k), N, R, prods);

                    % Compute Unp = sum_r alpha(r) * T_{n-e_k}(r)
                    Unp = 0;
                    for r = 1:R
                        Unp = Unp + alpha(r) * TN{idx_k}(1, r);
                    end

                    % T(k) = sum_r B(1,r) * T_{n-e_r}(k) + (1/alpha(k)) * B(1,k) * (1-Unp)
                    terms = [];
                    log_scales = [];

                    for r = 1:R
                        if n(r) > 0
                            idx_r = hashpop(oner(n, r), N, R, prods);
                            if BN_scaled{idx}(1, r) > 0 && TN{idx_r}(1, k) > 0
                                terms(end+1) = BN_scaled{idx}(1, r) * TN{idx_r}(1, k);
                                log_scales(end+1) = log_scale_BN{idx}(1, r);
                            end
                        end
                    end

                    % Add (1/alpha(k)) * B(1,k) * (1-Unp) term
                    if BN_scaled{idx}(1, k) > 0 && (1 - Unp) > 0
                        terms(end+1) = (1 / alpha(k)) * BN_scaled{idx}(1, k) * (1 - Unp);
                        log_scales(end+1) = log_scale_BN{idx}(1, k);
                    end

                    TN{idx}(1, k) = pfqn_lcfsqn_mva_logsumexp(terms, log_scales);
                end
            end
        end
        n = pprod(n, N);
    end

    Q = QN{idx};
    T = TN{idx};

    % Recover actual B values from scaled representation
    % B_actual = BN_scaled * exp(log_scale_BN)
    B = zeros(2, R);
    for r = 1:R
        if BN_scaled{idx}(1, r) > 0
            B(1, r) = BN_scaled{idx}(1, r) * exp(log_scale_BN{idx}(1, r));
        end
        if BN_scaled{idx}(2, r) > 0
            B(2, r) = BN_scaled{idx}(2, r) * exp(log_scale_BN{idx}(2, r));
        end
    end

    U = zeros(2, R);
    for r = 1:R
        U(1, r) = T(r) * alpha(r);
        U(2, r) = T(r) * beta(r);
    end
end
end

function result = pfqn_lcfsqn_mva_logsumexp(terms, log_scales)
% PFQN_LCFSQN_MVA_LOGSUMEXP Compute sum of scaled terms using log-sum-exp
%
% Computes: result = sum_i (terms(i) * exp(log_scales(i)))
%
% This is mathematically EXACT - we're just computing the sum in a
% numerically stable way by factoring out the maximum exponent.

if isempty(terms)
    result = 0;
    return;
end

% Filter out zero or negative terms (which would give -Inf or complex in log)
valid = terms > 0 & isfinite(log_scales);
if ~any(valid)
    result = 0;
    return;
end
terms = terms(valid);
log_scales = log_scales(valid);

% Compute log of each term: log(term * exp(log_scale)) = log(term) + log_scale
log_values = log(terms) + log_scales;

% Find maximum for numerical stability
max_log = max(log_values);

if ~isfinite(max_log)
    result = 0;
    return;
end

% log-sum-exp: log(sum(exp(x))) = max(x) + log(sum(exp(x - max(x))))
% So: sum(exp(x)) = exp(max(x)) * sum(exp(x - max(x)))
result = exp(max_log) * sum(exp(log_values - max_log));
end
