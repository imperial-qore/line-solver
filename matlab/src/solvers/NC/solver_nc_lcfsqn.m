function [Q,U,R,T,C,X,lG] = solver_nc_lcfsqn(sn, options, lcfsStat, lcfsprStat)
% [Q,U,R,T,C,X,LG] = SOLVER_NC_LCFSQN(SN, OPTIONS, LCFSSTAT, LCFSPRSTAT)
% Specialized NC solver for LCFS + LCFS-PR 2-station networks
%
% This function wraps the pfqn_lcfsqn_ca algorithm and computes performance
% metrics using the convolution approach.
%
% Parameters:
%   sn         - network structure
%   options    - solver options
%   lcfsStat   - index of the LCFS station
%   lcfsprStat - index of the LCFS-PR station
%
% Returns:
%   Q - queue length matrix (stations x classes)
%   U - utilization matrix (stations x classes)
%   R - response time matrix (stations x classes)
%   T - throughput matrix (stations x classes)
%   C - cycle time vector (1 x classes)
%   X - throughput vector (1 x classes)
%   lG - log of normalizing constant
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;
R = sn.nclasses;
njobs = sn.njobs;

% Extract service times for each class at each station
% alpha(r) = mean service time at LCFS station for class r
% beta(r) = mean service time at LCFS-PR station for class r
alpha = zeros(1, R);
beta = zeros(1, R);

rates = sn.rates;
for r = 1:R
    if njobs(r) > 0
        % Service rate at station for class r
        mu_lcfs = rates(lcfsStat, r);
        mu_lcfspr = rates(lcfsprStat, r);

        if mu_lcfs <= 0 || ~isfinite(mu_lcfs)
            line_error(mfilename, sprintf('Invalid service rate at LCFS station for class %d.', r));
        end
        if mu_lcfspr <= 0 || ~isfinite(mu_lcfspr)
            line_error(mfilename, sprintf('Invalid service rate at LCFS-PR station for class %d.', r));
        end

        alpha(r) = 1 / mu_lcfs;
        beta(r) = 1 / mu_lcfspr;
    end
end

% Get population vector
N = njobs;
K = sum(N);

% Call the LCFS convolution algorithm to get normalizing constant
[G, ~] = pfqn_lcfsqn_ca(alpha, beta, N);

% Compute log normalizing constant
if G > 0
    lG = log(G);
else
    lG = -Inf;
end

% Initialize output matrices
Q_lcfs = zeros(2, R);
U_lcfs = zeros(2, R);
T_lcfs = zeros(2, R);

% Compute throughputs and queue lengths for each class
% Based on the approach in lcfsqn_perf.m
% The permanent formulas assume one job per class (K single-job classes);
% classes with multiplicity N(r) > 1 are expanded into N(r) exchangeable
% single-job copies, with G_exp = G * prod_r N(r)! the distinguishable-jobs
% normalizing constant, and per-copy metrics are scaled back by N(r).
N = N(:)';
alphaE = repelem(alpha, N);
betaE = repelem(beta, N);
Gexp = G * prod(factorial(N));
ecls = 1 + cumsum([0, N(1:end-1)]);   % first expanded copy of each class
for r = 1:R
    if njobs(r) > 0
        e = ecls(r);
        Tcopy = 0; Qcopy = 0;
        for xt = 1:K
            % throughput using permanent calculations
            Tx = make_Tx(alphaE, betaE, xt, K, K, e);
            Tcopy = Tcopy + alphaE(e)^(xt-1) * perm(Tx) / Gexp;
            % queue length at station 1 (LCFS)
            Yx = make_Yx(alphaE, betaE, xt, K, K, e);
            Qcopy = Qcopy + perm(Yx) / Gexp;
        end
        T_lcfs(1:2,r) = N(r) * Tcopy;
        Q_lcfs(1,r) = N(r) * Qcopy;
        % Queue length at station 2 (LCFS-PR) by conservation
        Q_lcfs(2,r) = njobs(r) - Q_lcfs(1,r);
    end
end

% Compute utilizations
for i = 1:2
    for r = 1:R
        if i == 1
            U_lcfs(i,r) = T_lcfs(i,r) * alpha(r);
        else
            U_lcfs(i,r) = T_lcfs(i,r) * beta(r);
        end
    end
end

% Map results back to LINE format
% Initialize output matrices for all stations
Q = zeros(M, R);
U = zeros(M, R);
T = zeros(M, R);
R_resp = zeros(M, R);
X = zeros(1, R);
C = zeros(1, R);

% Map queue lengths
Q(lcfsStat, :) = Q_lcfs(1, :);
Q(lcfsprStat, :) = Q_lcfs(2, :);

% Map utilizations
U(lcfsStat, :) = U_lcfs(1, :);
U(lcfsprStat, :) = U_lcfs(2, :);

% Throughput is the same at all stations in a closed network
for r = 1:R
    if njobs(r) > 0
        X(r) = T_lcfs(1,r);
        T(lcfsStat, r) = T_lcfs(1,r);
        T(lcfsprStat, r) = T_lcfs(2,r);
    end
end

% Compute response times: R = Q / T (using Little's Law)
for k = [lcfsStat, lcfsprStat]
    for r = 1:R
        if T(k, r) > 0
            R_resp(k, r) = Q(k, r) / T(k, r);
        end
    end
end

% Compute cycle times: C = sum of response times at all stations
for r = 1:R
    if njobs(r) > 0
        C(r) = R_resp(lcfsStat, r) + R_resp(lcfsprStat, r);
    end
end

% Return R_resp in variable R (overriding function parameter)
R = R_resp;

end

function Tx = make_Tx(alpha, beta, xt, K, R, r)
% Make Tx matrix for throughput computation
% alpha : vector of length R
% beta  : vector of length R
% xt    : integer
% K     : total population
% R     : number of classes
% r     : class index to exclude

if issym(alpha)
    Tx = sym(zeros(K-1, K-1));
else
    Tx = zeros(K-1, K-1);
end

idx = 0;
for i = 1:R
    if i ~= r
        idx = idx + 1;
        for j = 1:(xt-1)
            Tx(idx, j) = alpha(i)^j;
        end
        for j = 1:(K-xt)
            Tx(idx, xt-1+j) = alpha(i)^(xt+j-1) * beta(i);
        end
    end
end
end

function Y = make_Yx(alpha, beta, xt, K, R, r)
% Make Yx matrix for queue length computation
% alpha : vector of length R
% beta  : vector of length R
% xt    : integer
% K     : total population
% R     : number of classes
% r     : class index

if issym(alpha)
    Y = sym(zeros(K, K));
else
    Y = zeros(K, K);
end

for i = 1:R
    for j = 1:xt
        Y(i, j) = alpha(i)^j;
    end
    for j = 1:(K-xt)
        if i~=r
            Y(i, xt+j) = alpha(i)^(xt+j-1) * beta(i);
        end
    end
end
end
