function [Q,U,R,T,C,X,lG] = solver_mva_lcfsqn(sn, options, lcfsStat, lcfsprStat)
% [Q,U,R,T,C,X,LG] = SOLVER_MVA_LCFSQN(SN, OPTIONS, LCFSSTAT, LCFSPRSTAT)
% Specialized MVA solver for LCFS + LCFS-PR 2-station networks
%
% This function wraps the pfqn_lcfsqn algorithm and maps LINE's data
% structures to/from the algorithm's expected format.
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
%   lG - log of normalizing constant (NaN for this method)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;
nclasses = sn.nclasses;
njobs = sn.njobs;

% Extract service times for each class at each station
% alpha(r) = mean service time at LCFS station for class r
% beta(r) = mean service time at LCFS-PR station for class r
alpha = zeros(1, nclasses);
beta = zeros(1, nclasses);

rates = sn.rates;
for r = 1:nclasses
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

% Call the LCFS MVA algorithm
% Returns: T_lcfs (throughput), Q_lcfs (2xR queue lengths),
%          U_lcfs (2xR utilizations), B_lcfs (2xR back probabilities)
[T_lcfs, Q_lcfs, U_lcfs, ~] = pfqn_lcfsqn_mva(alpha, beta, N);

% Map results back to LINE format
% Initialize output matrices for all stations
Q = zeros(M, nclasses);
U = zeros(M, nclasses);
T = zeros(M, nclasses);
R = zeros(M, nclasses);
X = zeros(1, nclasses);
C = zeros(1, nclasses);

% Map queue lengths
Q(lcfsStat, :) = Q_lcfs(1, :);
Q(lcfsprStat, :) = Q_lcfs(2, :);

% Map utilizations
U(lcfsStat, :) = U_lcfs(1, :);
U(lcfsprStat, :) = U_lcfs(2, :);

% Throughput is the same at all stations in a closed network
for r = 1:nclasses
    if njobs(r) > 0
        X(r) = T_lcfs(r);
        T(lcfsStat, r) = T_lcfs(r);
        T(lcfsprStat, r) = T_lcfs(r);
    end
end

% Compute response times: R = Q / T (using Little's Law)
for k = [lcfsStat, lcfsprStat]
    for r = 1:nclasses
        if T(k, r) > 0
            R(k, r) = Q(k, r) / T(k, r);
        end
    end
end

% Compute cycle times: C = sum of response times at all stations
for r = 1:nclasses
    if njobs(r) > 0
        C(r) = R(lcfsStat, r) + R(lcfsprStat, r);
    end
end

% Log of normalizing constant is not computed by this method
lG = NaN;

end
