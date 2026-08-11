function [T,rho] = qsys_mm1_dps(lambda, mu, w, tol, maxCutoff)
% [T,RHO] = QSYS_MM1_DPS(LAMBDA, MU, W, TOL, MAXCUTOFF)
%
% Numerically exact M/M/1 Discriminatory Processor Sharing (DPS) queue.
%
% Solves the multiclass DPS continuous-time Markov chain on the per-class
% population vector (n_1..n_K): arrivals lambda(k), class-k completion rate
% mu(k)*n_k*w(k)/sum_j n_j*w(j). The state space is truncated at a total
% population level chosen from the geometric tail bound and doubled until the
% per-class mean counts are stable, so the result is exact to solver precision
% and conserves the M/M/1 total for equal service rates by construction.
%
% Input:
%   lambda - per-class Poisson arrival rates (1,K)
%   mu     - per-class exponential service rates (1,K)
%   w      - per-class DPS weights (1,K), positive
%   tol    - convergence tolerance on the mean counts (default 1e-10)
%   maxCutoff - hard bound on the truncation level (default 2048)
%
% Output:
%   T   - per-class mean response times (1,K) via Little's law
%   rho - total utilization sum_k lambda(k)/mu(k)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 4 || isempty(tol), tol = 1e-10; end
if nargin < 5 || isempty(maxCutoff), maxCutoff = 2048; end

lambda = lambda(:)'; mu = mu(:)'; w = w(:)';
K = length(lambda);
if any(lambda <= 0) || any(mu <= 0) || any(w <= 0)
    line_error(mfilename,'lambda, mu, w must all be positive');
end
rho = sum(lambda ./ mu);
if rho >= 1
    line_error(mfilename,'System is unstable: utilization rho >= 1');
end

N = max(16, ceil(log(tol)/log(rho)));
N = min(N, maxCutoff);
ENprev = solveTrunc(N);
while N < maxCutoff
    N2 = min(2*N, maxCutoff);
    EN = solveTrunc(N2);
    if max(abs(EN - ENprev)) < tol
        ENprev = EN;
        break
    end
    ENprev = EN; N = N2;
    if N2 == maxCutoff, break; end
end
T = ENprev ./ lambda;

    function EN = solveTrunc(Ncut)
        % enumerate states with total population <= Ncut
        grids = cell(1,K);
        [grids{:}] = ndgrid(0:Ncut);
        S = zeros(numel(grids{1}), K);
        for kk=1:K, S(:,kk) = grids{kk}(:); end
        S = S(sum(S,2) <= Ncut, :);
        n = size(S,1);
        keys = S * ((Ncut+1).^(0:K-1))';
        idx = sparse(keys+1, 1, 1:n);
        rows = []; cols = []; vals = [];
        tot = sum(S,2);
        den = S * w';
        for kk=1:K
            % arrivals in class kk (blocked at the truncation level)
            en = tot < Ncut;
            src = find(en);
            dstKeys = keys(src) + (Ncut+1)^(kk-1);
            dst = full(idx(dstKeys+1));
            rows = [rows; src]; cols = [cols; dst]; vals = [vals; repmat(lambda(kk), numel(src), 1)]; %#ok<AGROW>
            % departures in class kk (DPS capacity split)
            en = S(:,kk) > 0;
            src = find(en);
            r = mu(kk) .* S(src,kk) .* w(kk) ./ den(src);
            dstKeys = keys(src) - (Ncut+1)^(kk-1);
            dst = full(idx(dstKeys+1));
            rows = [rows; src]; cols = [cols; dst]; vals = [vals; r]; %#ok<AGROW>
        end
        Q = sparse(rows, cols, vals, n, n);
        Q = Q - diag(sum(Q,2));
        A = Q';
        A(1,:) = 1;
        b = zeros(n,1); b(1) = 1;
        pi = A \ b;
        EN = pi' * S;
    end
end
