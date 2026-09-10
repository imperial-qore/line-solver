function [p, logp] = fluid_mvn_rectangle(m, C, a, b, npoints)
% [P, LOGP] = FLUID_MVN_RECTANGLE(M, C, A, B, NPOINTS)
%
% Rectangle probability P(A <= Y <= B) for Y ~ Normal(M, C), the multivariate
% normal orthant/cell integral behind @SolverFLD/getProbAggr under the
% moment-closure methods.
%
% The integral has no closed form beyond one dimension, so it is evaluated by
% the separation-of-variables transformation of Genz (1992): the Cholesky
% factor of C turns the rectangle into an iterated integral over the unit
% cube whose integrand is a product of normal-CDF differences, and the first
% coordinate is integrated exactly. The remaining cube is integrated with a
% DETERMINISTIC Richtmyer lattice rule, frac(k*sqrt(p_j)) over the first
% primes, averaged with its antithetic reflection. Determinism is required
% here, not merely convenient: the MATLAB, Java and Python twins must return
% the same number, and a randomized rule would make them agree only in
% distribution.
%
% C may be SINGULAR, which is the common case: a closed population fixes the
% sum of the station coordinates, so the covariance of a station holding a
% whole class is rank deficient. A coordinate whose CONDITIONAL variance
% vanishes is not integrated; it is a hard constraint, contributing 1 when the
% conditional mean falls inside its interval and 0 otherwise.
%
% Parameters:
%   m       - mean vector, d-by-1
%   C       - covariance matrix, d-by-d, symmetric positive SEMI-definite
%   a       - lower corner, d-by-1, -Inf allowed
%   b       - upper corner, d-by-1, +Inf allowed
%   npoints - lattice points per antithetic pair (default 4096)
%
% Returns:
%   p    - the rectangle probability in [0,1]
%   logp - log(p), -Inf when p is zero
%
% See also GETPROBAGGR, SOLVER_FLUID_MOMENTS, FLUID_MIN_CLOSURE.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5 || isempty(npoints)
    npoints = 4096;
end

m = m(:); a = a(:); b = b(:);
d = numel(m);
if d == 0
    p = 1; logp = 0;
    return
end
al = a - m;
bu = b - m;
if any(bu <= al)
    p = 0; logp = -Inf;
    return
end

% scale-relative tolerances: DTOL decides which coordinate carries noise,
% CTOL whether a deterministic coordinate satisfies its constraint
scale = max([1, max(abs(diag(C)))]);
dtol = 1e-12 * scale;
ctol = 1e-6 * sqrt(scale);

L = local_chol_psd(C, dtol);
isInt = diag(L) > 0;
intIdx = find(isInt);
nInt = numel(intIdx);
if nInt == 0
    lastInt = 0;
else
    lastInt = intIdx(end);
end
nw = max(0, nInt - 1);

if nw == 0
    W = zeros(1,0);
else
    if nw > 100
        line_error(mfilename, sprintf(['The lattice rule carries generators for at most 100 integration ' ...
            'dimensions, but this rectangle has %d. Aggregate classes before evaluating the cell.'], nw));
    end
    pr = local_primes(nw);
    k = (1:npoints)';
    Wbase = mod(k * sqrt(pr(:)'), 1);
    W = [Wbase; 1 - Wbase]; % antithetic reflection
end

Np = size(W,1);
y = zeros(Np, d);
f = ones(Np, 1);
kw = 0;
for i = 1:d
    if i > 1
        s = y(:,1:i-1) * L(i,1:i-1)';
    else
        s = zeros(Np,1);
    end
    if isInt(i)
        lo = (al(i) - s) / L(i,i);
        hi = (bu(i) - s) / L(i,i);
        dd = local_phi(lo);
        ee = local_phi(hi);
        f = f .* max(0, ee - dd);
        if i ~= lastInt
            kw = kw + 1;
            u = dd + W(:,kw) .* (ee - dd);
            % the inverse CDF is evaluated strictly inside the unit interval
            u = min(max(u, 1e-15), 1 - 1e-15);
            y(:,i) = local_phiinv(u);
        end
    else
        % zero conditional variance: the coordinate is pinned at S, so the
        % cell is either met or not
        f((s < al(i) - ctol) | (s > bu(i) + ctol)) = 0;
    end
end

p = mean(f);
p = min(max(p, 0), 1);
if p > 0
    logp = log(p);
else
    logp = -Inf;
end
end

function L = local_chol_psd(C, dtol)
% Cholesky factor of a symmetric positive SEMI-definite matrix. A vanishing
% pivot leaves a zero row/column, which the caller reads as a deterministic
% coordinate rather than as a failure.
d = size(C,1);
L = zeros(d);
for i = 1:d
    v = C(i,i) - L(i,1:i-1)*L(i,1:i-1)';
    if v > dtol
        L(i,i) = sqrt(v);
        for j = (i+1):d
            L(j,i) = (C(j,i) - L(j,1:i-1)*L(i,1:i-1)') / L(i,i);
        end
    else
        L(i,i) = 0;
        L((i+1):d,i) = 0;
    end
end
end

function y = local_phi(x)
% standard normal CDF without the Statistics Toolbox
y = 0.5*erfc(-x/sqrt(2));
end

function x = local_phiinv(u)
% standard normal quantile without the Statistics Toolbox
x = sqrt(2)*erfinv(2*u - 1);
end

function pr = local_primes(n)
% first N primes, listed rather than sieved so that the Java and Python twins
% generate the identical lattice
tbl = [2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71 ...
    73 79 83 89 97 101 103 107 109 113 127 131 137 139 149 151 157 163 167 173 ...
    179 181 191 193 197 199 211 223 227 229 233 239 241 251 257 263 269 271 277 281 ...
    283 293 307 311 313 317 331 337 347 349 353 359 367 373 379 383 389 397 401 409 ...
    419 421 431 433 439 443 449 457 461 463 467 479 487 491 499 503 509 521 523 541];
pr = tbl(1:n);
end
