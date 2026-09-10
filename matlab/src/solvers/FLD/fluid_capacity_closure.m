function [h, dh, d2h] = fluid_capacity_closure(n, c, s2, lldrow, isInf)
% [H, DH, D2H] = FLUID_CAPACITY_CLOSURE(N, C, S2, LLDROW, ISINF)
%
% Moment closure of the station capacity term psi(X) and of its derivative.
%
% Every scheduling branch of ODE_RATES_CLOSING_FACTORS scales the coordinates
% of a station by psi(n_i)/n_i, where psi is how much work the station clears
% at population n_i:
%
%   psi(n) = min(n,c) * alpha(n)   at a queueing station
%   psi(n) = n        * alpha(n)   at an infinite server
%
% and alpha is the limited load-dependent scaling (1 when the station has
% none). This function returns E[psi(X)] for X ~ Normal(n, s2) together with
% d/dn E[psi(X)], which is what the drift and the Jacobian need. S2 = 0 gives
% psi(n) itself, so the first-order closure is the same code path.
%
% Without load dependence the expectation is the closed form of
% FLUID_MIN_CLOSURE. With a tabulated alpha, FLUID_LLD_SCALING makes alpha
% piecewise linear on the integer lattice, so psi is piecewise QUADRATIC with
% breakpoints at the integers and at c, and the expectation is integrated
% segment by segment against the normal density using the truncated moments
% M0, M1, M2. The derivative is E[psi'(X)] by differentiation under the
% integral sign, valid because psi is Lipschitz.
%
% The integration must be exact, NOT a fixed quadrature rule. Gauss-Hermite
% with fixed nodes applied to a piecewise-linear integrand does not smooth its
% kinks, it relocates them: the resulting estimate of E[psi] is itself
% piecewise linear in n, so its second derivative is zero almost everywhere
% and FLUID_REFINE_MEANFIELD silently returns a null correction. That is what
% the segment-wise closed form below avoids.
%
% psi is extended by zero below n = 0, the only physically admissible
% continuation, since a station holding no jobs clears no work.
%
% Parameters:
%   n      - mean population at the station
%   c      - number of servers (ignored when ISINF)
%   s2     - population variance, 0 for the first-order closure
%   lldrow - (1 x lldlimit) load-dependent scaling, empty when absent
%   isInf  - true at an infinite-server station
%
% Returns:
%   h   - E[psi(X)]
%   dh  - d/dn E[psi(X)]
%   d2h - d2/dn2 E[psi(X)], which FLUID_DRIFT_JACOBIAN needs for the joint
%         product closure of the share and the capacity. psi is only piecewise
%         smooth, so this is E[psi''(X)] in the distributional sense: the
%         segment-wise quadratic term PLUS an atom at every breakpoint where
%         psi' jumps, weighted by the normal density there. Dropping the atoms
%         would report d2h = 0 for the pure min(), whose curvature is carried
%         by nothing else.
%
% See also FLUID_MIN_CLOSURE, FLUID_LLD_SCALING, ODE_RATES_CLOSING_FACTORS.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 4
    lldrow = [];
end
if nargin < 5 || isempty(isInf)
    isInf = false;
end

if isempty(lldrow)
    if isInf
        h = n;
        dh = 1;
        d2h = 0;
    else
        [h, dh, d2h] = fluid_min_closure(n, c, s2);
        % A population cannot be negative, but the normal marginal puts mass
        % below zero, and there min(X,c) = X < 0. Once n is small against the
        % standard deviation that drives E[min(X,c)] itself negative, giving
        % negative service rates and mass destroyed by the non-negativity
        % clamp of the ODE solver: on a closed class-switching model the total
        % population oscillated 4.0 -> 0.16 -> 4.0 between outer iterations.
        % A floor at zero is the minimal repair. Subtracting the whole
        % negative tail instead, E[min(X,c);X>0], is more self-consistent but
        % departs from the published closure everywhere rather than only where
        % it is unusable, and measurably lost accuracy in the mid-load range
        % (e.g. U=0.88 on the think-rate sweep: 0.55% error became 1.19%).
        if h < 0
            h = 0;
            dh = 0;
            d2h = 0;
        end
    end
    return
end

L = numel(lldrow);
if ~(s2 > 0)
    [h, dh, d2h] = local_psi(n, c, lldrow, L, isInf);
    return
end

s = sqrt(s2);
% breakpoints of psi: the lattice of the load-dependence table, the origin,
% and the saturation point c when it falls beyond the table
bps = 0:L;
if ~isInf && isfinite(c) && c > L
    bps = [bps, c];
end
bps = unique(bps);

h = 0; dh = 0; d2h = 0;
% psi is extended by zero below the first breakpoint, so psi' jumps there too
% and that atom belongs in psi'' exactly like the interior ones
Bprev = 0; Cprev = 0;
for k = 1:numel(bps)
    p = bps(k);
    if k < numel(bps)
        q = bps(k+1);
    else
        q = Inf;
    end
    [A, B, C] = local_segment(p, q, c, lldrow, L, isInf);
    [M0, M1, M2] = local_moments(p, q, n, s);
    h = h + A*M0 + B*M1 + C*M2;
    dh = dh + B*M0 + 2*C*M1;
    d2h = d2h + 2*C*M0;
    % psi' jumps across this breakpoint, so psi'' carries an atom there; the
    % segment sum above sees only the quadratic part and would miss it
    jump = (B + 2*C*p) - (Bprev + 2*Cprev*p);
    d2h = d2h + jump * exp(-0.5*((p-n)/s)^2)/(s*sqrt(2*pi));
    Bprev = B; Cprev = C;
end
end

function [A, B, C] = local_segment(p, q, c, lldrow, L, isInf)
% psi(u) = A + B*u + C*u^2 on [p,q], from base(u)*alpha(u)
% alpha is constant on [0,1] and on [L,inf), linear on each unit interval
if p >= L
    a0 = lldrow(L); a1 = 0;
elseif p < 1
    a0 = lldrow(1); a1 = 0;
else
    k = floor(p);
    a1 = lldrow(k+1) - lldrow(k);
    a0 = lldrow(k) - a1*k;
end
% base(u) is u below the saturation point and c above it; the breakpoints
% guarantee the segment lies entirely on one side
if isInf || ~isfinite(c) || q <= c
    A = 0; B = a0; C = a1;      % base = u
else
    A = c*a0; B = c*a1; C = 0;  % base = c
end
end

function [M0, M1, M2] = local_moments(p, q, n, s)
% truncated moments E[X^j * 1{p < X < q}] for X ~ Normal(n, s^2)
zp = (p-n)/s;
zq = (q-n)/s;
Pp = 0.5*erfc(-zp/sqrt(2));
Pq = 0.5*erfc(-zq/sqrt(2));
pp = exp(-0.5*zp^2)/sqrt(2*pi);
if isinf(q)
    pq = 0;
else
    pq = exp(-0.5*zq^2)/sqrt(2*pi);
end
M0 = Pq - Pp;
M1 = n*M0 + s*(pp - pq);
if isinf(q)
    M2 = (n^2+s^2)*M0 + s*((p+n)*pp);
else
    M2 = (n^2+s^2)*M0 + s*((p+n)*pp - (q+n)*pq);
end
end

function [p, dp, d2p] = local_psi(u, c, lldrow, L, isInf) %#ok<INUSD>
% psi and its derivatives at a continuous population, extended by zero below 0.
% alpha is piecewise linear and base piecewise linear, so psi is piecewise
% quadratic and psi'' = 2*base'*alpha' away from the breakpoints; the atoms at
% the breakpoints are not representable without a marginal to smooth them, and
% the first-order closure this branch serves does not smooth them either.
[a, da] = fluid_lld_scaling(lldrow, u);
if isInf
    base = u;
    dbase = 1;
else
    base = min(u, c);
    dbase = double(u < c);
end
p = base.*a;
dp = dbase.*a + base.*da;
d2p = 2*dbase.*da;
if u <= 0
    p = 0; dp = 0; d2p = 0;
end
end
