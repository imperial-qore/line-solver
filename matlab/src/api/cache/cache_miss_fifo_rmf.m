%{ @file cache_miss_fifo_rmf.m
 %  @brief Position-resolved mean-field miss rates for FIFO(m) caches
 %
 %  @author LINE Development Team
%}

%{
 % @brief Miss rates for FIFO(m) replacement via a position-resolved
 % density-dependent population process (DDPP) mean field.
 %
 % @details
 % FIFO(m) and RANDOM(m) share the exact stationary distribution (Gast and Van
 % Houdt, SIGMETRICS 2015, Thm 1: pi_FIFO(m) = pi_RAND(m)), so their steady-
 % state hit ratios coincide and SolverFLD serves FIFO steady state from the
 % cheaper RAND(m) refined mean field (cache_miss_rmf). Their mean-field
 % TRANSIENTS differ: FIFO evicts the deterministic tail (residence of exactly
 % m insertions) whereas RANDOM evicts a uniformly random victim (geometric
 % residence), so H(t) from a cold cache ramps differently even though H(inf)
 % agrees. This routine provides that dedicated FIFO transient.
 %
 % FIFO(m) differs from strict FIFO(m) only in the reinsertion position on a
 % hit: the demoted tail of list i+1 lands at the vacated position j of list i
 % (in place, no within-list shift), whereas strict FIFO reinserts at position 1.
 %
 % Reference:
 %   N. Gast and B. Van Houdt, "Transient and Steady-state Regime of a Family
 %   of List-based Cache Replacement Algorithms", ACM SIGMETRICS 2015.
 %
 % @par Syntax:
 % @code
 % [M, MU, MI, pi0] = cache_miss_fifo_rmf(gamma, m, lambda)
 % @endcode
%}
function [M,MU,MI,pi0,tout,pi0_t,MU_t,xtraj] = cache_miss_fifo_rmf(gamma, m, lambda, tspan, x0init, accost) %#ok<INUSL>
if nargin < 6, accost = []; end
if nargin < 4, tspan = []; end
if nargin < 5, x0init = []; end
u = size(lambda,1);
n = size(lambda,2);
h = length(m);
m = m(:)';

lam_i = zeros(1, n);
for v = 1:u
    row = lambda(v,:,1);
    row(~isfinite(row)) = 0;
    lam_i = lam_i + row;
end
tot = sum(lam_i);
if tot > 0, p = lam_i / tot; else, p = ones(1,n)/n; end

slots = zeros(sum(m), 2);
s = 0;
for i = 1:h
    for j = 1:m(i)
        s = s + 1;
        slots(s,:) = [i j];
    end
end
S = s;
sidx = zeros(h, max(m));
for s2 = 1:S
    sidx(slots(s2,1), slots(s2,2)) = s2;
end
dim = n * S;

[~, order] = sort(p, 'descend');
x0 = zeros(dim,1);
pos = 0;
for s2 = 1:S
    pos = pos + 1;
    if pos <= n
        x0((order(pos)-1)*S + s2) = 1.0;
    end
end

% A non-default access graph (accost) uses the general position-resolved drift
% from a COLD (empty) cache so non-admissible items drain; the linear default
% keeps the pre-filled path.
G = cache_build_item_graphs(accost, lambda, n, h);
if ~isempty(G)
    drift_h = @(t, x) cache_pos_drift_graph(x, p, G, m, n, h, slots, sidx, S, 'pos');
    x0s = zeros(dim,1);   % cold start so non-admissible items drain
else
    drift_h = @(t, x) fifo_drift(x, p, m, n, h, slots, sidx, S, dim);
    x0s = x0;
end
odeopt = odeset('AbsTol', 1e-10, 'RelTol', 1e-8);
[~, xvec] = ode15s(drift_h, [0, 20000], x0s, odeopt);
xss = xvec(end, :)';

pi0 = zeros(n,1);
for k = 1:n
    pi0(k) = fifo_out(xss, k, slots, S);
end
MI = lam_i(:) .* pi0;
MU = zeros(1, u);
for v = 1:u
    row = lambda(v,:,1);
    row(~isfinite(row)) = 0;
    MU(v) = row * pi0;
end
M = sum(MI);

tout = []; pi0_t = []; MU_t = []; xtraj = [];
if ~isempty(tspan)
    if isempty(x0init), x0t = x0s; else, x0t = x0init(:); end
    [tout, xtraj] = ode15s(drift_h, tspan, x0t, odeopt);
    xtraj = xtraj';
    nt = numel(tout);
    pi0_t = zeros(n, nt);
    for k = 1:n
        for c = 1:nt
            pi0_t(k,c) = fifo_out(xtraj(:,c), k, slots, S);
        end
    end
    MU_t = zeros(u, nt);
    for v = 1:u
        row = lambda(v,:,1);
        row(~isfinite(row)) = 0;
        MU_t(v,:) = row * pi0_t;
    end
end
end

function o = fifo_out(x, k, slots, S)
S2 = size(slots,1);
acc = 0;
for s = 1:S2
    acc = acc + x((k-1)*S + s);
end
o = max(0, min(1, 1 - acc));
end

function dX = fifo_drift(x, p, m, n, h, slots, sidx, S, dim)
% Mean-field drift F(x) for FIFO(m); x is dim x 1 over in-cache slots.
x = max(0, min(1, x));
K = @(k,i,j) (k-1)*S + sidx(i,j);

Hpos = zeros(h, max(m));
Hi = zeros(1, h+1);
for s = 1:S
    i = slots(s,1); j = slots(s,2);
    acc = 0;
    for k = 1:n
        acc = acc + p(k) * x(K(k,i,j));
    end
    Hpos(i,j) = acc;
    Hi(i) = Hi(i) + acc;
end
M = 0;
for k = 1:n
    M = M + p(k) * fifo_out(x, k, slots, S);
end

Sfull = zeros(1, h+1);
Sfull(1) = M;
for i = 2:h
    Sfull(i) = Hi(i-1);
end

dX = zeros(dim,1);
for k = 1:n
    for s = 1:S
        i = slots(s,1); j = slots(s,2);
        xk = x(K(k,i,j));
        % outflow: full shift toward j+1 (tail leaves); promotion up if i<h
        o = Sfull(i) * xk;
        if i < h
            o = o + p(k) * xk;
        end
        dX(K(k,i,j)) = dX(K(k,i,j)) - o;
        % inflow: full shift from j-1, or front insertion at j==1
        if j >= 2
            dX(K(k,i,j)) = dX(K(k,i,j)) + Sfull(i) * x(K(k,i,j-1));
        else
            if i == 1
                dX(K(k,1,1)) = dX(K(k,1,1)) + p(k) * fifo_out(x, k, slots, S);
            else
                occ = 0;
                for jprev = 1:m(i-1)
                    occ = occ + x(K(k,i-1,jprev));
                end
                dX(K(k,i,1)) = dX(K(k,i,1)) + p(k) * occ;
            end
        end
        % FIFO demotion: tail of list i+1 lands in place at the same position j
        if i < h
            dX(K(k,i,j)) = dX(K(k,i,j)) + Hpos(i,j) * x(K(k,i+1,m(i+1)));
        end
    end
end
end
