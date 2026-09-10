%{ @file cache_miss_sfifo_rmf.m
 %  @brief Position-resolved mean-field miss rates for strict FIFO(m) caches
 %
 %  @author LINE Development Team
%}

%{
 % @brief Miss rates for strict FIFO(m) replacement via a position-resolved
 % density-dependent population process (DDPP) mean field.
 %
 % @details
 % Strict FIFO(m) is NOT equivalent to RANDOM(m)/FIFO(m). Gast and Van Houdt
 % (SIGMETRICS 2015) prove pi_FIFO(m) = pi_RAND(m) exactly but show strict
 % FIFO(m) differs and give it no mean-field model (only trace simulation).
 % The difference is the within-list age ordering: on a hit in list i < h the
 % demoted tail of list i+1 is reinserted at position 1 of list i (positions
 % 1..j-1 shift back), which the per-item per-list occupancy of RANDOM(m)
 % cannot express. This routine tracks x[k,i,j] = P(item k in position j of
 % list i) with deterministic (age-based) demotion/eviction and returns the
 % plain mean-field fixed point; it reduces to the RANDOM(m)/FIFO(m) result
 % when m_1 = ... = m_{h-1} = 1 (the strict-FIFO == FIFO degeneracy).
 %
 % Strict FIFO(m) dynamics (aggregate IRM stream):
 %   - miss: insert missed item at position 1 of list 1; list 1 shifts back;
 %     tail (position m_1) is evicted.
 %   - hit at position j of list i < h: promote that item to position 1 of
 %     list i+1 (list i+1 shifts back, tail demoted); demoted tail enters
 %     position 1 of list i and positions 1..j-1 of list i shift back.
 %   - hit in the top list h: no change.
 %
 % Reference:
 %   N. Gast and B. Van Houdt, "Transient and Steady-state Regime of a Family
 %   of List-based Cache Replacement Algorithms", ACM SIGMETRICS 2015.
 %
 % @par Syntax:
 % @code
 % [M, MU, MI, pi0] = cache_miss_sfifo_rmf(gamma, m, lambda)
 % @endcode
%}
function [M,MU,MI,pi0,tout,pi0_t,MU_t,xtraj] = cache_miss_sfifo_rmf(gamma, m, lambda, tspan, x0init, accost) %#ok<INUSL>
if nargin < 6, accost = []; end
if nargin < 4, tspan = []; end
if nargin < 5, x0init = []; end
u = size(lambda,1);
n = size(lambda,2);
h = length(m);
m = m(:)';

% aggregate per-item request rates over users
lam_i = zeros(1, n);
for v = 1:u
    row = lambda(v,:,1);
    row(~isfinite(row)) = 0;
    lam_i = lam_i + row;
end
tot = sum(lam_i);
if tot > 0, p = lam_i / tot; else, p = ones(1,n)/n; end

% slot map: (list i, position j) -> flat slot index; item k occupies k*S+s
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

% popularity-ordered warm start
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
    drift_h = @(t, x) cache_pos_drift_graph(x, p, G, m, n, h, slots, sidx, S, 'head');
    x0s = zeros(dim,1);   % cold start so non-admissible items drain
else
    drift_h = @(t, x) sfifo_drift(x, p, m, n, h, slots, sidx, S, dim);
    x0s = x0;
end
odeopt = odeset('AbsTol', 1e-10, 'RelTol', 1e-8);
[~, xvec] = ode15s(drift_h, [0, 20000], x0s, odeopt);
xss = xvec(end, :)';

pi0 = zeros(n,1);
for k = 1:n
    pi0(k) = sfifo_out(xss, k, slots, S);
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
            pi0_t(k,c) = sfifo_out(xtraj(:,c), k, slots, S);
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

function o = sfifo_out(x, k, slots, S)
% Out-of-cache occupancy of item k (1 minus its total in-cache occupancy)
S2 = size(slots,1);
acc = 0;
for s = 1:S2
    acc = acc + x((k-1)*S + s);
end
o = max(0, min(1, 1 - acc));
end

function dX = sfifo_drift(x, p, m, n, h, slots, sidx, S, dim)
% Mean-field drift F(x) for strict FIFO(m); x is dim x 1 over in-cache slots.
x = max(0, min(1, x));
K = @(k,i,j) (k-1)*S + sidx(i,j);

% per-position and per-list hit rates, and the miss rate
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
    M = M + p(k) * sfifo_out(x, k, slots, S);
end

% full-shift rate of each list
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
        % outflow: shift toward j+1 (or leave list at tail); promotion up if i<h
        o = (Sfull(i) + sfifo_gi(i,j,Hpos,m,h)) * xk;
        if i < h
            o = o + p(k) * xk;
        end
        dX(K(k,i,j)) = dX(K(k,i,j)) - o;
        % inflow
        if j >= 2
            dX(K(k,i,j)) = dX(K(k,i,j)) + (Sfull(i) + sfifo_gi(i,j-1,Hpos,m,h)) * x(K(k,i,j-1));
        else
            if i == 1
                dX(K(k,1,1)) = dX(K(k,1,1)) + p(k) * sfifo_out(x, k, slots, S);
            else
                occ = 0;
                for jprev = 1:m(i-1)
                    occ = occ + x(K(k,i-1,jprev));
                end
                dX(K(k,i,1)) = dX(K(k,i,1)) + p(k) * occ;
            end
            if i < h
                dX(K(k,i,1)) = dX(K(k,i,1)) + Hi(i) * x(K(k,i+1,m(i+1)));
            end
        end
    end
end
end

function g = sfifo_gi(i, jp, Hpos, m, h)
% Partial-shift rate of a slot at position jp of list i: aggregate hit rate at
% deeper positions of list i. The top list h never moves on a hit.
if i == h
    g = 0; return;
end
g = 0;
for jj = (jp+1):m(i)
    g = g + Hpos(i,jj);
end
end
