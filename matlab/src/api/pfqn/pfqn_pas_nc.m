function [G, lG] = pfqn_pas_nc(Z, N, mu, prec, options)
% [G, LG] = PFQN_PAS_NC(Z, N, MU, PREC, OPTIONS)
%
% Normalizing constant of a closed pass-and-swap (P&S) queueing network that
% comprises a single aggregated infinite-server (delay) node and an arbitrary
% number of order-independent (OI) / P&S stations, restricted to the recurrent
% communicating class selected by the placement order PREC.
%
% With a non-empty swap graph the ordered-state chain is reducible (Comte &
% Dorsman, 2021, arXiv:2009.12299): the recurrent communicating classes are the
% placement-order-adhering sets and the product form pi(c) = prod_m Phi_m(c_m)/G_C
% holds per class. This routine returns that per-class constant G_C. It is a
% MICROSTATE routine: it walks the ordered chains position by position, because
% with a placement order the reachable set is a set of ORDERINGS and does not
% collapse onto the count lattice. For the plain OI case (empty PREC, every
% ordering feasible) the count lattice does suffice and PFQN_NCOI computes the
% same G at far lower cost; use this routine only when a placement order is
% present, or as a microstate reference.
%
% Method. Build station M's chain head-first: appending class r at chain
% position k = sum(occ)+1 is admissible iff the placement order allows it (no
% class already placed at that station must come after r), and contributes the
% reciprocal OI prefix rate 1/mu_M(occ+e_r); the chain may be finalized only
% when occ is a placement-order ideal at full multiplicity, whereupon the
% recursion moves to station M-1. When every P&S station has been peeled the
% residual population sits at the delay node with the multinomial weight
% prod_r Z_r^{N_r}/N_r!. With PREC empty this is exactly the balanced-fairness
% convolution of Bonald & Proutiere (2003) evaluated ordering by ordering,
%   Phi(0) = 1,  Phi(n) = (1/mu(n)) * sum_{r: n_r>0} Phi(n - e_r).
%
% COST. The recursion visits one node per feasible ordered prefix, so with an
% empty PREC the node count is sum_{b<=N} C(|b|+M-1,M-1) * |b|!/prod_r b_r!,
% i.e. factorial in the total population sum(N). A placement order prunes the
% orderings (a total order leaves a single one per count split), which is what
% makes the microstate walk affordable in the P&S case.
%
% Parameters:
%   Z    - (1 x R) think-time demand vector of the aggregated delay node.
%          Z(r) = 1/sigma_r for a delay with per-class rate sigma_r.
%   N    - (1 x R) closed population vector, finite.
%   mu   - cell array {1 x M} of function handles, one per P&S station. Each
%          mu{m}(n) returns the total service rate of station m given the
%          per-class occupancy (count) vector n (1 x R). For an OI station the
%          rate depends only on the support of n (which classes are present),
%          i.e. the sum of the capacities of the compatible servers. May be
%          empty to model a pure delay network.
%   prec - placement order. Either a cell {1 x M} of (R x R) precedence
%          matrices, one per station, or a single (R x R) matrix broadcast to
%          every station, with prec{m}(i,j) ~= 0 iff class i must be placed
%          before class j at station m (the closure returned by PAS_PLACEMENT,
%          fed by the global DAG of PAS_SWAP2ORDER). Empty or all-zero means no
%          order: every ordering is feasible and G is the plain OI constant.
%          NOTE the orientation: around a cycle 1->2->...->M->1 the chain of
%          each downstream station is traversed in the opposite direction, so
%          the downstream stations take the TRANSPOSE of the upstream order
%          (prec = {P, P'} for a two-station cycle). Passing the same P to both
%          stations of a cycle silently returns a smaller, wrong G.
%   options - solver options (optional, currently unused; accepted for
%             signature parity with the other pfqn_* routines).
%
% Returns:
%   G  - Normalizing constant G_C of the communicating class.
%   lG - log(G_C).
%
% Example (two P&S stations in a cycle with swap graph SWAP, no delay):
%   H = pas_swap2order(swap, {@(c) mu1, @(c) mu2});
%   P = pas_placement(H);
%   G = pfqn_pas_nc([], N, {@(n) mu1, @(n) mu2}, {P, P'});
%
% See also PFQN_NCOI, PFQN_PAS_IS, PAS_PLACEMENT, PAS_SWAP2ORDER.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5
    options = struct(); %#ok<NASGU>
end
if nargin < 4
    prec = {};
end
if nargin < 3 || isempty(mu)
    mu = {};
end
if ~iscell(mu)
    mu = {mu};
end

R = numel(N);
if isempty(Z)
    Z = zeros(1, R);
end
if numel(Z) ~= R
    line_error(mfilename, 'Z and N must have the same number of classes.');
end
if any(~isfinite(N))
    line_error(mfilename, 'pfqn_pas_nc requires finite (closed) populations.');
end

N = round(N(:)');
Z = Z(:)';
M = numel(mu);

% Normalize the placement order to one logical (R x R) matrix per station.
if ~iscell(prec)
    if isempty(prec)
        prec = {};
    else
        prec = repmat({prec}, 1, M);
    end
end
if isempty(prec)
    prec = repmat({false(R)}, 1, M);
elseif numel(prec) == 1 && M > 1
    prec = repmat(prec, 1, M);
end
if numel(prec) ~= M
    line_error(mfilename, 'prec must supply one precedence matrix per station.');
end
for m = 1:M
    if isempty(prec{m})
        prec{m} = false(R);
    else
        if any(size(prec{m}) ~= [R, R])
            line_error(mfilename, 'each precedence matrix must be R x R.');
        end
        prec{m} = logical(prec{m} ~= 0);
    end
end

G = pas_nc_rec(Z, N, N, mu, prec, M, zeros(1, R));
lG = log(G);
end

function G = pas_nc_rec(Z, N, Norig, mu, prec, m, occ)
% Head-peeling recursion. N is the population still available to stations
% m,m-1,...,1 and to the delay node; Norig is the total population (for the
% order-ideal test); occ is the per-class chain already placed at station m.
R = numel(N);

if m == 0
    % Base case: the residual population N sits at the delay node with
    % unnormalized weight prod_r Z_r^{N_r} / N_r!.
    logf = 0;
    for r = 1:R
        if N(r) > 0
            if Z(r) <= 0
                G = 0; % class r has population but no delay demand: infeasible
                return
            end
            logf = logf + N(r) * log(Z(r)) - gammaln(N(r) + 1);
        end
    end
    G = exp(logf);
    return
end

% Step A: finalize station m's chain here and peel to station m-1, but only if
% the accumulated occupancy is a placement-order ideal of station m (a class
% present here requires all its order-predecessors present at full multiplicity;
% otherwise the split is unreachable within the communicating class).
if pas_nc_isideal(occ, prec{m}, Norig)
    G = pas_nc_rec(Z, N, Norig, mu, prec, m - 1, zeros(1, R));
else
    G = 0;
end

% Step B: append one more class r at the next chain position of station m,
% contributing the reciprocal OI prefix rate.
active = mu{m};
precm = prec{m};
placed = occ > 0;
for r = 1:R
    if N(r) > 0 && ~any(precm(r, placed))
        e_r = zeros(1, R);
        e_r(r) = 1;
        mu_r = active(occ + e_r);
        if mu_r <= 0
            continue
        end
        Nm = N;
        Nm(r) = Nm(r) - 1;
        G = G + (1 / mu_r) * pas_nc_rec(Z, Nm, Norig, mu, prec, m, occ + e_r);
    end
end
end

function tf = pas_nc_isideal(occ, precm, Norig)
% True iff occ is a placement-order ideal at full multiplicity: for every
% i prec j, occ(j) > 0 requires occ(i) = Norig(i). Reduces to support
% downward-closure when Norig == 1, and to "always true" for an empty order.
tf = true;
R = numel(occ);
for i = 1:R
    for j = 1:R
        if precm(i, j) && occ(j) > 0 && occ(i) < Norig(i)
            tf = false;
            return
        end
    end
end
end
