function [P, placeable] = pas_placement(H)
% [P, PLACEABLE] = PAS_PLACEMENT(H)
%
% Placement-order logic of a pass-and-swap (P&S) / order-independent network
% with swap graph H. Isolates the check that decides which class orderings are
% feasible (adhere to the placement partial order) for PFQN_PAS_IS / PFQN_PAS_NC.
%
% An ordering c = (c_1, ..., c_ell) is FEASIBLE iff it is non-decreasing with
% respect to H, i.e. class a never appears before class b whenever H(b,a) ~= 0
% (Comte & Dorsman, 2021, arXiv:2009.12299). Equivalently H(b,a) ~= 0 means b
% must be placed before a. Collecting these constraints and taking the
% transitive closure yields the precedence matrix
%   P(i,j) = 1  iff class i must be placed before class j,
% so an ordering is feasible iff every class is placed only after all of its
% P-predecessors. H may be given as a mere Hasse diagram; the closure makes the
% full order explicit.
%
% Parameters:
%   H - (R x R) swap-graph adjacency (H(b,a) ~= 0 forces b before a). The empty
%       or all-zero graph yields P = 0 (no constraints, pure OI: every ordering
%       feasible).
%
% Returns:
%   P         - (R x R) precedence closure; P(i,j)=1 iff i must precede j.
%   placeable - function handle placeable(x) returning the row vector of class
%               indices that may be placed next given the remaining per-class
%               count vector x (1 x R): those present classes with no remaining
%               predecessor still to be placed. Used to enumerate/sample the
%               feasible orderings and to check placement-order adherence.
%
% See also PFQN_PAS_IS, PFQN_PAS_NC, PAS_SWAP2ORDER.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

R = size(H, 1);
if isempty(H)
    P = [];
    placeable = @(x) find(x > 0);
    return
end
H = (H ~= 0);

% Transitive closure of the "must precede" relation (P(i,j)=1 iff i precedes j,
% i.e. edge i->j in H). Iterate R times so paths of any length are captured.
P = H;
for it = 1:R %#ok<NASGU>
    Pnext = (P | (P * H)) > 0;
    if isequal(Pnext, P)
        break
    end
    P = Pnext;
end
P = double(P);

placeable = @(x) local_placeable(x, P);
end

function idx = local_placeable(x, P)
% Class j is placeable iff it is present (x(j)>0) and no still-present class i
% must precede it: sum_i x(i) P(i,j) == 0.
x = x(:)';
if isempty(P)
    idx = find(x > 0);
else
    idx = find((x * P) == 0 & x > 0);
end
end
