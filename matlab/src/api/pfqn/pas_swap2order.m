function H = pas_swap2order(swap, listRate, N0)
% H = PAS_SWAP2ORDER(SWAP, LISTRATE, N0)
%
% Derive the GLOBAL placement-order DAG H of a closed two-station pass-and-swap
% (P&S) tandem 1->2->1 directly from its swap graph, for use with PFQN_PAS_IS.
%
% With a non-empty swap graph the ordered-state chain is reducible (Comte &
% Dorsman, 2021, arXiv:2009.12299): the recurrent communicating class is the set
% of splits (c_{1..k}; c_{ell..k+1}) of the orderings c that are the linear
% extensions of a single placement partial order on the classes (their Prop.).
% PFQN_PAS_IS samples those orderings from H, so it needs exactly this global
% order -- NOT the per-station placement orders consumed by the exact
% convolution PFQN_PAS_NC (which carry the same information one station at a
% time, transposed around the cycle, and are individually over-constrained for
% the single-order IS formulation). Feed H to PAS_PLACEMENT to obtain the
% precedence closure that PFQN_PAS_NC expects.
%
% The placement order is a class-level property, independent of the per-class
% multiplicity, so it is extracted from the minimal single-job-per-class
% instance (N0 = ones(1,R)): enumerate the reachable communicating class from
% the all-in-queue-1 initial state; each reachable state (l1;l2) exposes the
% full ordering c = [l1, reverse(l2)] in D; then set H(i,j)=1 iff class i
% precedes class j in EVERY c in D (forced precedence). Incomparable pairs are
% left 0 (antichain). The result is valid for any population N.
%
% Parameters:
%   swap     - (R x R) swap graph, or cell {G1,G2} of the two per-queue graphs;
%              G(a,b)~=0 means class a chases class b. Empty/all-zero => plain
%              OI: H = 0 (every ordering feasible).
%   listRate - cell {1 x 2} of OI microstate rate functions; listRate{m}(c) is
%              the total service rate of queue m on the ordered prefix c. Used
%              to prune zero-rate (non-head) completions.
%   N0       - (1 x R) minimal probing population; defaults to ones(1,R).
%
% Returns:
%   H - (R x R) global placement-order DAG; H(i,j)=1 iff i must precede j.
%
% See also PFQN_PAS_IS, PAS_PLACEMENT, PFQN_PAS_NC.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = 2;
if ~iscell(swap), swap = repmat({swap}, 1, M); end
if nargin < 3 || isempty(N0)
    R = size(swap{1}, 1);
    N0 = ones(1, R);
end
R = numel(N0);

if (isempty(swap{1}) || ~any(swap{1}(:) ~= 0)) && ...
        (isempty(swap{2}) || ~any(swap{2}(:) ~= 0))
    H = zeros(R, R);            % pure OI: no placement constraint
    return
end

% minimal single-job-per-class initial placement, all jobs at queue 1
initState = [];
for r = 1:R, initState = [initState, repmat(r, 1, N0(r))]; end %#ok<AGROW>
init = {initState, []};

% breadth-first enumeration of the reachable (communicating) class
key = containers.Map('KeyType', 'char', 'ValueType', 'logical');
classStates = {}; frontier = {init}; key(pas_enc(init)) = true;
while ~isempty(frontier)
    st = frontier{end}; frontier(end) = [];
    classStates{end+1} = st; %#ok<AGROW>
    for m = 1:M
        c = st{m}; next = mod(m, M) + 1;
        for p = 1:numel(c)
            % Marginal OI completion rate of the p-th customer. The empty-prefix
            % rate is 0 by definition (no jobs -> no service); some svcRateFun
            % handles return a nonzero constant on [], so guard p==1 explicitly.
            if p == 1
                prevRate = 0;
            else
                prevRate = listRate{m}(c(1:p-1));
            end
            rate = listRate{m}(c(1:p)) - prevRate;
            if rate <= 1e-12, continue; end
            [cnew, dep] = pas_swap_local(c, p, swap{m});
            stn = st; stn{m} = cnew; stn{next} = [st{next}, dep];
            k = pas_enc(stn);
            if ~isKey(key, k), key(k) = true; frontier{end+1} = stn; end %#ok<AGROW>
        end
    end
end

% full orderings D: every reachable state (l1;l2) exposes c = [l1, reverse(l2)]
Dmap = containers.Map('KeyType', 'char', 'ValueType', 'logical');
D = {};
for s = 1:numel(classStates)
    c = [classStates{s}{1}, fliplr(classStates{s}{2})];
    kc = sprintf('%d,', c);
    if ~isKey(Dmap, kc), Dmap(kc) = true; D{end+1} = c; end %#ok<AGROW>
end

% forced precedence: i precedes j iff in every c in D, every copy of i is older
% than every copy of j (position-wise), whenever both are present.
H = false(R);
for a = 1:R
    for b = 1:R
        if a == b, continue; end
        both = false; forced = true;
        for s = 1:numel(D)
            c = D{s};
            pa = find(c == a); pb = find(c == b);
            if isempty(pa) || isempty(pb), continue; end
            both = true;
            if ~(max(pa) < min(pb)), forced = false; break; end
        end
        H(a, b) = both && forced;
    end
end
H = double(H);
end

% ------------------------------------------------------------------------
function s = pas_enc(st)
parts = cell(1, numel(st));
for m = 1:numel(st), parts{m} = sprintf('%d,', st{m}); end
s = strjoin(parts, '|');
end

function [cnew, dep] = pas_swap_local(c, p, G)
n = numel(c); chain = p; moving = c(p); cur = p;
while true
    q = -1;
    for j = cur+1:n, if ~isempty(G) && G(moving, c(j)) ~= 0, q = j; break; end; end
    if q == -1, break; end
    chain(end+1) = q; moving = c(q); cur = q; %#ok<AGROW>
end
dep = c(chain(end)); tmp = c;
for i = 1:numel(chain)-1, tmp(chain(i+1)) = c(chain(i)); end
tmp(chain(1)) = []; cnew = tmp;
end
