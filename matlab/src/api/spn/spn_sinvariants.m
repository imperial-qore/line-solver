function inv = spn_sinvariants(sn, init)
% INV = SPN_SINVARIANTS(SN)
% INV = SPN_SINVARIANTS(SN, INIT)
% Minimal-support S-invariants (P-invariants) of a stochastic Petri net, and
% the load vector V = S m0.
%
% An S-invariant is a non-negative left null vector of the incidence matrix,
% U' C = 0, so U' m is conserved by every firing. The minimal-support ones form
% a basis of all of them and are what the convolution algorithm SPN_CONV
% decomposes the reachability set along; SPN_MDD uses a single positive
% invariant for a much weaker purpose, to bound each place a priori.
%
% FARKAS' ALGORITHM, on [C | I]: for each transition column in turn, keep the
% rows that already annihilate it and add, for every pair of rows of opposite
% sign in it, the positive combination that cancels it; then drop every row
% whose support strictly contains another's, which is what leaves the minimal
% supports. Rows are kept in integer arithmetic and divided by their gcd, so a
% multiplicity is never lost to rounding and two invariants that differ only by
% a positive scale are the same row.
%
% ARC MULTIPLICITIES MUST BE INTEGRAL. A fractional arc has no Petri-net
% meaning and would make the gcd normalisation and the ILP-free convolution
% both wrong, so it is refused rather than rounded.
%
% Levels are the (place, class) pairs of SPN_MDD, place-major, so the
% invariants come out in the coordinates the decision diagram and SPN_CONV both
% use.
%
% -- Input
% SN   : a NetworkStruct holding Places and Transitions
% INIT : initial tokens per place level, place-major; omitted or [] takes them
%        from the reference station of each closed class, as SPN_MDD does
% -- Output
% INV : struct with fields places (node indices), S (one row per invariant, one
%       column per place level), V = S m0 and m0
%
% -- Reference
% S. Balsamo, A. Marin, I. Stojic, "Computation of the normalising constant for
% product-form models of distributed systems with synchronisation", Future
% Generation Computer Systems 111 (2020) 475-490, Sec. 3.1.
%
% See also SPN_CONV, SPN_MDD.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2, init = []; end
places = find(sn.nodetype == NodeType.Place);
transitions = find(sn.nodetype == NodeType.Transition);
if isempty(places) || isempty(transitions)
    line_error(mfilename, 'the model holds no Place or no Transition node');
end
places = places(:)';
R = sn.nclasses;
P = numel(places);
n = P * R;

% ---- incidence matrix C(level, mode) = post - pre, one column per (transition, mode)
C = [];
for ind = transitions(:)'
    param = sn.nodeparam{ind};
    nmodes = param.nmodes;
    for m = 1:nmodes
        pre = zeros(sn.nnodes, R);
        post = zeros(sn.nnodes, R);
        if numel(param.enabling) >= m && ~isempty(param.enabling{m})
            pre = reshape(param.enabling{m}, sn.nnodes, R);
        end
        if numel(param.firing) >= m && ~isempty(param.firing{m})
            post = reshape(param.firing{m}, sn.nnodes, R);
        end
        col = zeros(n, 1);
        for pp = 1:P
            for k = 1:R
                col((pp - 1) * R + k) = ...
                    i_asinteger(max(0, post(places(pp), k)), 'a firing arc') - ...
                    i_asinteger(max(0, pre(places(pp), k)), 'an enabling arc');
            end
        end
        C = [C, col]; %#ok<AGROW>
    end
end
ncols = size(C, 2);

% ---- Farkas on [C | I]: row p starts as (C(p,:), e_p)
rows = [C, eye(n)];
for c = 1:ncols
    nxt = rows(rows(:, c) == 0, :);
    pos = find(rows(:, c) > 0);
    neg = find(rows(:, c) < 0);
    for a = pos(:)'
        for b = neg(:)'
            pa = rows(a, c); nb = -rows(b, c);
            d = gcd(pa, nb);
            combo = (nb / d) * rows(a, :) + (pa / d) * rows(b, :);
            gall = 0;
            for x = combo, gall = gcd(gall, abs(x)); end
            if gall > 1, combo = combo / gall; end
            if any(combo(ncols + 1:end) ~= 0)
                nxt = [nxt; combo]; %#ok<AGROW>
            end
        end
    end
    % support-minimality filter, applied at every step so the row set cannot
    % grow combinatorially on the way to the answer
    keep = true(size(nxt, 1), 1);
    supp = nxt(:, ncols + 1:end) ~= 0;
    for r = 1:size(nxt, 1)
        for s = 1:size(nxt, 1)
            if s == r || ~keep(r), continue; end
            if any(supp(s, :) & ~supp(r, :)), continue; end   % supp(s) not within supp(r)
            same = all(supp(s, :) == supp(r, :));
            if ~same || s < r, keep(r) = false; end            % keep the first of equal supports
        end
    end
    rows = nxt(keep, :);
end

S = rows(:, ncols + 1:end);
S = S(all(S >= 0, 2), :);          % an S-invariant is non-negative by definition

% ---- initial marking and the load vector V = S m0
if ~isempty(init)
    if numel(init) ~= n
        line_error(mfilename, 'init must hold one token count per place level');
    end
    m0 = zeros(1, n);
    for i = 1:n, m0(i) = i_asinteger(init(i), 'an initial marking'); end
else
    m0 = zeros(1, n);
    for k = 1:R
        if ~isfinite(sn.njobs(k))
            line_error(mfilename, sprintf(['class %d is open, so the net has no finite load ' ...
                'vector'], k));
        end
        ref = sn.refstat(k);
        for pp = 1:P
            if sn.nodeToStation(places(pp)) == ref
                m0((pp - 1) * R + k) = m0((pp - 1) * R + k) + ...
                    i_asinteger(sn.njobs(k), 'a class population');
            end
        end
    end
end
V = (S * m0(:))';

inv = struct('places', places, 'S', S, 'V', V, 'm0', m0);
end

% -------------------------------------------------------------------------
function v = i_asinteger(x, what)
% An arc multiplicity, refused unless integral.
r = floor(x + 0.5);
if abs(x - r) > 1e-9
    line_error('spn_sinvariants', sprintf(['%s is not integral; a fractional arc ' ...
        'multiplicity has no Petri-net meaning and no invariant basis over the integers'], what));
end
v = r;
end
