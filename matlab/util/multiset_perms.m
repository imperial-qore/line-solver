function P = multiset_perms(v)
% P = MULTISET_PERMS(v)
% All distinct permutations of the multiset held in vector v, one per row.
% Equivalent to unique(perms(v),'rows') but without materialising the n!
% intermediate, so a vector with repeated entries costs its multinomial
% count n!/prod(c_i!) rather than n!.
%
% ROW ORDER IS PART OF THE CONTRACT, not an implementation detail. State
% enumeration feeds these rows straight into sn.space, whose row 1 (after the
% end:-1:1 flip the builders apply) is the default initial state, so reordering
% them moves which state a chain starts in. Two orders are produced, and which
% one applies depends on the vector:
%   - all entries distinct: perms(v), i.e. MATLAB's own REVERSE-lexicographic
%     listing over positions (largest first);
%   - at least one repeat: grouped by leading value, taking the distinct values
%     in ASCENDING order, and recursively likewise for the remainder. The
%     recursion bottoms out in the two cases above, so a tail that happens to
%     be all-distinct is listed reverse-lexicographically inside its group.
% The mixture is inherited from the state space this function has always built;
% see cpp/include/line/lang/qn/state.h (pas_multiset_perms) and the python twin
% line_solver/api/state/multiset_perms.py, which reproduce the same listing.
%
% Callers in +State pass v assembled in ascending class order.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

v = v(:);
n = numel(v);

if n == 0
    P = [];
    return
end

u = unique(v); % distinct values, ascending
nu = numel(u);

if nu == 1
    % every arrangement of a single repeated value is the same one
    P = v.';
    return
end

if nu == n
    % no repeats: defer to MATLAB's own listing, reverse-lexicographic
    P = perms(v);
    return
end

% At least one value repeats. Count the rows up front so the result can be
% written into place: with c_i copies of the i-th distinct value the count is
% the multinomial n!/prod(c_i!), evaluated as a product of binomials so it
% stays an exact integer instead of overflowing through n!.
c = zeros(nu,1);
for i = 1:nu
    c(i) = sum(v == u(i));
end
nrows = 1;
left = n;
for i = 1:nu
    nrows = nrows * nchoosek(left, c(i));
    left = left - c(i);
end

P = zeros(nrows, n, 'like', v);
r = 0;
for i = 1:nu
    % strike the first copy of u(i); the rest keeps the order it came in
    rest = v;
    rest(find(rest == u(i), 1)) = [];
    tail = multiset_perms(rest);
    k = size(tail, 1);
    P(r+(1:k), 1) = u(i);
    P(r+(1:k), 2:n) = tail;
    r = r + k;
end
end
