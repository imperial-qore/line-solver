function perm_require_support(A, caller)
% PERM_REQUIRE_SUPPORT  Refuse a matrix the permanent approximations cannot take.
%
% PERM_REQUIRE_SUPPORT(A, CALLER) raises an error if A has an entry that is not
% strictly positive. The four approximate permanent engines -- PERM_BETHE,
% PERM_HEUR, PERM_HUBERLAW and PERM_ADAPART -- all rest, directly or through
% the Sinkhorn scaling they share, on a strictly positive matrix.
%
% WHY THIS IS REFUSED RATHER THAN FLOORED. Each of these routines used to
% replace a zero by a small constant EPSILON before doing anything else. That
% substitution is not invertible. Every permutation picks exactly one entry
% from each row, so a matrix with an identically zero row has permanent 0 while
% the floored matrix has permanent N!*EPSILON times the permanent of the rest.
% N! outruns EPSILON quickly: with EPSILON = 2.22e-16 the fabricated value
% passes 1% at N = 17 and 1 at N = 18, and at N = 20 the floor alone
% manufactures a permanent of about 540 where the truth is exactly zero. The
% order of the replicated demand matrix in PFQN_JOINTMARG is SUM(N), so a
% closed model with 18 jobs and one structurally zero demand is already in that
% regime.
%
% POSITIVITY IS SUFFICIENT BUT NOT NECESSARY. The sharp precondition of the
% Sinkhorn scaling is TOTAL SUPPORT -- every positive entry lies on a positive
% permutation -- which a matrix with a strictly positive permanent can still
% fail: [J3 0; J3 J3] has permanent 36, no total support, and the scaling exits
% on its tolerance rather than on convergence. Positivity is the test used here
% because it is O(N^2), it is the contract the C++ header already states, and
% it is the guard PFQN_JOINTMARG already applies one level up.
%
% Use PERM or PFQN_PERM instead: the exact engine handles zeros correctly, a
% station holding no jobs contributes no row, a class holding no jobs no
% column, and the permanent of the empty matrix is 1.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2 || isempty(caller)
    caller = 'perm_require_support';
end

bad = find(~(A > 0), 1);
if ~isempty(bad)
    [i0, j0] = ind2sub(size(A), bad);
    line_error(caller, ['The ''%s'' permanent approximation requires a strictly ' ...
        'positive matrix: entry (%d,%d) is %g, so the matrix has no full support. ' ...
        'Flooring it would change the permanent by n!*eps, which is O(1) by n=18. ' ...
        'Use the exact engine (perm or pfqn_perm).'], caller, i0, j0, full(A(bad)));
end
end
