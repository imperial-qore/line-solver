%{
%{
 % @file pfqn_perm.m
 % @brief Permanent of a demand matrix, with optional column multiplicities.
%}
%}

%{
%{
 % @brief Permanent of a demand matrix, with optional column multiplicities.
 % @fn pfqn_perm(A, m)
 % @param A Matrix whose permanent is required.
 % @param m Multiplicity of each column of A (optional).
 % @return val Permanent value.
%}
%}
function val = pfqn_perm(A, m)
% VAL = PFQN_PERM(A)     permanent of the square matrix A
% VAL = PFQN_PERM(A, M)  permanent of the matrix whose column J is column J of
%                        A repeated M(J) times, so that SUM(M) == SIZE(A,1)
%
% The pfqn_ entry point of the permanent library. It exists because the
% product-form joint queue-length probability of the per-station TOTAL
% populations is a permanent of the demand matrix replicated once per job,
% which is a normalizing-constant quantity rather than a general-purpose
% linear algebra one; see pfqn_jointmarg.
%
% Orientation is chosen before repeated lines are grouped. That is a
% correctness concern, not an optimisation: perm(A) is transpose-invariant but
% the Ryser SUM is not, and exploiting repeated rows silently expands the
% transpose. See matlab/util/perm.m and _kb/03-api-layer.md.
%
% Reference:
%   H. J. Ryser, "Combinatorial Mathematics", Carus Mathematical Monographs
%   14, Mathematical Association of America, 1963.

if nargin < 2
    val = perm(A);
else
    val = perm(A, m);
end
end
