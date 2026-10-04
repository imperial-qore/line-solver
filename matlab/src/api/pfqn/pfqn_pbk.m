%{
%{
 % @file pfqn_pbk.m
 % @brief Iterative PB(k) proportional bounds (Eager-Sevcik) for single-class
 %        closed networks with delay.
%}
%}

function [Xlo,Xhi] = pfqn_pbk(L,N,Z,k)
%{
%{
 % @brief PB(k): the Eager-Sevcik proportional (performance) bound at
 %        iteration count k, i.e. the level-k member of the Performance Bound
 %        Hierarchy computed at populations N, N-1, ..., N-k (Casale-Muntz-
 %        Serazzi 2008 cite this as the iterative extension used in Tables 5,7).
 %        Reduces to the noniterative asymptotic bound at k=0 and tightens to
 %        exact as k->N. Backed by the validated PBH recursion (pfqn_pbh).
 % @fn pfqn_pbk(L, N, Z, k)
 % @param L Service demand vector (M x 1).
 % @param N Population (scalar).
 % @param Z Think time (scalar, default 0).
 % @param k Iteration count >= 0 (default 1).
 % @return Xlo Lower throughput bound.
 % @return Xhi Upper throughput bound.
%}
%}
if nargin < 3 || isempty(Z), Z = 0; end
if nargin < 4 || isempty(k), k = 1; end
[Xlo,Xhi] = pfqn_pbh(L,N,Z,k);
end
