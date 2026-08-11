%{
%{
 % @file pfqn_bjbk.m
 % @brief Iterative BJB(k) balanced job bounds (Zahorjan et al.) for
 %        single-class closed networks with delay.
%}
%}

function [Xlo,Xhi] = pfqn_bjbk(L,N,Z,k)
%{
%{
 % @brief BJB(k): the iterative Balanced Job Bound at iteration count k
 %        (Casale-Muntz-Serazzi 2008, Tables 5,7). BJB(1) recovers the
 %        noniterative Balanced Job Bound; each additional iteration performs
 %        one exact MVA step from the balanced seed at N-k, and the bracket
 %        tightens to exact as k->N. Realized through the validated PBH
 %        recursion (pfqn_pbh), whose level-1 optimistic bound equals the
 %        noniterative BJB optimistic bound (Eager-Sevcik 1983, eq. 10).
 % @fn pfqn_bjbk(L, N, Z, k)
 % @param L Service demand vector (M x 1).
 % @param N Population (scalar).
 % @param Z Think time (scalar, default 0).
 % @param k Iteration count >= 1 (default 1).
 % @return Xlo Lower throughput bound.
 % @return Xhi Upper throughput bound.
%}
%}
if nargin < 3 || isempty(Z), Z = 0; end
if nargin < 4 || isempty(k), k = 1; end
[Xlo,Xhi] = pfqn_pbh(L,N,Z,k);
end
