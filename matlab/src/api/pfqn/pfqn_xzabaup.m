%{
%{
 % @file pfqn_xzabaup.m
 % @brief Upper asymptotic bound on throughput (Zahorjan-Balanced).
%}
%}

%{
%{
 % @brief Upper asymptotic bound on throughput (Zahorjan-Balanced).
 % @fn pfqn_xzabaup(L, N, Z)
 % @param L Service demand vector.
 % @param N Population.
 % @param Z Think time.
 % @return XN Upper bound on throughput.
%}
%}
function [XN]=pfqn_xzabaup(L,N,Z)
    XN=min([1/max(L) N/(sum(L)+Z)]);
end