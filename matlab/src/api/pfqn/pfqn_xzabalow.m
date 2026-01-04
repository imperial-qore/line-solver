%{
%{
 % @file pfqn_xzabalow.m
 % @brief Lower asymptotic bound on throughput (Zahorjan-Balanced).
%}
%}

%{
%{
 % @brief Lower asymptotic bound on throughput (Zahorjan-Balanced).
 % @fn pfqn_xzabalow(L, N, Z)
 % @param L Service demand vector.
 % @param N Population.
 % @param Z Think time.
 % @return XN Lower bound on throughput.
%}
%}
function [XN]=pfqn_xzabalow(L,N,Z)
    Ltot=sum(L);
    XN=N/(Z+Ltot*N);
end