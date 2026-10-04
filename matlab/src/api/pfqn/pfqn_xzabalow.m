%{
%{
 % @file pfqn_xzabalow.m
 % @brief Lower ABA (asymptotic bound analysis) bound on throughput. Not the Zahorjan-Balanced bound (that is pfqn_xzgsblow); single-class only.
%}
%}

%{
%{
 % @brief Lower ABA (asymptotic bound analysis) bound on throughput. Not the Zahorjan-Balanced bound (that is pfqn_xzgsblow); single-class only.
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