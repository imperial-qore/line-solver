%{
%{
 % @file pfqn_qzgblow.m
 % @brief Lower asymptotic bound on queue length (Zahorjan-Gittelsohn-Bryant).
%}
%}

%{
%{
 % @brief Lower asymptotic bound on queue length (Zahorjan-Gittelsohn-Bryant).
 % @fn pfqn_qzgblow(L, N, Z, i)
 % @param L Service demand vector.
 % @param N Population.
 % @param Z Think time.
 % @param i Station index.
 % @return Qgb Lower bound on mean queue length at station i.
%}
%}
function Qgb=pfqn_qzgblow(L,N,Z,i)
yi=N*L(i)/(Z+sum(L)+max(L)*N);
Qgb=yi/(1-yi) - (yi^(N+1))/(1-yi);
end