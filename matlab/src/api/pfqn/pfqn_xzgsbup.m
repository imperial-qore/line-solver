%{
%{
 % @file pfqn_xzgsbup.m
 % @brief Upper asymptotic bound on throughput (Zahorjan-Gittelsohn-Schweitzer-Bryant).
%}
%}

%{
%{
 % @brief Upper asymptotic bound on throughput (Zahorjan-Gittelsohn-Schweitzer-Bryant).
 % @fn pfqn_xzgsbup(L, N, Z)
 % @param L Service demand vector.
 % @param N Population.
 % @param Z Think time.
 % @return X Upper bound on throughput.
%}
%}
function X=pfqn_xzgsbup(L,N,Z)
M=length(L);
Lmax = max(L);
R=Z+sum(L)+Lmax*(N-1);
for ist=1:M
    if L(ist) < Lmax
        R=R+(L(ist)-Lmax)*pfqn_qzgbup(L,N-1,Z,ist);
    end
end
discriminant = R^2 - 4*Z*Lmax*N;
if discriminant < 0
    discriminant = 0;
end
X=2*N/(R+sqrt(discriminant));
end
