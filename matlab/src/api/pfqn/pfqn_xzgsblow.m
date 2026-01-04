%{
%{
 % @file pfqn_xzgsblow.m
 % @brief Lower asymptotic bound on throughput (Zahorjan-Gittelsohn-Schweitzer-Bryant).
%}
%}

%{
%{
 % @brief Lower asymptotic bound on throughput (Zahorjan-Gittelsohn-Schweitzer-Bryant).
 % @fn pfqn_xzgsblow(L, N, Z)
 % @param L Service demand vector.
 % @param N Population.
 % @param Z Think time.
 % @return X Lower bound on throughput.
%}
%}
function X=pfqn_xzgsblow(L,N,Z)
M=length(L);
R=Z+sum(L)+max(L)*(N-1);
for ist=1:M
    if L(ist) < max(L)
        R=R+(L(ist)-max(L))*pfqn_qzgblow(L,N-1,Z,ist); 
    end
end
X=2*N*(1/(R+sqrt(R^2-4*Z*max(L)*(N-1))));

end