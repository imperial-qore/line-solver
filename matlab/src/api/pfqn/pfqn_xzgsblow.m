function X=pfqn_xzgsblow(L,N,Z)
M=length(L);
R=Z+sum(L)+max(L)*(N-1);
for i=1:M
    if L(i) < max(L)
        R=R+(L(i)-max(L))*pfqn_qzgblow(L,N-1,Z,i); 
    end
end
X=2*N*(1/(R+sqrt(R^2-4*Z*max(L)*(N-1))));

end