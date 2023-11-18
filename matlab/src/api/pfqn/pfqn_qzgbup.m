function Qgb=pfqn_qzgbup(L,N,Z,i)
sigma=sum(L.^2)/sum(L);
Yi=L(i)*min([1/max(L),N/(Z+sum(L)+sigma*(N-1-Z*pfqn_xzabaup(L,N-1,Z)))]);
if Yi<1
    Qgb=Yi/(1-Yi) - (Yi^(N+1))/(1-Yi);
else
    Qgb=N;
end
end