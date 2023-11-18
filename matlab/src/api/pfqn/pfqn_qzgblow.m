function Qgb=pfqn_qzgblow(L,N,Z,i)
yi=N*L(i)/(Z+sum(L)+max(L)*N);
Qgb=yi/(1-yi) - (yi^(N+1))/(1-yi);
end