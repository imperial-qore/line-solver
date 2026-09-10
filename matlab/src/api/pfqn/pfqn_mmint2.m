%{
%{
 % @file pfqn_mmint2.m
 % @brief McKenna-Mitra integral form for repairman models using MATLAB integral.
%}
%}

%{
%{
 % @brief McKenna-Mitra integral form for repairman models using MATLAB integral.
 % @fn pfqn_mmint2(L, N, Z)
 % @param L Service demand vector.
 % @param N Population vector.
 % @param Z Think time vector.
 % @return G Normalizing constant.
 % @return lG Logarithm of normalizing constant.
%}
%}
function [G,lG]= pfqn_mmint2(L,N,Z)
% [G,LOGG] = PFQN_MMINT2(L,N,Z)

nnzClasses = find(N);
% repairmen integration
order = 12;
f= @(u) (exp(-u').*prod((Z(nnzClasses)+L(nnzClasses).*repmat(u(:),1,length(nnzClasses))).^N(nnzClasses),2))';

p = 1-10^-order;
exp1prctile = -log(1-p)/1; % cutoff for exponential term
w = warning ;
warning off;
lG = log(integral(f,0,exp1prctile,'AbsTol',10^-order)) - sum(factln(N));
G = exp(lG);
warning(w);
end
