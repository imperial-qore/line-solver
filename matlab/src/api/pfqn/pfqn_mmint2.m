function [G,lG]= pfqn_mmint2(L,N,Z,method)
% [G,LOGG] = PFQN_PNC2(L,N,Z)

if nargin<4
    method = 'default';
end

nnzClasses = find(N);
% repairmen integration
order = 12;
% below we use a variable substitution u->u^2 as it is numerically better
switch method
    case 'quadratic'
        f= @(u) (2*u'.*exp(-(u.^2)').*prod((Z(nnzClasses)+L(nnzClasses).*repmat(u(:).^2,1,length(nnzClasses))).^N(nnzClasses),2))';
    case 'default'
        f= @(u) (exp(-u').*prod((Z(nnzClasses)+L(nnzClasses).*repmat(u(:),1,length(nnzClasses))).^N(nnzClasses),2))';
end

p = 1-10^-order;
exp1prctile = -log(1-p)/1; % cutoff for exponential term
w = warning ;
warning off;
lG = log(integral(f,0,exp1prctile,'AbsTol',10^-order)) - sum(factln(N));
G = exp(lG);
warning(w);
end
