function [G,lG]= pfqn_mmint2_gausslaguerre(L,N,Z,m)
% [G,LOGG] = PFQN_MMINT2_GAUSSLAGUERRE(L,N,Z,m)
%
% Integrate with Gauss-Laguerre

if nargin<4
    m=1;
end

persistent gausslaguerreNodes;
persistent gausslaguerreWeights;

if isempty(gausslaguerreNodes)
    [ gausslaguerreNodes,gausslaguerreWeights ] = gengausslegquadrule(300,10^-5);
end

lambda = 0*N;
nonzeroClasses = find(N);

% repairmen integration
f= @(u) N(nonzeroClasses)*log(Z(nonzeroClasses)+L(nonzeroClasses)*u)';
x = gausslaguerreNodes;
w = gausslaguerreWeights;
n = min(300,2*sum(N)+1);
F = zeros(size(x));
for i=1:length(x)
    F(i)=(m-1)*log(x(i))+f(x(i));
end
g = log(w) + F - sum(factln(N))- factln(m-1);
lG = log(sum(exp(g)));
if ~isfinite(lG) % if numerical difficulties switch to logsumexp trick
    lG = logsumexp(g);
end
G = exp(lG);
end