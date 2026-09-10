%{
%{
 % @file pfqn_mmint2_gausslaguerre.m
 % @brief McKenna-Mitra integral with Gauss-Laguerre quadrature.
%}
%}

%{
%{
 % @brief McKenna-Mitra integral with Gauss-Laguerre quadrature.
 % @fn pfqn_mmint2_gausslaguerre(L, N, Z, m)
 % @param L Service demand vector.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param m Replication factor (default: 1).
 % @return G Normalizing constant.
 % @return lG Logarithm of normalizing constant.
%}
%}
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
    [gausslaguerreNodes, gausslaguerreWeights] = load_gausslaguerre_data();
end

x = gausslaguerreNodes;
w = gausslaguerreWeights;
npts = length(x);
F = zeros(1,npts);
for i=1:npts
    F(i)=(m-1)*log(x(i));
    for r=1:length(N)
        if N(r) ~= 0
            F(i) = F(i) + N(r) * log(Z(r)+L(r)*x(i));
        end
    end
end
g = log(w) + F - sum(factln(N))- factln(m-1);
lG = log(sum(exp(g)));
if ~isfinite(lG) % if numerical difficulties switch to logsumexp trick
    lG = logsumexp(g);
end
G = exp(lG);
end

function [nodes, weights] = load_gausslaguerre_data()
if coder.target('MATLAB')
    data = load('gausslaguerre-data.mat', 'gausslaguerreNodes', 'gausslaguerreWeights');
else
    data = coder.load('gausslaguerre-data.mat', 'gausslaguerreNodes', 'gausslaguerreWeights');
end
nodes = data.gausslaguerreNodes;
weights = data.gausslaguerreWeights;
end
