%{
%{
 % @file pfqn_mmint2_gausslegendre.m
 % @brief McKenna-Mitra integral with Gauss-Legendre quadrature.
%}
%}

%{
%{
 % @brief McKenna-Mitra integral with Gauss-Legendre quadrature.
 % @fn pfqn_mmint2_gausslegendre(L, N, Z, m)
 % @param L Service demand vector.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param m Replication factor (default: 1).
 % @return G Normalizing constant.
 % @return lG Logarithm of normalizing constant.
%}
%}
function [G,lG]= pfqn_mmint2_gausslegendre(L,N,Z,m)
% [G,LOGG] = PFQN_MMINT2_GAUSSLEGENDRE(L,N,Z,m)
%
% Integrate McKenna-Mitra integral form with Gauss-Legendre in [0,1e6]
if nargin<4
    m=1; % multiplicity
end

persistent gausslegendreNodes;
persistent gausslegendreWeights;

% see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)

if isempty(gausslegendreNodes)
    [gausslegendreNodes, gausslegendreWeights] = load_gausslegendre_data();
end

% use at least 300 points
n = max(300,min(length(gausslegendreNodes),2*(sum(N)+m-1)-1));
y = zeros(1,n);
for i=1:n
    y(i)=N*log(Z+L*gausslegendreNodes(i))';
end
g = log(gausslegendreWeights(1:n))-gausslegendreNodes(1:n)+y(:);
coeff = - sum(factln(N))- factln(m-1) + (m-1)*sum(log(gausslegendreNodes(1:n)));
lG = log(sum(exp(g))) + coeff;
if ~isfinite(lG) % if numerical difficulties switch to logsumexp trick
    lG = logsumexp(g) + coeff;
end
G = exp(lG);
end

function [nodes, weights] = load_gausslegendre_data()
if coder.target('MATLAB')
    try
        data = load('gausslegendre-data.mat', 'gausslegendreNodes', 'gausslegendreWeights');
        nodes = data.gausslegendreNodes;
        weights = data.gausslegendreWeights;
    catch
        nodes = load(which('gausslegendre-nodes.txt'));
        weights = load(which('gausslegendre-weights.txt'));
    end
else
    data = coder.load('gausslegendre-data.mat', 'gausslegendreNodes', 'gausslegendreWeights');
    nodes = data.gausslegendreNodes;
    weights = data.gausslegendreWeights;
end
end
