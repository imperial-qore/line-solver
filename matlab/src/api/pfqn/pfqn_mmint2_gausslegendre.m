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

% nodes and weights generated with tridiagonal eigenvalues method in 
% high-precision using Julia:
%
% using LinearAlgebra
%
% function gauss(a, b, N)
%    λ, Q = eigen(SymTridiagonal(zeros(N), [n / sqrt(4n^2 - 1) for n = 1:N-1]))
%    @. (λ + 1) * (b - a) / 2 + a, [2Q[1, i]^2 for i = 1:N] * (b - a) / 2
% end

if isempty(gausslegendreNodes)
    gausslegendreNodes=load(which('gausslegendre-nodes.txt'));
    gausslegendreWeights=load(which('gausslegendre-weights.txt'));
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