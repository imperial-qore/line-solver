%{
%{
 % @file infradius_h.m
 % @brief Helper function for infinite radius computation with logistic transformation.
%}
%}

%{
%{
 % @brief Helper function for infinite radius computation with logistic transformation.
 % @fn infradius_h(x, L, N, alpha)
 % @param x Logistic transformation parameters.
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param alpha Load-dependent rate matrix.
 % @return y Evaluated function value for integration.
%}
%}
function y=infradius_h(x,L,N,alpha)
M = size(L,1);
Nt = sum(N);
beta = N/Nt;
t = exp(x)./(1+exp(x));
tb = sum(beta.*t,2);
h = @(x) pfqn_gld(sum(L.*repmat(exp(2*pi*1i*(t-tb)),M,1),2),Nt,alpha) * prod(exp(x)./(1+exp(x)).^2,2);

y = zeros(size(x,1),1);
for i=1:size(x,1)
    y(i) = real(h(x(i,:)));
end
end
