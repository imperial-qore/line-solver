%{
%{
 % @file explicit_gdistinct.m
 % @brief Single-class fixed-rate normalizing constant at pairwise distinct demands.
%}
%}

%{
%{
 % @brief Single-class fixed-rate normalizing constant at pairwise distinct demands.
 %
 % Eq. (14) of Casale, "Accelerating Performance Inference over Closed Systems
 % by Asymptotic Methods", ACM SIGMETRICS 2017: Gordon's partial-fraction form
 %
 %    g(Nt) = sum_k th_k^(Nt+K-1) / prod_{i~=k} (th_k - th_i)
 %
 % A zero demand contributes nothing, which also realizes the 0/0 = 0 convention
 % of Eq. (15) when the zero is repeated. Shared by pfqn_explicit, which reads it
 % at the induced demands of a fixed-rate model, and by pfqn_explicit_ld, which
 % reads it at the SCALED demands th_k/alpha_k(s_k) of a load-dependent one.
 %
 % @fn explicit_gdistinct(th, Nt, K)
 % @param th Demands (Kx1), pairwise distinct.
 % @param Nt Total population.
 % @param K Number of queues.
 % @return lg Logarithm of |g|.
 % @return sg Sign of g.
 % @return lossDigits Decimal digits lost to cancellation.
%}
%}
function [lg,sg,lossDigits] = explicit_gdistinct(th,Nt,K)
lin = -Inf(K,1);
sgv = zeros(K,1);
for k=1:K
    if th(k)<=0
        continue
    end
    d = th(k) - th([1:(k-1),(k+1):K]);
    lin(k) = (Nt+K-1)*log(th(k)) - sum(log(abs(d)));
    sgv(k) = prod(sign(d));
end
[lg,sg,lossDigits] = explicit_signedlogsumexp(lin,sgv);
end
