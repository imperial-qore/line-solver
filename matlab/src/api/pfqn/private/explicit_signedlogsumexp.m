%{
%{
 % @file explicit_signedlogsumexp.m
 % @brief Signed log-sum-exp shared by the explicit closed forms.
%}
%}

%{
%{
 % @brief Signed log-sum-exp of S = sum_i sterm(i)*exp(lterm(i)).
 %
 % Shared kernel of pfqn_explicit and pfqn_explicit_ld: both evaluate sums that
 % alternate in sign with terms far larger than the result, so the magnitudes
 % are carried in the log domain and the signs alongside them.
 %
 % @fn explicit_signedlogsumexp(lterm, sterm)
 % @param lterm Logarithms of the term magnitudes.
 % @param sterm Signs of the terms, zero for a term that drops out.
 % @return lS Logarithm of |S|.
 % @return sgnS Sign of S, zero when every term dropped out.
 % @return lossDigits Decimal digits lost to cancellation.
%}
%}
function [lS,sgnS,lossDigits] = explicit_signedlogsumexp(lterm,sterm)
keep = isfinite(lterm) & sterm~=0;
if ~any(keep)
    lS = -Inf; sgnS = 0; lossDigits = 0;
    return
end
lterm = lterm(keep);
sterm = sterm(keep);
a = max(lterm);
s = sum(sterm .* exp(lterm - a));
sgnS = sign(s);
if s==0
    lS = -Inf; lossDigits = Inf;
    return
end
lS = a + log(abs(s));
% max(exp(lterm-a)) is 1, so -log10|s| is the shortfall of the sum against its
% largest term. Every one of the n terms carries a rounding error of order
% eps*max_term, so the accumulated absolute error is n*eps*max_term and the
% digits actually lost are that shortfall PLUS log10(n). Dropping the count
% understates the loss by log10(n) and lets a wrong answer past the guard.
lossDigits = max(0,log10(numel(sterm)/abs(s)));
end
