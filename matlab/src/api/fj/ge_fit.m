%{ @file ge_fit.m
 %  @brief Two-moment fit of a generalized exponential law
 %
 %  @author LINE Development Team
%}

%{
 % @brief Two-moment fit of a generalized exponential law
 %
 % @details
 % Matches the generalized exponential CDF
 %
 %   F(x) = (1 - exp(-x/beta))^alpha
 %
 % on a mean and a variance,
 %
 %   E[T] = beta*(psi(alpha+1) - psi(1)),  V[T] = beta^2*(psi'(1) - psi'(alpha+1)).
 %
 % The squared coefficient of variation depends on the SHAPE ALONE and
 % decreases monotonically in it, so the shape is recovered by a scalar
 % root-find on a logarithmic scale and the scale then follows in closed
 % form. SCV = 1 is the exponential case alpha = 1, kept exact.
 %
 % It is the per-branch law of the ForkTail approximation, shared by
 % FJ_TAIL_FORKTAIL (the AND-join) and FJ_TAIL_ORDSTAT (the k-of-n quorum),
 % which is why it lives in a file of its own rather than inside either.
 %
 % @par Syntax:
 % @code
 % [alpha, beta] = ge_fit(ET, VT)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>ET<td>Mean of the law, positive
 % <tr><td>VT<td>Variance of the law, positive
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>alpha<td>Fitted shape parameter
 % <tr><td>beta<td>Fitted scale parameter
 % </table>
 %
 % @par References:
 % M. Nguyen, S. Alesawi, N. Li, H. Che, H. Jiang, "ForkTail: A Black-Box
 % Fork-Join Tail Latency Prediction Model for User-Facing Datacenter
 % Workloads", ACM HPDC 2018, pp. 206-217.
%}
function [alpha, beta] = ge_fit(ET, VT)

scv = VT / ET^2;
if abs(scv - 1) < GlobalConstants.FineTol
    alpha = 1;
else
    scvOf = @(a) (psi(1,1) - psi(1,a+1)) / (psi(0,a+1) - psi(0,1))^2;
    residual = @(la) scvOf(exp(la)) - scv;
    lo = -30; hi = 30;
    while residual(lo) < 0 && lo > -700
        lo = lo - 30;   % smaller shape -> larger SCV
    end
    while residual(hi) > 0 && hi < 700
        hi = hi + 30;   % larger shape -> smaller SCV
    end
    alpha = exp(fzero(residual, [lo, hi]));
end
beta = ET / (psi(0,alpha+1) - psi(0,1));
end
