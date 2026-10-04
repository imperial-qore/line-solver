%{
%{
 % @file pfqn_usumbound.m
 % @brief Upper bound on the sum of device utilizations of a closed R-class
 %        network, from Dowdy et al. (1992), J. ACM 39(1), Expression (6).
%}
%}

function Umax = pfqn_usumbound(R,K,N)
%{
%{
 % @brief Largest value the sum of device utilizations sum_k U_k,R can take in
 %        any closed product-form network with R classes, K devices and N
 %        customers (their Theorem 6):
 %            sum_k U_k,R <= (H-1) + (K-H+1)(N-H+1)/(K+N-2H+1),  H = min(R,K).
 %        The bound is demand-free and nondecreasing in R, which is what makes
 %        it invertible into a lower bound on the number of necessary customer
 %        classes; see pfqn_minclasses. At R >= min(N,K) it reaches min(N,K),
 %        the trivial cap of one busy server per device.
 %
 %        The paper's worked case is K = 2, N = 3, R = 1, giving 2N/(N+1) = 1.5:
 %        a measured sum of 1.6 then refutes the single-class assumption.
 % @fn pfqn_usumbound(R, K, N)
 % @param R Number of single-customer classes (scalar, 1 <= R <= N).
 % @param K Number of devices (scalar).
 % @param N Total number of customers (scalar).
 % @return Umax Upper bound on sum_k U_k,R.
%}
%}
R = round(R); K = round(K); N = round(N);
if N < 1 || K < 1
    line_error(mfilename,'pfqn_usumbound requires N >= 1 and K >= 1.');
end
if R < 1 || R > N
    line_error(mfilename,'pfqn_usumbound requires 1 <= R <= N (R=%d, N=%d).', R, N);
end
H = min(R,K);
Umax = (H-1) + (K-H+1)*(N-H+1) / (K+N-2*H+1);
end
