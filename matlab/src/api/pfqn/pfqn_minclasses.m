%{
%{
 % @file pfqn_minclasses.m
 % @brief Lower bound on the number of customer classes needed to explain a
 %        measured sum of device utilizations, from Dowdy et al. (1992),
 %        J. ACM 39(1), Section 4.7.
%}
%}

function Rmin = pfqn_minclasses(Usum,K,N)
%{
%{
 % @brief Smallest number of customer classes R consistent with an observed sum
 %        of device utilizations, obtained by inverting the demand-free
 %        Expression (6) bound of pfqn_usumbound, which is nondecreasing in R.
 %        Only measured quantities are needed -- the utilizations, the device
 %        count and the population -- so the answer is available BEFORE any
 %        class-specific demand has been characterized, which is the point: it
 %        tells a clustering analysis how many classes it must at least find.
 %        An upper bound on R is meaningless (extra classes can always be
 %        introduced by splitting), so none is returned.
 %
 %        The paper's example: K = 2 devices, N = 3 customers, measured
 %        sum_k U_k = 1.6. A single class admits at most 2N/(N+1) = 1.5, so the
 %        single-class assumption is unjustified and Rmin = 2.
 % @fn pfqn_minclasses(Usum, K, N)
 % @param Usum Measured sum of device utilizations sum_k U_k (scalar).
 % @param K Number of devices (scalar).
 % @param N Total number of customers (scalar).
 % @return Rmin Least R in 1..N with pfqn_usumbound(R,K,N) >= Usum; NaN when
 %          Usum exceeds min(N,K) and so is unattainable by ANY class structure,
 %          which signals a measurement or bookkeeping error rather than a
 %          workload that needs more classes.
%}
%}
K = round(K); N = round(N);
if N < 1 || K < 1
    line_error(mfilename,'pfqn_minclasses requires N >= 1 and K >= 1.');
end
if Usum < 0
    line_error(mfilename,'pfqn_minclasses requires a nonnegative utilization sum.');
end
tol = 1e-12 * max(1,abs(Usum));
Rmin = NaN;
for R = 1:N
    if pfqn_usumbound(R,K,N) >= Usum - tol
        Rmin = R;
        return
    end
end
end
