%{
%{
 % @file pfqn_scbgap.m
 % @brief Maximum relative throughput error of aggregating r customer classes,
 %        from Dowdy et al. (1992), J. ACM 39(1), Expressions (3)-(5).
%}
%}

function e = pfqn_scbgap(N,K,r,undominated)
%{
%{
 % @brief Demand-free bound on the relative throughput error incurred when r of
 %        the N single-customer classes of a closed product-form network are
 %        merged into one class. With r = N (the default) this is the full
 %        single-class aggregation error of their Theorem 3, and the error is
 %        then at most 50%; with r < N it is the partial-aggregation error of
 %        their Theorem 4.
 %
 %        The bound depends only on N, K and r, never on the demands, so it can
 %        be attached as a certified error bar to any result computed on merged
 %        chains -- LINE merges classes into chains routinely through
 %        sn_get_demands_chain.
 %
 %        General case, dominating classes allowed (Expression 4, and with
 %        r = N Expression 3):
 %            e = (min(r,K)-1) / (r + min(r,K) - 1).
 %        Undominated case, every customer placing the same total demand
 %        (Theorem 5 and its comment (3), which lifts the N = R restriction):
 %            e = r(r-1) / (min(N,K)(2r-1)),      valid for r <= K only,
 %        smaller than the general case by the factor r/min(N,K) and equal to it
 %        at r = K. THE DOMAIN IS NOT COSMETIC: Theorem 5 gives each of its R
 %        classes a dedicated device, so r never exceeds K there, and comment (3)
 %        states the generalization for r < K. Evaluated at r > K the expression
 %        climbs past the general bound and past the 50% cap of Theorem 3, i.e.
 %        it stops being a bound, so r > K is refused rather than returned.
 % @fn pfqn_scbgap(N, K, r, undominated)
 % @param N Total number of customers (scalar). One customer per class, so N is
 %          also the number of classes before merging.
 % @param K Number of queueing devices (scalar).
 % @param r Number of classes merged into one (default N, full aggregation).
 % @param undominated True to use the tighter Theorem-5 form, valid only when no
 %          class dominates, i.e. every customer's total device demand is equal,
 %          and only for r <= K (default false).
 % @return e Maximum relative throughput error, in [0,1/2].
%}
%}
if nargin < 3 || isempty(r), r = N; end
if nargin < 4 || isempty(undominated), undominated = false; end
N = round(N); K = round(K); r = round(r);
if N < 1 || K < 1
    line_error(mfilename,'pfqn_scbgap requires N >= 1 and K >= 1.');
end
if r < 1 || r > N
    line_error(mfilename,'pfqn_scbgap requires 1 <= r <= N (r=%d, N=%d).', r, N);
end
if r == 1
    e = 0;                                  % merging one class changes nothing
    return
end
if undominated
    if r > K
        line_error(mfilename,'The undominated (Theorem 5) form is defined for r <= K only (r=%d, K=%d); beyond it the expression exceeds the general bound and the 50%% cap.', r, K);
    end
    e = r*(r-1) / (min(N,K) * (2*r-1));     % Theorem 5, comment (3)
else
    m = min(r,K);
    e = (m-1) / (r + m - 1);                % Expression (4); r=N gives (3)
end
end
