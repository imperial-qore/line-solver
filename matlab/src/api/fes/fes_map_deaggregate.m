%{ @file fes_map_deaggregate.m
 %  @brief Per-station metrics behind a MAP flow-equivalent server
 %
 %  @author LINE Development Team
%}

%{
 % @brief Recovers the per-station metrics of an aggregated subnetwork by
 % conditioning on the population held by the flow-equivalent server
 %
 % @details
 % fes_map_solve returns the distribution pk of the jobs held by the
 % aggregate. The metrics of the stations behind it follow by conditioning,
 % E[Y_i] = sum_k pk(k) Y_i(k), with Y_i(k) the metric of station i when
 % the isolated subnetwork holds k jobs. This is the decomposition step of
 % the hierarchical analysis of Chandy, Herzog and Woo, IBM J. Res. Dev.
 % 19(1), 1975, and it is exact for a product-form subnetwork. It is an
 % approximation whenever the burstiness that the MAP flow-equivalent
 % server carries also matters inside the subnetwork, because the
 % conditional solve is the product-form one; the aggregate metrics
 % returned by fes_map_solve do not rely on it.
 %
 % @par Syntax:
 % @code
 % [QN,UN,XN,RN] = fes_map_deaggregate(pk, L, mi, isDelay)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>pk<td>Distribution of the jobs held by the aggregate, pk(k+1) = P(k)
 % <tr><td>L<td>(M_sub x 1) service demands of the isolated subnetwork
 % <tr><td>mi<td>(1 x M_sub) servers per station, Inf for a delay
 % <tr><td>isDelay<td>(1 x M_sub) true where the station is a pure delay
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>QN<td>(M_sub x 1) mean queue length per station
 % <tr><td>UN<td>(M_sub x 1) utilization per station
 % <tr><td>XN<td>(M_sub x 1) throughput per station
 % <tr><td>RN<td>(M_sub x 1) mean residence time per station
 % </table>
 %
 % @see fes_map_solve, fes_map_aggregate
%}
function [QN,UN,XN,RN] = fes_map_deaggregate(pk, L, mi, isDelay)

M = numel(L);
N = numel(pk) - 1;
QN = zeros(M,1);
UN = zeros(M,1);
XN = zeros(M,1);

queueIdx = find(~isDelay);
delayIdx = find(isDelay);
Lq = L(queueIdx);
miq = mi(queueIdx);
Z = sum(L(delayIdx));

for k = 1:N
    if pk(k+1) <= GlobalConstants.Zero
        continue
    end
    % MI is the additive C=L*(mi+Qarv) term, not a server count: multiservers
    % go through PFQN_MVAMS, as in FES_COMPUTE_THROUGHPUTS. Its UN is per
    % STATION on the multiserver branch and reports 1-P(0) rather than the
    % [0,1] load, so utilization is recomputed here as U=X*L/S.
    [Xk,Qk] = pfqn_mvams(0, Lq(:), k, Z, ones(numel(Lq),1), miq(:));
    Uk = Xk(1) * Lq(:) ./ max(1, miq(:));
    QN(queueIdx) = QN(queueIdx) + pk(k+1)*Qk(:);
    UN(queueIdx) = UN(queueIdx) + pk(k+1)*Uk(:);
    XN(queueIdx) = XN(queueIdx) + pk(k+1)*Xk;
    for d = delayIdx
        QN(d) = QN(d) + pk(k+1)*Xk*L(d);
        UN(d) = UN(d) + pk(k+1)*Xk*L(d);
        XN(d) = XN(d) + pk(k+1)*Xk;
    end
end

RN = zeros(M,1);
nz = XN > GlobalConstants.Zero;
RN(nz) = QN(nz)./XN(nz);
end
