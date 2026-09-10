%{
%{
 % @file sjn_cap.m
 % @brief Capacity constraint on the utilization of a shortest-job-next station.
%}
%}

function [C, X, kappa, bound] = sjn_cap(caller, C, L, N, Z, sjnset, umax)
%{
%{
 % @brief Enforce U <= umax at every shortest-job-next station by inflating
 %        its waiting time, and return the resulting residence times and
 %        throughputs.
 %
 %        The SJN response time equation is an open-system one: its
 %        denominator is 1 - U(x), the fraction of the server taken by jobs no
 %        longer than x, and it has no solution once that reaches one. A
 %        closed network never reaches it in reality, but the approximation
 %        can, because it underestimates the residence time at a congested SJN
 %        station and the resulting throughput then exceeds the station
 %        capacity 1/max_r L. What is imposed here is the utilization law
 %        sum_r X_r L_mr <= umax, an exact property of the network and not a
 %        property of the approximation.
 %
 %        The constraint is imposed on the waiting time rather than on the
 %        throughput, i.e. by scaling the excess C - L of the station by a
 %        factor kappa >= 1 found by bisection. Throughput is then recomputed
 %        from the residence times, so that X (Z + sum_m C) = N still holds
 %        exactly and no jobs are lost, which capping X directly would break.
 %        The same kappa scales the conditional waiting time profile of the
 %        station, every quantity derived from it being linear in it.
 %
 %        A binding cap means the station is in the starvation regime, where
 %        long jobs are held back indefinitely and the arrival theorem is
 %        badly violated. The caller is expected to warn: the returned values
 %        are stable but their accuracy is not warranted there.
 % @fn sjn_cap(caller, C, L, N, Z, sjnset, umax)
 % @param caller Name of the calling function, used in error messages.
 % @param C Residence times (M x R).
 % @param L Service demands (M x R).
 % @param N Population vector (1 x R).
 % @param Z Think times (1 x R).
 % @param sjnset Indices of the SJN stations.
 % @param umax Utilization cap, strictly below one.
 % @return C Residence times after the cap.
 % @return X Throughputs consistent with the returned C (1 x R).
 % @return kappa Waiting time inflation applied at each SJN station (1 x nsjn).
 % @return bound True if the cap was binding at any station.
%}
%}
nsjn = length(sjnset);
kappa = ones(1,nsjn);
bound = false;
X = local_thru(C, N, Z);
if nsjn == 0
    return
end
if umax >= 1
    line_error(caller,'the utilization cap must be strictly below one, the response time equation is singular at one');
end
for sweep = 1:20
    viol = false;
    for q = 1:nsjn
        m = sjnset(q);
        if sum(X .* L(m,:)) <= umax
            continue
        end
        viol = true;
        bound = true;
        Wq = C(m,:) - L(m,:);
        hi = 2;
        while local_rho(C, L, m, N, Z, Wq, hi) > umax
            hi = 2*hi;
            if hi > 1e12
                line_error(caller,sprintf(['station %d cannot be brought under the utilization cap by any waiting time:\n' ...
                    'its service demands alone saturate it at this population.'],m));
            end
        end
        lo = 1;
        for b = 1:200
            mid = (lo + hi)/2;
            if local_rho(C, L, m, N, Z, Wq, mid) > umax
                lo = mid;
            else
                hi = mid;
            end
        end
        kappa(q) = kappa(q) * hi;
        C(m,:) = L(m,:) + hi * Wq;
        X = local_thru(C, N, Z);
    end
    if ~viol
        return
    end
end
line_error(caller,'the utilization cap did not settle across the SJN stations');
end

function X = local_thru(C, N, Z)
den = Z + sum(C,1);
X = zeros(1,length(N));
act = N > 0;
X(act) = N(act) ./ den(act);
end

function rho = local_rho(C, L, m, N, Z, Wq, kappa)
C(m,:) = L(m,:) + kappa * Wq;
rho = sum(local_thru(C, N, Z) .* L(m,:));
end
