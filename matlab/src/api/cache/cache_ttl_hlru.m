%{ @file cache_ttl_hlru.m
 %  @brief TTL (characteristic-time) approximation for h-LRU / LRU(m) caches
 %
 %  @author LINE Development Team
%}

%{
 % @brief Steady-state list occupancy probabilities for an h-LRU cache
 %
 % @details
 % Computes the characteristic-time (TTL) approximation of the list-based
 % h-LRU (LRU(m)) replacement policy: h LRU lists of capacities m(1..h), a
 % miss inserts the item at the head of list 1, a hit in list l exchanges
 % the item with the tail of list l+1 (Gast and Van Houdt, SIGMETRICS 2015).
 % Under the approximation each list l has a characteristic time T(l); the
 % level process of an item with request rate lam is a birth-death chain
 % with up-probability 1-e(l) and down-probability e(l), e(l)=exp(-lam*T(l)),
 % giving pi(l) proportional to prod_{s<=l} (1-e(s))/e(s). For h=1 this
 % reduces exactly to the Che approximation for LRU.
 %
 % @par Syntax:
 % @code
 % pij = cache_ttl_hlru(lambda, m)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lambda<td>(u x n x h+1) per-class per-item request rates (as built
 %                   by solver_mva_cache_analyzer; identical across lists)
 % <tr><td>m<td>(1 x h) list capacities
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>pij<td>(n x h+1) probabilities; column 1 = not cached, column
 %                1+l = in list l
 % </table>
%}
function pij = cache_ttl_hlru(lambda, m)
% aggregate the request rate of each item over the user classes; the
% analyzer replicates the same rate on every list slice, so read slice 2
lam = zeros(size(lambda,2),1);
for v = 1:size(lambda,1)
    lam = lam + reshape(lambda(v,:,min(2,size(lambda,3))), [], 1);
end
n = length(lam);
h = length(m);
m = m(:)';

% characteristic times solved by per-list bisection (Gauss-Seidel sweeps):
% sum_k pi_l(k; T) = m(l), with sum_k pi_l increasing in T(l)
T = ones(1,h) / max(mean(lam), GlobalConstants.FineTol);
maxSweeps = 200;
for sweep = 1:maxSweeps
    Told = T;
    for l = 1:h
        lo = 0;
        hi = max(T(l), 1/max(mean(lam), GlobalConstants.FineTol));
        % grow hi until the list-l occupancy reaches its capacity
        while cache_hlru_occ(lam, T, l, hi, h) < m(l) && hi < 1e12
            hi = 2*hi;
        end
        for it = 1:100
            mid = (lo+hi)/2;
            if cache_hlru_occ(lam, T, l, mid, h) < m(l)
                lo = mid;
            else
                hi = mid;
            end
        end
        T(l) = (lo+hi)/2;
    end
    if max(abs(T-Told)./max(Told,GlobalConstants.Zero)) < GlobalConstants.FineTol
        break
    end
end

pij = cache_hlru_levelprobs(lam, T, h); % (n x h+1): [pi_0, pi_1..pi_h]
end

function occ = cache_hlru_occ(lam, T, l, Tl, h)
% total occupancy of list l when its characteristic time is Tl
T(l) = Tl;
P = cache_hlru_levelprobs(lam, T, h);
occ = sum(P(:,1+l));
end

function P = cache_hlru_levelprobs(lam, T, h)
% birth-death level probabilities: pi_l ~ prod_{s<=l} (1-e_s)/e_s
n = length(lam);
P = zeros(n, h+1);
for k = 1:n
    w = zeros(1, h+1);
    w(1) = 1;
    for l = 1:h
        e = exp(-lam(k)*T(l));
        w(1+l) = w(l) * (1-e)/max(e, GlobalConstants.Zero);
    end
    P(k,:) = w / sum(w);
end
end
