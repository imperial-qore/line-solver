%{ @file cache_t_hlru.m
 %  @brief Characteristic times for h-LRU / LRU(m) cache lists
 %
 %  @author LINE Development Team
%}

%{
 % @brief Characteristic time of each list of an h-LRU cache
 %
 % @details
 % Solves the TTL (characteristic-time) fixed point of the list-based h-LRU
 % (LRU(m)) policy: sum_k pi_l(k;T) = m(l) for each list l, where the level
 % probabilities follow the birth-death form pi_l ~ prod_{s<=l} (1-e_s)/e_s
 % with e_s = exp(-gamma_k*T(s)) (Gast and Van Houdt, SIGMETRICS 2015).
 % Solved by per-list bisection with Gauss-Seidel sweeps; no Optimization
 % Toolbox dependency.
 %
 % @par Syntax:
 % @code
 % t = cache_t_hlru(gamma, m)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>(n x 1) per-item request rates; an (n x h) matrix is
 %                  accepted for backward compatibility (first column used)
 % <tr><td>m<td>(1 x h) list capacities
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>t<td>(1 x h) characteristic time of each list
 % </table>
%}
function t = cache_t_hlru(gamma, m)
lam = gamma(:,1);
n = length(lam); %#ok<NASGU>
h = length(m);
m = m(:)';

t = ones(1,h) / max(mean(lam), GlobalConstants.FineTol);
maxSweeps = 200;
for sweep = 1:maxSweeps
    told = t;
    for l = 1:h
        lo = 0;
        hi = max(t(l), 1/max(mean(lam), GlobalConstants.FineTol));
        while occ_l(lam, t, l, hi, h) < m(l) && hi < 1e12
            hi = 2*hi;
        end
        for it = 1:100
            mid = (lo+hi)/2;
            if occ_l(lam, t, l, mid, h) < m(l)
                lo = mid;
            else
                hi = mid;
            end
        end
        t(l) = (lo+hi)/2;
    end
    if max(abs(t-told)./max(told,GlobalConstants.Zero)) < GlobalConstants.FineTol
        break
    end
end
end

function occ = occ_l(lam, t, l, tl, h)
t(l) = tl;
n = length(lam);
occ = 0;
for k = 1:n
    w = zeros(1, h+1);
    w(1) = 1;
    for s = 1:h
        e = exp(-lam(k)*t(s));
        w(1+s) = w(s) * (1-e)/max(e, GlobalConstants.Zero);
    end
    occ = occ + w(1+l)/sum(w);
end
end
