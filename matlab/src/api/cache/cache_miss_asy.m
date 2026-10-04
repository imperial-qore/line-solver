%{ @file cache_miss_asy.m
 %  @brief Asymptotic (large-cache) miss ratio by a rank-threshold fixed point
 %
 %  @author LINE Development Team
%}

%{
 % @brief Asymptotic miss ratio of a multi-list cache
 %
 % @details
 % The deterministic limit of a multi-list cache: as the item count grows,
 % list l holds exactly the m(l) items of largest effective popularity, so an
 % item's membership becomes a threshold test rather than a probability.
 % Writing pi(k) for the miss probability of item k, the effective popularity
 % of item j in list l is gamma(l,j)*(1-pi(j)) and the fixed point is
 %
 %   pi(k) = sum_l gamma(l,k) 1{item k outside the top m(l)} / sum_l gamma(l,k),
 %
 % iterated to a sup-norm tolerance from the uniform start pi = 1/n. The
 % returned scalar is the request-weighted miss ratio
 % sum_{l,k} gamma(l,k) pi(k) / sum_{l,k} gamma(l,k).
 %
 % INDEX CONVENTION, AND IT IS THE REVERSE OF EVERY OTHER CACHE FUNCTION HERE:
 % gamma is (h,n), LIST-major, whereas cache_spm, cache_erec and cache_miss all
 % take gamma as (n,h), ITEM-major. Callers holding an item-major gamma must
 % transpose. The threshold is strict, so an item exactly at the cutoff is
 % admitted. A degenerate capacity (zero total, or any negative entry) returns
 % 1, i.e. every request misses.
 %
 % Reference:
 %   N. Gast, B. Van Houdt, "Transient and steady-state regime of a family of
 %   list-based cache replacement algorithms", Queueing Syst. 83, 2016.
 %
 % @par Syntax:
 % @code
 % missratio = cache_miss_asy(gamma, m)
 % missratio = cache_miss_asy(gamma, m, maxiter, tol)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>(h,n) list-major access factors
 % <tr><td>m<td>(1,h) list capacities
 % <tr><td>maxiter<td>Cap on fixed-point sweeps (default 1000)
 % <tr><td>tol<td>Sup-norm stopping tolerance on pi (default 1e-8)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>missratio<td>Request-weighted asymptotic miss ratio
 % </table>
%}
function missratio = cache_miss_asy(gamma, m, maxiter, tol)
if nargin < 3 || isempty(maxiter)
    maxiter = 1000;
end
if nargin < 4 || isempty(tol)
    tol = 1e-8;
end

h = size(gamma,1);
n = size(gamma,2);

if sum(m(:)) == 0 || min(m(:)) < 0
    missratio = 1;
    return
end

pi_k = ones(1,n)/n;
for iter = 1:maxiter %#ok<NASGU>
    prev = pi_k;
    newpi = zeros(1,n);
    for k = 1:n
        numer = 0;
        denom = 0;
        for l = 1:h
            cap = floor(m(l));
            if cap <= 0
                continue
            end
            other = [1:(k-1), (k+1):n];
            pop = gamma(l,other) .* (1 - prev(other));
            pop = sort(pop, 'descend');
            take = min(cap, numel(pop));
            if take < cap
                % fewer competitors than slots: item k is always cached
                notin = 0;
            else
                notin = 1;
                if gamma(l,k) * (1 - prev(k)) > pop(take)
                    notin = 0;
                end
            end
            numer = numer + gamma(l,k) * notin;
            denom = denom + gamma(l,k);
        end
        if denom > GlobalConstants.Zero
            newpi(k) = numer / denom;
        else
            newpi(k) = 1;
        end
    end
    pi_k = newpi;
    if max(abs(pi_k - prev)) < tol
        break
    end
end

missrate = sum(sum(gamma .* repmat(pi_k, h, 1)));
totrate = sum(gamma(:));
if totrate > GlobalConstants.Zero
    missratio = missrate / totrate;
else
    missratio = 1;
end
end
