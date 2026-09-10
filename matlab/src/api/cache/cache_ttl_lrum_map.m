%{ @file cache_ttl_lrum_map.m
 %  @brief TTL approximation of LRU(m) caches under MAP request streams
 %
 %  @author LINE Development Team
%}

%{
 % @brief Request-weighted hit/miss probabilities for LRU(m) with MAP requests
 %
 % @details
 % Front-end of the Gast-Van Houdt (Performance Evaluation 2017) TTL
 % approximation for LRU(m) caches whose items have Markovian arrival
 % process request streams. Intended for items with genuinely distinct or
 % correlated request processes (e.g. marked MAP arrivals); when items are
 % i.i.d. marks of a common stream the request sequence is IRM and the
 % Poisson-based TTL approximations already apply.
 %
 % @par Syntax:
 % @code
 % [pij, pijtime] = cache_ttl_lrum_map(D0c, D1c, m)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>D0c<td>Cell (1,n); D0c{k} is the (d,d) hidden-transition matrix of item k
 % <tr><td>D1c<td>Cell (1,n); D1c{k} is the (d,d) arrival matrix of item k
 % <tr><td>m<td>Cache capacity vector (1,h)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>pij<td>(n,h+1) request-weighted probabilities; column 1 is the
 %                probability that a request for item k misses, column 1+l
 %                the probability that it hits in list l
 % <tr><td>pijtime<td>(n,h+1) time-stationary level occupancy probabilities
 % </table>
%}
function [pij, pijtime] = cache_ttl_lrum_map(D0c, D1c, m)
n = numel(D0c);
h = numel(m);

t = cache_t_lrum_map(D0c, D1c, m);

pij = zeros(n, h+1);
pijtime = zeros(n, h+1);
for k = 1:n
    [prob, ~, hitfrac] = cache_lrum_map_levelstats(D0c{k}, D1c{k}, t);
    pijtime(k, :) = prob;
    pij(k, 2:end) = hitfrac;
    pij(k, 1) = max(0, 1 - sum(hitfrac));
end
end
