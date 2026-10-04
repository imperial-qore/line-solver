%{ @file cache_cost.m
 %  @brief Mean per-list storage cost of a cache with item sizes
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes the mean storage cost held by each cache list
 %
 % @details
 % Evaluates K_j = sum_i sigma_i pi_ij, the expected storage cost of the
 % items resident in list j at steady state, as defined in Casale-Gast
 % (IEEE/ACM ToN, 2021), Sec. IX. With no cost caps supplied the same
 % expression still applies and reports the mean cost of the unconstrained
 % model.
 %
 % @par Syntax:
 % @code
 % [K,pij] = cache_cost(gamma, m, sigma)
 % [K,pij] = cache_cost(gamma, m, sigma, k)
 % [K,pij] = cache_cost(gamma, m, sigma, k, pij)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>Item popularity probabilities (n x h)
 % <tr><td>sigma<td>Item storage costs (sizes), 1 x n
 % <tr><td>m<td>Cache capacity vector, 1 x h
 % <tr><td>k<td>(Optional) per-list storage cost caps, 1 x h
 % <tr><td>pij<td>(Optional) precomputed occupancy matrix, n x (h+1) with column 1 the miss probability
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Mean storage cost of each list, 1 x h
 % <tr><td>pij<td>Occupancy matrix used, n x (h+1)
 % </table>
 %
 % @see cache_prob_erec
%}
function [K,pij] = cache_cost(gamma,m,sigma,k,pij)
[n,h] = size(gamma);
sigma = sigma(:).';
if numel(sigma)~=n
    line_error(mfilename,'The item size vector must have one entry per item.');
end
if nargin<4
    k = [];
else
    k = k(:).';
end
if nargin<5 || isempty(pij)
    pij = cache_prob_erec(gamma,m,sigma,k);
end
K = zeros(1,h);
for j=1:h
    K(j) = sigma * pij(:,1+j);
end
end
