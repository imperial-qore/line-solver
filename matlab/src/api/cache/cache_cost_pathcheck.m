%{ @file cache_cost_pathcheck.m
 %  @brief Screens a cost-capped cache for lists that no item can ever reach
 %
 %  @author LINE Development Team
%}

%{
 % @brief Detects promotion paths blocked by storage cost caps
 %
 % @details
 % The constrained normalizing constant E(m,k) of Casale-Gast (IEEE/ACM ToN,
 % 2021), Sec. IX, sums the product form over every size-feasible cache
 % state. Under RR-C(m) an item only reaches list j by being promoted along
 % the path from the miss list to j, one list at a time, so a cap on an
 % intermediate list can make size-feasible states unreachable. When that
 % happens the size-feasible set is no longer a single recurrent class and
 % E(m,k) normalizes over states the cache never visits.
 %
 % This routine reports each (item, list) pair that is size-feasible for the
 % list but blocked on the way to it. An empty report is a necessary, not
 % sufficient, condition: a list of capacity above one may still be
 % unreachable when its cap admits no combination containing the item.
 %
 % @par Syntax:
 % @code
 % viol = cache_cost_pathcheck(gamma, sigma, k, parent)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>Item access factors (n x h)
 % <tr><td>sigma<td>Item storage costs (sizes), 1 x n
 % <tr><td>k<td>Per-list storage cost caps, 1 x h
 % <tr><td>parent<td>Parent list of each list, 1 x h, 0 for lists rooted in the miss list
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>viol<td>Matrix of blocked pairs, one row [item, list, blockinglist]
 % </table>
 %
 % @see cache_erec
%}
function viol = cache_cost_pathcheck(gamma,sigma,k,parent)
[n,h] = size(gamma);
sigma = sigma(:).';
k = k(:).';
parent = parent(:).';
viol = zeros(0,3);
for i=1:n
    for j=1:h
        if gamma(i,j)==0 || sigma(i)>k(j)
            continue % item i never resides in list j anyway
        end
        l = parent(j);
        while l>0
            if sigma(i)>k(l)
                viol(end+1,:) = [i, j, l]; %#ok<AGROW>
                break
            end
            l = parent(l);
        end
    end
end
end
