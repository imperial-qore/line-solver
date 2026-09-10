%{ @file retrieval_nc.m
 %  @brief Exact recurrence for the normalizing constant of a delayed-hit (list-based) cache
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes the delayed-hit cache normalizing constant E(v,m) by the exact recurrence
 %
 % @details
 % This function recursively computes the normalizing constant E(v,m) of a
 % list-based cache with delayed hits, whose retrieval system has one IS station
 % (index s=0) and r PS stations (s=1,...,r), using the exact recurrence relation
 %
 %   E(v,m) = (1 + lambda_k*eta_{0,k}) E_k(v,m)
 %          + sum_{s=1}^{r} lambda_k*eta_{s,k}*(v_s+1) E_k(v+1_s, m)
 %          + sum_{j=1}^{h} m_j*gamma_{k,j} E_k(v, m-1_j)
 %
 % where E_k is the same constant for the system without item k. Boundary
 % conditions are E = 1 if there are no items left and E = 0 if sum_j m_j exceeds
 % the number of items or any m_j < 0. The plain normalizing constant is recovered
 % as E(m) = E(0,m).
 %
 % @par Syntax:
 % @code
 % E = retrieval_nc(v, m, lambda, eta, gamma)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>v<td>Moment-order vector for the PS stations (zeros(1,r) for the plain constant)
 % <tr><td>m<td>Cache list capacities
 % <tr><td>lambda<td>Per-item arrival rates
 % <tr><td>eta<td>Fetching demands eta(i,s+1)=eta_{s,i} (column 1 = IS station s=0, columns 2..r+1 = PS stations). A PS column models any symmetric/identical-rate single-server discipline (PS, SIRO, FCFS, LCFS-PR), insensitive to the per-list service order.
 % <tr><td>gamma<td>Access factors gamma(i,j)=gamma_{i,j}
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>E<td>Normalizing constant E(v,m)
 % </table>
%}
function E=retrieval_nc(v,m,lambda,eta,gamma)
E = sub_retrieval_nc(v(:).',m(:).',lambda,eta,gamma,numel(lambda));
end

function E=sub_retrieval_nc(v,m,lambda,eta,gamma,k)
r=numel(v);
h=numel(m);
if sum(m)>k || min(m)<0
    E=0;
    return
end
if k==0
    E=1;
    return
end

% item k either out of the system or fetched at the IS station s=0
E = (1+lambda(k)*eta(k,1))*sub_retrieval_nc(v,m,lambda,eta,gamma,k-1);

% item k fetched at PS station s=1,...,r
for s=1:r
    vp=v;
    vp(s)=vp(s)+1;
    E = E + lambda(k)*eta(k,s+1)*(v(s)+1)*sub_retrieval_nc(vp,m,lambda,eta,gamma,k-1);
end

% item k stored in cache list j=1,...,h
for j=1:h
    if m(j)>0
        E = E + gamma(k,j)*m(j)*sub_retrieval_nc(v,oner(m,j),lambda,eta,gamma,k-1);
    end
end
end
