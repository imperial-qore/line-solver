%{ @file retrieval_metrics.m
 %  @brief Performance metrics of a delayed-hit (list-based) cache
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes miss, hit and delayed-hit metrics for a delayed-hit cache
 %
 % @details
 % Evaluates the performance measures of Proposition "prop:performance_measures"
 % for a list-based cache with delayed hits whose retrieval system has one IS
 % station (s=0) and r PS stations (s=1,...,r). Letting E(m)=E(0,m) be the
 % normalizing constant and E_i the constant of the system without item i (both
 % computed with retrieval_nc), the metrics are
 %
 %   miss ratio          pi_{i,0} = E_i(m)/E(m)
 %   hit ratio (list j)  pi_{i,j} = m_j*gamma_{i,j}*E_i(m-1_j)/E(m)
 %   delayed hit (IS)    phi_{0,i} = lambda_i*eta_{0,i}*E_i(m)/E(m)
 %   delayed hit (PS s)  phi_{s,i} = lambda_i*eta_{s,i}*E_i(1_s,m)/E(m)
 %
 % These satisfy the balance pi_{i,0} + sum_s phi_{s,i} + sum_j pi_{i,j} = 1 for
 % every item i (eq. balance).
 %
 % @par Syntax:
 % @code
 % [pmiss, phit, pdh] = retrieval_metrics(m, lambda, eta, gamma)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>m<td>Cache list capacities (1 x h)
 % <tr><td>lambda<td>Per-item arrival rates (1 x n)
 % <tr><td>eta<td>Fetching demands eta(i,s+1)=eta_{s,i} (column 1 = IS station s=0, columns 2..r+1 = PS stations), (n x (r+1)). A PS column models any symmetric/identical-rate single-server discipline (PS, SIRO, FCFS, LCFS-PR), insensitive to the per-list service order.
 % <tr><td>gamma<td>Access factors gamma(i,j)=gamma_{i,j}, (n x h)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>pmiss<td>Miss ratios pi_{i,0}, (1 x n)
 % <tr><td>phit<td>Hit ratios pi_{i,j}, (h x n)
 % <tr><td>pdh<td>Delayed-hit probabilities phi_{s,i} for s=0,...,r, ((r+1) x n)
 % </table>
%}
function [pmiss,phit,pdh]=retrieval_metrics(m,lambda,eta,gamma)
n = numel(lambda);
m = m(:).';
h = numel(m);
r = size(eta,2)-1;
v0 = zeros(1,r);

E = retrieval_nc(v0,m,lambda,eta,gamma);

pmiss = zeros(1,n);
phit = zeros(h,n);
pdh = zeros(r+1,n);

for i=1:n
    keep = [1:i-1, i+1:n];
    lambda_i = lambda(keep);
    eta_i = eta(keep,:);
    gamma_i = gamma(keep,:);

    % miss ratio: pi_{i,0} = E_i(m)/E(m)
    Ei = retrieval_nc(v0,m,lambda_i,eta_i,gamma_i);
    pmiss(i) = Ei/E;

    % delayed hit at the IS station s=0: phi_{0,i} = lambda_i*eta_{0,i}*E_i(m)/E(m)
    pdh(1,i) = lambda(i)*eta(i,1)*Ei/E;

    % delayed hit at PS station s=1,...,r: phi_{s,i} = lambda_i*eta_{s,i}*E_i(1_s,m)/E(m)
    for s=1:r
        vs = v0;
        vs(s) = 1;
        Eis = retrieval_nc(vs,m,lambda_i,eta_i,gamma_i);
        pdh(s+1,i) = lambda(i)*eta(i,s+1)*Eis/E;
    end

    % hit ratio at list j=1,...,h: pi_{i,j} = m_j*gamma_{i,j}*E_i(m-1_j)/E(m)
    for j=1:h
        if m(j)>0
            Eij = retrieval_nc(v0,oner(m,j),lambda_i,eta_i,gamma_i);
            phit(j,i) = m(j)*gamma(i,j)*Eij/E;
        end
    end
end
end
