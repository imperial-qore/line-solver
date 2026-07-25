%{ @file retrieval_fpi.m
 %  @brief Fixed-point iteration (FPI) heuristic for delayed-hit (list-based) caches
 %
 %  @author LINE Development Team
%}

%{
 % @brief Approximates miss, hit and delayed-hit metrics of a delayed-hit cache via FPI
 %
 % @details
 % Implements the fixed-point heuristic of the paper (Sec. "Fixed-Point
 % Heuristic", eqs. theta/xi/pi/pi_0/phi_FPI) for a list-based cache with delayed
 % hits whose retrieval system has one IS station (s=0) and r PS stations
 % (s=1,...,r). The heuristic truncates the regular perturbation expansion at 0th
 % order (pi_{i,l}(m-1_j) ~ pi_{i,l}(m)) and solves the resulting nonlinear system
 % by successive substitution. Iteration t -> t+1:
 %
 %   D_i        = 1 + lambda_i*eta_{0,i} + sum_{s=1}^r lambda_i*eta_{s,i}*(1 + sum_{k!=i} phi_{s,k})
 %   theta_{ij} = gamma_{ij} / D_i
 %   xi_j       = m_j / sum_k theta_{kj} (1 - sum_l pi_{kl})
 %   pi_{ij}    = theta_{ij} xi_j / (1 + sum_l theta_{il} xi_l)
 %   pi_{i0}    = (1 - sum_j pi_{ij}) / D_i
 %   phi_{si}   = lambda_i*eta_{s,i}*(1 + sum_{k!=i} phi_{s,k}) pi_{i0}   (s=1..r)
 %   phi_{0i}   = lambda_i*eta_{0,i} pi_{i0}
 %
 % where phi_{s,k} is the delayed-hit probability of item k at PS station s. (The
 % paper writes phi_{r,k} in eqs. theta_FPI/phi_FPI, but the exact recursion it
 % approximates, eq. theta RS, uses the per-station index s; the two coincide for
 % a single PS station.) The solution satisfies the balance
 % pi_{i0} + sum_s phi_{si} + sum_j pi_{ij} = 1. Output orientation matches
 % retrieval_metrics for direct comparison.
 %
 % @par Syntax:
 % @code
 % [pmiss, phit, pdh] = retrieval_fpi(m, lambda, eta, gamma)
 % [pmiss, phit, pdh] = retrieval_fpi(m, lambda, eta, gamma, max_iter, tol)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>m<td>Cache list capacities (1 x h)
 % <tr><td>lambda<td>Per-item arrival rates (1 x n)
 % <tr><td>eta<td>Fetching demands eta(i,s+1)=eta_{s,i} (column 1 = IS station s=0, columns 2..r+1 = PS stations), (n x (r+1)). A PS column models any symmetric/identical-rate single-server discipline (PS, SIRO, FCFS, LCFS-PR), all treated alike via the mean-field sharing slowdown.
 % <tr><td>gamma<td>Access factors gamma(i,j)=gamma_{i,j}, (n x h)
 % <tr><td>max_iter<td>Maximum number of iterations (default 1e5)
 % <tr><td>tol<td>Convergence tolerance on relative change (default 1e-6)
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
function [pmiss,phit,pdh]=retrieval_fpi(m,lambda,eta,gamma,max_iter,tol)
if nargin < 5 || isempty(max_iter), max_iter = 1000; end
if nargin < 6 || isempty(tol), tol = 1e-6; end

n = numel(lambda);
m = m(:).';
h = numel(m);
r = size(eta,2) - 1;
lambda = lambda(:);          % n x 1
eta0 = eta(:,1);             % n x 1 (IS station s=0)
etaPS = eta(:,2:end);        % n x r (PS stations s=1..r)

% initial guess (paper Sec. Fixed-Point Heuristic)
phi = ones(r+1,n) / ((h+1)*(r+2));   % phi(s+1,i), s=0..r
pij = ones(h,n) / (h+1);             % pij(j,i)
pi0 = ones(1,n) / ((h+1)*(r+2));     % pi0(i)

for t = 1:max_iter
    % per-station PS-sharing factor F(s,i) = 1 + sum_{k!=i} phi_{s,k}
    % (the paper writes phi_{r,k}, but the exact recursion eq. theta RS uses the
    %  station index s; for a single PS station the two coincide.)
    F = ones(r, n);
    for s = 1:r
        phis = phi(s+1,:);
        F(s,:) = 1 + (sum(phis) - phis);
    end

    D = 1 + lambda.*eta0;                  % n x 1
    for s = 1:r
        D = D + lambda .* etaPS(:,s) .* F(s,:).';
    end

    theta = gamma ./ D;                    % n x h

    oneminus = (1 - sum(pij,1)).';         % n x 1 (uses pij^[t])
    xi = zeros(1,h);
    for j = 1:h
        xi(j) = m(j) / sum(theta(:,j).*oneminus);
    end

    txi = theta .* xi;                     % n x h
    pij_new = (txi ./ (1 + sum(txi,2))).'; % h x n

    pi0_new = (1 - sum(pij_new,1)) ./ D.'; % 1 x n

    phi_new = zeros(r+1,n);
    phi_new(1,:) = (lambda.*eta0).' .* pi0_new;        % phi_{0,i}
    for s = 1:r
        phi_new(s+1,:) = (lambda.*etaPS(:,s)).' .* F(s,:) .* pi0_new;  % phi_{s,i}
    end

    % convergence: max relative change in miss, hit, delayed-hit ratios
    delta = max([reldiff(pi0_new,pi0), reldiff(pij_new,pij), reldiff(phi_new,phi)]);

    pij = pij_new;
    pi0 = pi0_new;
    phi = phi_new;

    if ~isfinite(delta)
        line_warning(mfilename,'FPI diverged (non-finite iterate) at iteration %d.\n', t);
        break
    end
    if delta < tol
        break
    end
    if t == max_iter
        line_warning(mfilename,'FPI did not converge within %d iterations (last rel. change %.2e).\n', max_iter, delta);
    end
end

pmiss = pi0;
phit = pij;
pdh = phi;
end

function d = reldiff(a,b)
denom = max(abs(b(:)));
if denom == 0
    denom = 1;
end
d = max(abs(a(:)-b(:))) / denom;
end
