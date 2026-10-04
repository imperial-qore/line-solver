%{ @file retrieval_mva.m
 %  @brief Exact MVA-style recursion for delayed-hit (list-based) cache metrics
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes exact miss, hit and delayed-hit metrics via the MVA recursion
 %
 % @details
 % Implements the exact recursive characterization of the paper (Theorems
 % "thm:arvthm" and "thm:pi_xi", and eq. motivation_to_pi0_FPI) that the FPI
 % heuristic (retrieval_fpi) approximates. For a list-based cache with delayed
 % hits whose retrieval system has one IS station (s=0) and r PS stations
 % (s=1,...,r), with phi^{(k)}_{s,i} the delayed-hit probability of item i at
 % station s in the system WITHOUT item k:
 %
 %   theta_{ij}(m) = gamma_{ij} / (1 + lambda_i eta_{0,i}
 %                     + sum_{s=1}^r lambda_i eta_{s,i}(1 + sum_{k!=i} phi^{(i)}_{s,k}(m-1_j)))
 %   xi_j(m)       = m_j / sum_i theta_{ij}(m)(1 - pihit_i(m-1_j))
 %   pi_{ij}(m)    = theta_{ij}(m) xi_j(m) (1 - pihit_i(m-1_j))
 %   pihit_i(m)    = sum_j pi_{ij}(m)
 %   pi_{i0}(m)    = (1 - pihit_i(m)) / (1 + lambda_i eta_{0,i}
 %                     + sum_{s=1}^r lambda_i eta_{s,i}(1 + sum_{k!=i} phi^{(i)}_{s,k}(m)))
 %   phi_{sk}(m)   = lambda_k pi_{k0}(m) eta_{s,k}(1 + sum_{i!=k} phi^{(k)}_{s,i}(m))   (s=1..r)
 %   phi_{0k}(m)   = lambda_k eta_{0,k} pi_{k0}(m)
 %
 % The recursion terminates at the empty item set and at the empty cache
 % (pihit_i = 0). This is exact and agrees with retrieval_nc / retrieval_metrics;
 % it is memoized over (item-subset, capacity) with O(2^n n^2 h r prod_j(1+m_j))
 % time, so it is feasible only for small systems (use retrieval_fpi otherwise).
 %
 % @par Syntax:
 % @code
 % [pmiss, phit, pdh] = retrieval_mva(m, lambda, eta, gamma)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>m<td>Cache list capacities (1 x h)
 % <tr><td>lambda<td>Per-item arrival rates (1 x n)
 % <tr><td>eta<td>Fetching demands eta(i,s+1)=eta_{s,i} (column 1 = IS station s=0, columns 2..r+1 = PS stations), (n x (r+1))
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
function [pmiss,phit,pdh]=retrieval_mva(m,lambda,eta,gamma)
n = numel(lambda);
m = m(:).';
h = numel(m);
r = size(eta,2) - 1;
lambda = lambda(:).';        % 1 x n
eta0 = eta(:,1).';           % 1 x n
etaPS = eta(:,2:end);        % n x r

radix = m + 1;
ncap = prod(radix);
nmask = 2^n;

done  = false(nmask, ncap);
PI0   = zeros(nmask, ncap, n);
PIHIT = zeros(nmask, ncap, n);
PIJ   = zeros(nmask, ncap, n, h);
PHI   = zeros(nmask, ncap, n, r+1);   % station index s=0..r -> page 1..r+1

fullmask = nmask - 1;
solve(fullmask, m);

ci = capidx(m);
pmiss = zeros(1,n);
phit = zeros(h,n);
pdh = zeros(r+1,n);
for i = 1:n
    pmiss(i) = PI0(fullmask+1, ci, i);
    for j = 1:h
        phit(j,i) = PIJ(fullmask+1, ci, i, j);
    end
    for s = 1:r+1
        pdh(s,i) = PHI(fullmask+1, ci, i, s);
    end
end

    function idx = capidx(c)
        idx = 1; mul = 1;
        for jj = 1:h
            idx = idx + c(jj)*mul;
            mul = mul*radix(jj);
        end
    end

    function solve(mask, c)
        ci_ = capidx(c);
        if done(mask+1, ci_), return; end
        if mask == 0
            done(mask+1, ci_) = true;
            return
        end
        active = find(bitget(mask, 1:n));

        % boundary: when the cache can hold every active item (sum(c) >= |active|),
        % all items are permanently cached -> pihit=1, pi0=0, phi=0 (no fetching).
        if sum(c) >= numel(active)
            for i = active
                PIHIT(mask+1, ci_, i) = 1;
            end
            done(mask+1, ci_) = true;
            return
        end

        % --- ensure dependencies are solved ---
        for jj = 1:h
            if c(jj) > 0
                cj = c; cj(jj) = cj(jj) - 1;
                solve(mask, cj);
                for i = active
                    solve(bitset(mask, i, 0), cj);
                end
            end
        end
        for k = active
            solve(bitset(mask, k, 0), c);
        end

        % --- theta, xi, pi_{ij} (cache recursion on m-1_j) ---
        for jj = 1:h
            if c(jj) > 0
                cj = c; cj(jj) = cj(jj) - 1; cjx = capidx(cj);
                theta = zeros(1,n);
                for i = active
                    maski = bitset(mask, i, 0);
                    acc = 0;
                    for s = 1:r
                        sphi = 0;
                        for k = active
                            if k ~= i
                                sphi = sphi + PHI(maski+1, cjx, k, s+1);
                            end
                        end
                        acc = acc + lambda(i)*etaPS(i,s)*(1 + sphi);
                    end
                    theta(i) = gamma(i,jj) / (1 + lambda(i)*eta0(i) + acc);
                end
                sden = 0;
                for i = active
                    sden = sden + theta(i)*(1 - PIHIT(mask+1, cjx, i));
                end
                xi_j = c(jj) / sden;
                for i = active
                    PIJ(mask+1, ci_, i, jj) = theta(i)*xi_j*(1 - PIHIT(mask+1, cjx, i));
                end
            end
        end

        % --- pihit_i = sum_j pi_{ij} ---
        for i = active
            ph = 0;
            for jj = 1:h
                ph = ph + PIJ(mask+1, ci_, i, jj);
            end
            PIHIT(mask+1, ci_, i) = ph;
        end

        % --- pi_{i0} (uses phi^{(i)}(m) at same capacity) ---
        for i = active
            maski = bitset(mask, i, 0);
            acc = 0;
            for s = 1:r
                sphi = 0;
                for k = active
                    if k ~= i
                        sphi = sphi + PHI(maski+1, ci_, k, s+1);
                    end
                end
                acc = acc + lambda(i)*etaPS(i,s)*(1 + sphi);
            end
            PI0(mask+1, ci_, i) = (1 - PIHIT(mask+1, ci_, i)) / (1 + lambda(i)*eta0(i) + acc);
        end

        % --- phi_{sk}(m) ---
        for k = active
            maskk = bitset(mask, k, 0);
            pi0k = PI0(mask+1, ci_, k);
            PHI(mask+1, ci_, k, 1) = lambda(k)*eta0(k)*pi0k;   % s = 0 (IS)
            for s = 1:r
                sphi = 0;
                for i = active
                    if i ~= k
                        sphi = sphi + PHI(maskk+1, ci_, i, s+1);
                    end
                end
                PHI(mask+1, ci_, k, s+1) = lambda(k)*pi0k*etaPS(k,s)*(1 + sphi);
            end
        end

        done(mask+1, ci_) = true;
    end
end
