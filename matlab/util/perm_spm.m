function [P, lP, xi, lcap] = perm_spm(A, m)
% PERM_SPM  Saddle-point approximation of the permanent of a positive matrix.
%
% [P, LP, XI, LCAP] = PERM_SPM(A, M) approximates the permanent of the matrix
% obtained from the N-by-H matrix A by repeating column L exactly M(L) times.
% SUM(M) must equal N. With M omitted, A must be square and PERM_SPM
% approximates PERM(A).
%
% THE HOMOGENEOUS VARIANT OF CACHE_SPM. Both routines evaluate the same Cauchy
% integral by Laplace's method and differ only in the generating function whose
% coefficient they extract:
%
%   cache_spm   E(m) = prod_l m_l! [prod_l z_l^m_l] prod_k (1 + sum_l g_kl z_l)
%   perm_spm    P    = prod_l m_l! [prod_l z_l^m_l] prod_k (    sum_l A_kl z_l)
%
% The cache factor carries a "+1" because an item may stay out of the cache, so
% the coefficient it extracts is a RECTANGULAR permanent over N items and
% SUM(M) < N slots. Dropping the "+1" forces every row to be matched, which is
% exactly the permanent and requires SUM(M) = N. That is the one case CACHE_SPM
% cannot serve: at N = SUM(M) its multipliers diverge, and it falls back on
% CACHE_EREC rather than expanding. Here the integrand is homogeneous instead,
% and the saddle point is interior in the H-1 directions that survive.
%
% METHOD. Write z_l = xi_l exp(i th_l). The saddle point in xi solves
%
%   sum_k A_kl xi_l / (sum_j A_kj xi_j) = m_l,   l = 1..H,
%
% that is, PP(k,l) = A(k,l) xi(l) / S(k) with S = A*xi is the diagonal scaling
% of A to row sums 1 and column sums M (Sinkhorn scaling; doubly stochastic
% when M is all ones). There
%
%   PHI = sum_k log S(k) - sum_l m_l log xi_l
%
% is the log of the Gurvits capacity, returned as LCAP, and PERM <= EXP(LCAP)
% holds for every nonnegative matrix. The Gaussian correction uses
%
%   H = diag(M) - PP'*PP,
%
% a weighted graph Laplacian on the columns: H*ONES = 0, which is the
% invariance of the integrand under th -> th + c*ONES that homogeneity creates.
% That direction is a full period rather than a Gaussian, so it contributes
% 2*pi and leaves an (H-1)-dimensional Laplace integral. Any principal
% (H-1)-by-(H-1) submatrix serves, because all cofactors of a Laplacian are
% equal, and the estimate is
%
%   LP = sum_l log(m_l!) - (H-1)/2 log(2 pi) + PHI - 1/2 log det(H_red).
%
% ACCURACY, AND WHAT IT IS NOT. The estimate is EXACT for H = 1, where the
% permanent is N! prod_k A(k,1). It is a genuine asymptotic expansion as MIN(M)
% grows with H held fixed: measured against the exact permanent on random
% positive matrices the ratio falls from 1.11 at M = (2,2,2) to 1.02 at
% M = (3,3). At M = ONES the dimension of the integral grows with the expansion
% parameter and the leading term keeps a systematic BIAS. For the N-by-N matrix
% of ones it returns (2 pi)^(-(N-1)/2) N^(N+1/2) against the exact N!, a ratio
% tending to (e/sqrt(2 pi))^N = 1.084^N, and random positive matrices track that
% closely (1.31 at N = 4, 1.87 at N = 8). So at M = ONES this OVERESTIMATES,
% with a spread across matrices far tighter than the bias itself; it is not a
% bound in either direction. For scale, the raw capacity EXP(LCAP) is off by
% 352x on the same N = 8 instances, and PERM_BETHE is a genuine lower bound.
%
% Input:
%   A - N-by-H strictly positive matrix. A zero entry is refused by
%       PERM_REQUIRE_SUPPORT rather than floored, as in PERM_HEUR
%   M - 1-by-H non-negative integer column multiplicities with SUM(M) = N
%       (default: ONES(1,N), which requires A square)
%
% Output:
%   P    - approximate permanent
%   LP   - its logarithm, correct even when P overflows
%   XI   - 1-by-H saddle point, scaled to unit geometric mean, zero on a
%          column of multiplicity zero
%   LCAP - log of the Gurvits capacity at the saddle point, an upper bound on
%          the log permanent
%
% Example:
%   A = rand(6) + 0.05;
%   [P, lP] = perm_spm(A);       % approximates perm(A)
%   perm_spm(rand(6,2)+0.05, [3 3]);  % perm of that matrix with both columns tripled
%
% References:
%   G. Casale, "Accelerating Performance Inference over Closed Systems by
%   Asymptotic Methods", ACM SIGMETRICS, 2017 (the saddle-point expansion)
%   L. Gurvits, "Hyperbolic Polynomials Approach to Van der Waerden /
%   Schrijver-Valiant Like Conjectures", STOC, 2006 (the capacity bound)
%   R. Sinkhorn, "A Relationship Between Arbitrary Positive Matrices and Doubly
%   Stochastic Matrices", Ann. Math. Statist. 35(2), 1964 (the scaling)
%
% See also PERM, PERM_HEUR, PERM_BETHE, CACHE_SPM.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 1 || isempty(A)
    P = 1; lP = 0; xi = []; lcap = 0;   % the permanent of the empty matrix is 1
    return
end
if any(A(:) < 0)
    error('perm_spm:InvalidInput', 'Matrix must be non-negative.');
end

[n, h] = size(A);
if nargin < 2 || isempty(m)
    if h ~= n
        error('perm_spm:InvalidInput', ['Without column multiplicities the matrix must be ' ...
            'square; A is %d-by-%d.'], n, h);
    end
    m = ones(1, n);
end
m = m(:)';
if numel(m) ~= h
    error('perm_spm:InvalidInput', ['The multiplicity vector has %d entries against %d ' ...
        'columns of A.'], numel(m), h);
end
if any(m < 0) || any(abs(m - round(m)) > 0)
    error('perm_spm:InvalidInput', 'Column multiplicities must be non-negative integers.');
end
if sum(m) ~= n
    error('perm_spm:InvalidInput', ['The column multiplicities must sum to the number of ' ...
        'rows (sum(m) = %g against %d rows). The integrand is homogeneous of degree %d, so ' ...
        'every other coefficient of it is exactly zero.'], sum(m), n, n);
end
perm_require_support(A, mfilename);

% A column repeated zero times leaves the permanent unchanged, and its xi is a
% boundary of the Laplace integral rather than a direction of it, so it must
% leave the expansion. Dropping it is exact: setting z_l = 0 in the generating
% function removes the column, and prod_l m_l! is unchanged because 0! = 1.
keep = find(m > 0);
Ak = A(:, keep);
mk = m(keep);
hk = numel(keep);

xi = zeros(1, h);

% Saddle point: scale A to row sums 1 and column sums mk. Row sums are 1 by
% construction of PP, so only the column sums are iterated on.
tol = 1e-11;
maxiter = 10000;
xik = ones(hk, 1);
converged = false;
margin = Inf;
for it = 1:maxiter %#ok<NASGU>
    S = Ak * xik;
    csum = xik .* (Ak' * (1 ./ S));
    margin = max(abs(csum - mk'));
    if margin < tol
        converged = true;
        break
    end
    xik = xik .* (mk' ./ csum);
    xik = xik / exp(mean(log(xik)));    % the saddle is a ray; pin its scale
end
if ~converged
    line_error(mfilename, ['The scaling to row sums 1 and column sums m did not converge in ' ...
        '%d sweeps (margin error %g against a tolerance of %g). The expansion assumes the ' ...
        'saddle point, so no value is returned. The usual cause is a matrix without total ' ...
        'support.'], maxiter, margin, tol);
end
xi(keep) = xik;

S = Ak * xik;
PP = (Ak .* xik') ./ S;
lcap = sum(log(S)) - mk * log(xik);

% Reduced Hessian. H = diag(mk) - PP'*PP is a Laplacian, so it is singular
% along ONES and all of its principal cofactors are equal; the last index is
% dropped only because one has to be. Strict positivity of A makes the column
% graph complete, hence H_red positive definite and CHOL the right factor.
if hk > 1
    H = diag(mk) - PP' * PP;
    [R, p] = chol(H(1:(hk-1), 1:(hk-1)));
    if p > 0
        line_error(mfilename, ['The reduced Hessian is not positive definite (leading minor ' ...
            '%d), so the saddle point is degenerate and the Gaussian factor does not exist.'], p);
    end
    ldet = 2 * sum(log(diag(R)));
else
    ldet = 0;   % no direction survives the homogeneity, and det of the empty matrix is 1
end

lP = sum(gammaln(mk + 1)) - 0.5 * (hk - 1) * log(2 * pi) + lcap - 0.5 * ldet;
P = exp(lP);
end
