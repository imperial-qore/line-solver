function [x,flag,relres,iter] = gmres_iterate(A,L,U,b,x0,tol,restart,maxit)
% [X,FLAG,RELRES,ITER]=GMRES_ITERATE(A,L,U,B,X0,TOL,RESTART,MAXIT)
%
% Right-preconditioned restarted GMRES(m) (Saad and Schultz, 1986) on an ALREADY
% equilibrated, reordered and factorized system. Internal to CTMC_GMRES and
% CTMC_GMRES_MULTI, which share it so that one factorization can serve a whole
% block of right-hand sides.
%
% THE BUILT-IN GMRES IS NOT USED, for the same reason BICGSTAB_ITERATE replaces
% the built-in BICGSTAB. MATLAB applies an (L,U) pair on the LEFT, so it expands
% the Krylov space of M\A and minimizes the PRECONDITIONED residual, while the
% Java, Python and C++ ports of this kernel expand the space of A*M^-1 and
% minimize the TRUE one. Left preconditioning is also what forced CTMC_GMRES to
% keep a PIVOTED incomplete factorization: with the lean non-pivoting Crout
% factor the built-in stagnates (flag 3, true deviation 1.6e-04), so the kernel
% had to pay for an 'ilutp' factor carrying 4,981,518 nonzeros against the Crout
% factor's 30,015 on the M/M/1/K generator at K = 5000 -- a factor of 166, and
% 20%% of a dense 5001x5001, on exactly the large models the Krylov path exists
% for. Written out here and applied on the right, the same recurrence converges
% on the Crout factor, so both MATLAB kernels now use one lean factorization and
% run the algorithm the other three codebases run.
%
% The preconditioner is applied as U\(L\v); L and U are the incomplete factors,
% or the Jacobi pair (diag(1./d), I) when the factorization broke down.
%
% Orthogonalization is classical Gram-Schmidt run TWICE. Once is unstable on a
% generator whose rows span several orders of magnitude, and twice is enough
% (Giraud, Langou and Rozloznik, 2005); modified Gram-Schmidt would cost the same
% products without vectorizing.
%
% FLAG follows the MATLAB GMRES convention: 0 converged, 1 iteration limit
% reached, 2 preconditioner ill-conditioned, 3 stagnation. ITER counts
% matrix-vector products with A, i.e. total inner iterations across all restart
% cycles, which is what the other three codebases report.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Relative threshold below which the Arnoldi vector is treated as a lucky
% breakdown: the Krylov space is already invariant and cannot be extended.
BREAKDOWN_TOL = 1e-14;

n = size(A,1);
x = reshape(x0,n,1);
b = reshape(b,n,1);

bnorm = norm(b);
if bnorm == 0
    bnorm = 1;
end

r = b - A*x;
beta = norm(r);
relres = beta/bnorm;
iter = 0;
if relres <= tol
    flag = 0;
    return
end

V = zeros(n,restart+1);
H = zeros(restart+1,restart);
cs = zeros(restart,1);
sn = zeros(restart,1);
g = zeros(restart+1,1);

flag = 1;
initres = relres;
prevrelres = Inf;

for cycle = 1:maxit
    beta = norm(r);
    if beta == 0
        flag = 0;
        relres = 0;
        break
    end
    V(:,1) = r/beta;
    g(:) = 0;
    g(1) = beta;
    H(:) = 0;

    k = 0;
    for j = 1:restart
        % Right preconditioning: the Krylov space is that of A*M^-1, so the
        % basis vector is preconditioned before, not after, the product.
        z = U\(L\V(:,j));
        w = A*z;
        iter = iter + 1;

        wnorm0 = norm(w);
        % Classical Gram-Schmidt, twice; see the header.
        for pass = 1:2
            h = V(:,1:j)'*w;
            H(1:j,j) = H(1:j,j) + h;
            w = w - V(:,1:j)*h;
        end
        hnext = norm(w);
        H(j+1,j) = hnext;

        k = j;
        breakdown = hnext <= BREAKDOWN_TOL*wnorm0;
        if ~breakdown
            V(:,j+1) = w/hnext;
        end

        % Apply the accumulated Givens rotations to the new Hessenberg column,
        % then annihilate its subdiagonal entry with a fresh rotation.
        for i = 1:(j-1)
            t1 = cs(i)*H(i,j) + sn(i)*H(i+1,j);
            H(i+1,j) = -sn(i)*H(i,j) + cs(i)*H(i+1,j);
            H(i,j) = t1;
        end
        denom = hypot(H(j,j),H(j+1,j));
        if denom == 0
            cs(j) = 1;
            sn(j) = 0;
        else
            cs(j) = H(j,j)/denom;
            sn(j) = H(j+1,j)/denom;
        end
        H(j,j) = cs(j)*H(j,j) + sn(j)*H(j+1,j);
        H(j+1,j) = 0;
        g(j+1) = -sn(j)*g(j);
        g(j) = cs(j)*g(j);

        % abs(g(j+1)) is the residual norm of the least-squares problem, so the
        % inner loop stops without forming the iterate.
        relres = abs(g(j+1))/bnorm;
        if relres <= tol || breakdown
            break
        end
    end

    % Back-substitute on the rotated Hessenberg system, then map the correction
    % back through the preconditioner.
    y = zeros(k,1);
    for i = k:-1:1
        s = g(i);
        if i < k
            s = s - H(i,(i+1):k)*y((i+1):k);
        end
        if H(i,i) == 0
            y(i) = 0;
        else
            y(i) = s/H(i,i);
        end
    end
    x = x + U\(L\(V(:,1:k)*y));

    % The true residual, recomputed rather than carried: the rotated one drifts
    % from it once the factorization is inexact.
    r = b - A*x;
    relres = norm(r)/bnorm;

    if relres <= tol
        flag = 0;
        break
    end
    % Growing two orders of magnitude past the initial residual is divergence.
    if ~(relres < 1e2*initres)
        flag = 3;
        break
    end
    % A restart cycle that does not reduce the residual cannot be repaired by
    % running the same cycle again, since GMRES(m) restarts from the same space.
    if relres >= prevrelres*(1-1e-12)
        flag = 3;
        break
    end
    prevrelres = relres;
end

if any(~isfinite(x))
    flag = 3;
    relres = Inf;
    return
end
if relres <= tol
    flag = 0;
end
end
