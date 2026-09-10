function [x,flag,relres,iter] = bicgstab_iterate(A,L,U,b,x0,tol,maxit)
% [X,FLAG,RELRES,ITER]=BICGSTAB_ITERATE(A,L,U,B,X0,TOL,MAXIT)
%
% Right-preconditioned stabilized biconjugate gradients (van der Vorst, 1992) on
% an ALREADY equilibrated, reordered and factorized system. Internal to
% CTMC_BICGSTAB and CTMC_BICGSTAB_MULTI, which share it so that one
% factorization can serve a whole block of right-hand sides.
%
% THE BUILT-IN BICGSTAB IS NOT USED, and the reason is not style. MATLAB applies
% an (L,U) pair on the LEFT, so it iterates on M\A and reports the
% PRECONDITIONED residual, while the Java, Python and C++ ports of this kernel
% expand the Krylov space of A*M^-1 and report the TRUE residual. On a generator
% the difference is not cosmetic: the left-preconditioned recurrence breaks down
% (flag 4, relative residual 6e-05) on the M/M/1/K balance system at K = 5000,
% where the right-preconditioned one converges to 1e-13. Writing the recurrence
% here makes all four codebases run the same algorithm and report the same FLAG,
% RELRES and ITER on the same system.
%
% The preconditioner is applied as U\(L\v); L and U are the ILUT factors, or the
% Jacobi pair (diag(1./d), I) when the incomplete factorization broke down.
%
% FLAG follows the MATLAB BICGSTAB convention: 0 converged, 1 iteration limit
% reached, 3 stagnation or divergence, 4 a scalar quantity became too small or
% too large to continue. ITER counts matrix-vector products with A: two per
% complete iteration, and one when the iteration converges at its half step.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Threshold below which rho or omega is treated as a Lanczos breakdown, taken
% relative to the norms whose product formed it.
BREAKDOWN_TOL = 1e-14;

n = size(A,1);
x = reshape(x0,n,1);
b = reshape(b,n,1);

bnorm = norm(b);
if bnorm == 0
    bnorm = 1;
end

r = b - A*x;
relres = norm(r)/bnorm;
iter = 0;
if relres <= tol
    flag = 0;
    return
end

% The shadow residual is fixed at the initial residual, the standard choice: any
% vector not orthogonal to r would do, and this one cannot be orthogonal to it.
rhat = r;
p = zeros(n,1);
v = zeros(n,1);
rho = 1;
alpha = 1;
omega = 1;
flag = 1;
bestrelres = relres;

for it = 1:maxit
    rhoNew = rhat'*r;
    % rho vanishing is the biorthogonality breakdown of the underlying Lanczos
    % process, not slow convergence: restarting with a fresh shadow vector would
    % discard the iterate, so the caller is told to use another method instead.
    if abs(rhoNew) <= BREAKDOWN_TOL*norm(rhat)*norm(r)
        flag = 4;
        break
    end
    if it == 1
        p = r;
    else
        if omega == 0
            flag = 4;
            break
        end
        beta = (rhoNew/rho)*(alpha/omega);
        p = r + beta*(p - omega*v);
    end
    rho = rhoNew;

    ph = U\(L\p);
    v = A*ph;
    iter = iter + 1;

    rhatv = rhat'*v;
    if rhatv == 0 || ~isfinite(rhatv)
        flag = 4;
        break
    end
    alpha = rho/rhatv;

    s = r - alpha*v;

    % Half-step convergence: s is the residual of x + alpha*ph, so a converged s
    % reaches the answer without the second matrix-vector product.
    snorm = norm(s);
    if snorm/bnorm <= tol
        x = x + alpha*ph;
        relres = snorm/bnorm;
        flag = 0;
        break
    end

    sh = U\(L\s);
    t = A*sh;
    iter = iter + 1;

    tt = t'*t;
    if tt == 0 || ~isfinite(tt)
        flag = 4;
        break
    end
    omega = (t'*s)/tt;

    x = x + alpha*ph + omega*sh;
    r = s - omega*t;

    relres = norm(r)/bnorm;
    if relres <= tol
        flag = 0;
        break
    end
    % omega vanishing stalls the update of x while leaving r finite, so the
    % iteration would spin without progress.
    if abs(omega) <= BREAKDOWN_TOL
        flag = 4;
        break
    end
    % BiCGSTAB residuals are non-monotone by construction, so an increase is not
    % by itself stagnation and the test is against the BEST residual seen rather
    % than the previous one. Growing two orders of magnitude past that best is
    % divergence, and continuing from such an iterate is not worth the products.
    if relres > 1e2*bestrelres
        flag = 3;
        break
    end
    bestrelres = min(bestrelres, relres);
end

if any(~isfinite(x))
    flag = 4;
    relres = Inf;
    return
end
if relres <= tol
    flag = 0;
end
end
