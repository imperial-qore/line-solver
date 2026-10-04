function [x,flag,relres,iter] = ctmc_bicgstab(A,b,tol,maxit,x0)
% [X,FLAG,RELRES,ITER]=CTMC_BICGSTAB(A,B,TOL,MAXIT,X0)
%
% Solve the sparse nonsymmetric linear system A*x=b by the stabilized biconjugate
% gradient method of van der Vorst (1992), with the same equilibration and
% reordering as CTMC_GMRES and a non-pivoting incomplete factorization (see the
% ILU note below, which is where the two kernels part company). It is the short-recurrence counterpart of
% that kernel: work and storage per iteration are constant rather than growing
% with the Krylov dimension, so it does not restart and does not lose the
% optimality that restarting costs GMRES. Where GMRES(m) stagnates because the
% useful subspace is wider than m, this converges; where it does not, GMRES(m)
% is the more robust of the two, hence the order in which CTMC_SOLVE tries them.
%
% Preparation is not optional on a generator and is identical to CTMC_GMRES:
% rows are equilibrated to unit max norm, so the O(1) normalization row does not
% mix with rows carrying rates of a different magnitude, and the states are then
% reordered by reverse Cuthill-McKee, since in the natural ordering of a
% birth-death chain the unpivoted elimination has growth factor (mu/lambda)^n.
%
% A is the already-assembled coefficient matrix; no CTMC-specific processing is
% performed here, so the same kernel serves the stochastic complementation and
% aggregation kernels.
%
% Defaults: TOL=1e-12, MAXIT=min(n,200), X0=ones(n,1)/n. TOL is a linear-solve
% residual and is therefore much tighter than the fixed-point tolerance
% options.iter_tol.
%
% FLAG follows the MATLAB BICGSTAB convention: 0 converged, 1 iteration limit
% reached, 2 preconditioner ill-conditioned, 3 stagnation, 4 a scalar quantity
% became too small or too large to continue. Callers must check it and fall back
% to another solve when it is nonzero.
%
% ITER counts matrix-vector products with A: two per complete iteration, and one
% when the iteration converges at its half step. Reporting the product count
% rather than the iteration count is what makes it comparable with the ITER of
% CTMC_GMRES and across the four codebases, whose iteration bookkeeping differs.
%
% The recurrence is in the private BICGSTAB_ITERATE rather than the built-in
% BICGSTAB, which applies the preconditioner on the LEFT and breaks down on the
% systems a generator produces; see that function for the measurement.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Relative threshold below which the incomplete factorization discards a fill-in
% entry.
%
% ILU TYPE: CROUT, NOT ILUTP, AND THE DIFFERENCE IS MEASURED. CTMC_GMRES asks
% MATLAB's `ilu` for 'ilutp', which pivots by column; on a banded generator that
% pivoting destroys the band the reverse Cuthill-McKee ordering just created, and
% the factor of the M/M/1/K balance system at K = 5000 comes out with 4,981,518
% nonzeros against 30,015 for 'crout' at the same drop tolerance -- a factor 166,
% and 20%% of a dense 5001x5001. With that factor the BiCGSTAB recurrence breaks
% down at the second product (flag 4, relative residual 1.7e-04); with the Crout
% factor it converges in four products to 2.7e-16, and the answer agrees with the
% closed-form truncated geometric to 2.3e-13. 'crout' is also what the Java, C++
% and Python kernels use: their ILUT does not pivot either (scipy is asked for
% permc_spec='NATURAL'), so this brings MATLAB closer to them, not further.
ILUT_DROP_TOL = 1e-4;

% Default cap on complete iterations. BiCGSTAB storage is O(n) regardless of the
% count, so the cap exists to bound time, not memory.
BICGSTAB_DEFAULT_MAXIT = 200;

n = size(A,1);
if nargin<3 || isempty(tol) || tol<=0
    tol = 1e-12;
end
if nargin<4 || isempty(maxit) || maxit<=0
    maxit = min(n,BICGSTAB_DEFAULT_MAXIT);
end
maxit = max(1,min(maxit,n));
if nargin<5 || isempty(x0)
    x0 = ones(n,1)/n;
end
x0 = reshape(x0,n,1);
b = reshape(b,n,1);

if ~issparse(A)
    A = sparse(A);
end

% Row equilibration. Row scaling leaves the solution unchanged.
rownorm = full(max(abs(A),[],2));
rownorm(rownorm==0) = 1;
A = spdiags(1./rownorm,0,n,n)*A;
b = b./rownorm;

% Reverse Cuthill-McKee on the symmetrized pattern.
p = symrcm(spones(A)+spones(A)');
A = A(p,p);
b = b(p);
x0 = x0(p);

L = [];
U = [];
try
    [L,U] = ilu(A,struct('type','crout','droptol',ILUT_DROP_TOL,'udiag',1));
    if any(~isfinite(nonzeros(L))) || any(~isfinite(nonzeros(U)))
        L = [];
        U = [];
    end
catch
    L = [];
    U = [];
end

if isempty(L)
    % Jacobi fallback. A zero diagonal entry would make the preconditioner
    % singular, so those rows are left unscaled rather than inverted.
    d = full(diag(A));
    d(d==0) = 1;
    L = spdiags(1./d,0,n,n);
    U = speye(n);
end

[x,flag,relres,iter] = bicgstab_iterate(A,L,U,b,x0,tol,maxit);

% Undo the reverse Cuthill-McKee permutation.
xp = x;
x = zeros(n,1);
x(p) = xp;
end
