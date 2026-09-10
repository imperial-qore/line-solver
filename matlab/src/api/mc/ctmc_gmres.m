function [x,flag,relres,iter] = ctmc_gmres(A,b,tol,restart,maxit,x0)
% [X,FLAG,RELRES,ITER]=CTMC_GMRES(A,B,TOL,RESTART,MAXIT,X0)
%
% Solve the sparse nonsymmetric linear system A*x=b by restarted GMRES with an
% ILUT right preconditioner, falling back to a Jacobi preconditioner when the
% incomplete factorization breaks down. This is the iterative counterpart of the
% direct sparse solve used by CTMC_SOLVE, intended for generators whose LU
% fill-in exceeds available memory.
%
% Two preparation steps are not optional on a generator. Rows are equilibrated
% to unit max norm, so the O(1) normalization row does not mix with rows
% carrying rates of a different magnitude. The states are then reordered by
% reverse Cuthill-McKee: in the natural ordering of a birth-death chain the
% unpivoted elimination has growth factor (mu/lambda)^n, which overflows by a
% few thousand states, and a fill-reducing ordering rather than pivoting is what
% removes it.
%
% A is the already-assembled coefficient matrix; no CTMC-specific processing is
% performed here, so the same kernel serves the stochastic complementation and
% aggregation kernels.
%
% Defaults: TOL=1e-12, RESTART=min(n,50), MAXIT=ceil(n/RESTART), X0=ones(n,1)/n.
% TOL is a linear-solve residual and is therefore much tighter than the
% fixed-point tolerance options.iter_tol.
%
% FLAG follows the MATLAB GMRES convention: 0 converged, 1 iteration limit
% reached, 2 preconditioner ill-conditioned, 3 stagnation. Callers must check
% it and fall back to the direct solve when it is nonzero.
%
% ITER counts matrix-vector products with A, i.e. total inner iterations across
% all restart cycles, so that the four codebases report one comparable scalar.
%
% The recurrence is in the private GMRES_ITERATE rather than the built-in GMRES,
% which applies the preconditioner on the LEFT; see that function for what that
% costs on a generator.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Relative threshold below which the incomplete factorization discards a fill-in
% entry.
%
% ILU TYPE: CROUT, NOT ILUTP. Pivoting by column destroys the band that reverse
% Cuthill-McKee just created: on the M/M/1/K generator at K = 5000 the 'ilutp'
% factor carries 4,981,518 nonzeros against 30,015 for 'crout' at the same drop
% tolerance, a factor of 166 and 20%% of a dense 5001x5001. The lean factor is
% usable here only because GMRES_ITERATE applies it on the RIGHT; the built-in
% GMRES applies it on the left and stagnates on it. This is also what the Java,
% C++ and Python kernels use, whose ILUT does not pivot either (scipy is asked
% for permc_spec='NATURAL').
ILUT_DROP_TOL = 1e-4;

n = size(A,1);
if nargin<3 || isempty(tol) || tol<=0
    tol = 1e-12;
end
if nargin<4 || isempty(restart) || restart<=0
    restart = min(n,50);
end
restart = min(restart,n);
if nargin<5 || isempty(maxit) || maxit<=0
    maxit = ceil(n/restart);
end
maxit = max(1,min(maxit,n));
if nargin<6 || isempty(x0)
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

[xp,flag,relres,iter] = gmres_iterate(A,L,U,b,x0,tol,restart,maxit);

% Undo the reverse Cuthill-McKee permutation.
x = zeros(n,1);
x(p) = xp;
end
