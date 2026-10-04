function [X,flag] = ctmc_bicgstab_multi(A,B,tol,maxit)
% [X,FLAG]=CTMC_BICGSTAB_MULTI(A,B,TOL,MAXIT)
%
% Solve A*X=B for every column of B by stabilized biconjugate gradients, reusing
% one ILUT factorization across all of them and starting each column from the
% previous solution. This is the shape of the stochastic complement, whose
% right-hand side is a whole block of the generator: refactorizing per column
% would cost more than the direct solve it replaces.
%
% Preparation follows CTMC_BICGSTAB: rows are equilibrated to unit max norm and
% the states reordered by reverse Cuthill-McKee, without which the unpivoted
% elimination overflows on a chain of a few thousand states. The iteration is
% the private BICGSTAB_ITERATE, the same one CTMC_BICGSTAB runs.
%
% FLAG is 0 only if every column converged. On any other value X must be
% discarded and the caller must fall back to another solve; returning a partial
% block would leave the fallback ambiguous.

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

% Default cap on complete iterations, as in CTMC_BICGSTAB.
BICGSTAB_DEFAULT_MAXIT = 200;

n = size(A,1);
nrhs = size(B,2);
if nargin<3 || isempty(tol) || tol<=0
    tol = 1e-12;
end
if nargin<4 || isempty(maxit) || maxit<=0
    maxit = min(n,BICGSTAB_DEFAULT_MAXIT);
end
maxit = max(1,min(maxit,n));

if ~issparse(A)
    A = sparse(A);
end

rownorm = full(max(abs(A),[],2));
rownorm(rownorm==0) = 1;
A = spdiags(1./rownorm,0,n,n)*A;
B = B./rownorm;

p = symrcm(spones(A)+spones(A)');
A = A(p,p);
B = B(p,:);

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
    d = full(diag(A));
    d(d==0) = 1;
    L = spdiags(1./d,0,n,n);
    U = speye(n);
end

Xp = zeros(n,nrhs);
guess = ones(n,1)/n;
for c = 1:nrhs
    [xc,fc] = bicgstab_iterate(A,L,U,full(B(:,c)),guess,tol,maxit);
    if fc ~= 0
        X = [];
        flag = fc;
        return
    end
    Xp(:,c) = xc;
    guess = xc;
end

X = zeros(n,nrhs);
X(p,:) = Xp;
flag = 0;
end
