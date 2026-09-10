function [X,flag] = ctmc_gmres_multi(A,B,tol,restart,maxit)
% [X,FLAG]=CTMC_GMRES_MULTI(A,B,TOL,RESTART,MAXIT)
%
% Solve A*X=B for every column of B by restarted GMRES, reusing one ILUT
% factorization across all of them and starting each column from the previous
% solution. This is the shape of the stochastic complement, whose right-hand
% side is a whole block of the generator: refactorizing per column would cost
% more than the direct solve it replaces.
%
% Preparation follows CTMC_GMRES: rows are equilibrated to unit max norm and the
% states reordered by reverse Cuthill-McKee, without which the unpivoted
% elimination overflows on a chain of a few thousand states.
%
% FLAG is 0 only if every column converged. On any other value X must be
% discarded and the caller must fall back to the direct solve; returning a
% partial block would leave the fallback ambiguous.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Relative threshold below which the incomplete factorization discards a fill-in
% entry. CROUT, NOT ILUTP, and with the recurrence taken from the private
% GMRES_ITERATE rather than the built-in GMRES; see CTMC_GMRES for the fill
% measurement that decides it, which the block right-hand side only magnifies.
ILUT_DROP_TOL = 1e-4;

n = size(A,1);
nrhs = size(B,2);
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
    [xc,fc] = gmres_iterate(A,L,U,full(B(:,c)),guess,tol,restart,maxit);
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
