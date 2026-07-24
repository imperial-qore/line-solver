% Analytic validation of ctmc_gmres.
%
% The oracle is the M/M/1/K stationary distribution, a truncated geometric,
% rather than a recorded baseline: a recorded value cannot catch an error that
% the reference implementation shares.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

failures = 0;

%% Analytic oracle. K=20000 is past the size at which a direct factorization is
%% the intended method, and rho<1 is the case in which an unpreconditioned or
%% naturally ordered elimination overflows.
for rho = [0.7 1.0]
    for K = [10 20000]
        [A,b] = local_normalized_system(local_mm1k(rho,1.0,K));
        [x,flag,relres] = ctmc_gmres(A,b);
        pi_exact = (rho.^(0:K))';
        pi_exact = pi_exact/sum(pi_exact);
        err = max(abs(x-pi_exact));
        ok = (flag==0) && (err < 1e-9);
        failures = failures + ~ok;
        fprintf('mm1k rho=%.2f K=%d flag=%d relres=%.2e err=%.2e %s\n', ...
            rho,K,flag,relres,err,local_verdict(ok));
    end
end

%% The dispatch in ctmc_solve must be numerically invisible, so the two methods
%% have to agree far below any tolerance a caller would notice.
cfgs = [3 0.5; 10 0.85; 200 0.85];
for c = 1:size(cfgs,1)
    K = cfgs(c,1); rho = cfgs(c,2);
    [A,b] = local_normalized_system(local_mm1k(rho,1.0,K));
    xd = full(A\b);
    [x,flag] = ctmc_gmres(A,b);
    dev = max(abs(x-xd));
    ok = (flag==0) && (dev < 1e-9);
    failures = failures + ~ok;
    fprintf('vs direct K=%d rho=%.2f flag=%d dev=%.2e %s\n',K,rho,flag,dev,local_verdict(ok));
end

%% A caller may only trust the answer when the flag is zero. One restart cycle of
%% dimension one cannot converge on a 500-state chain, and the kernel has to say
%% so rather than return the iterate it happens to hold.
[A,b] = local_normalized_system(local_mm1k(0.9,1.0,500));
[~,flag,relres] = ctmc_gmres(A,b,1e-12,1,1);
ok = (flag~=0) && (relres > 1e-12);
failures = failures + ~ok;
fprintf('non-convergence flag=%d relres=%.2e %s\n',flag,relres,local_verdict(ok));

%% The two methods must agree far below any tolerance a caller would notice, and
%% the answer must not jump as a model grows past the dispatch threshold.
opts_d = struct('method','direct','verbose',0,'iter_max',1000);
opts_g = struct('method','gmres','verbose',0,'iter_max',1000);
for K = [50 400]
    Q = full(local_mm1k(0.8,1.0,K));
    pd = ctmc_solve(Q,opts_d);
    pg = ctmc_solve(Q,opts_g);
    pi_exact = (0.8.^(0:K));
    pi_exact = pi_exact/sum(pi_exact);
    dev = max(abs(pd(:)-pg(:)));
    err = max(abs(pg(:)-pi_exact(:)));
    ok = (dev < 1e-9) && (err < 1e-9);
    failures = failures + ~ok;
    fprintf('dispatch K=%d dev=%.2e err=%.2e %s\n',K,dev,err,local_verdict(ok));
end

%% A generator with no diagonal entries breaks the incomplete factorization. The
%% kernel must fall back rather than propagate the breakdown.
Z = sparse([1 2 3],[2 3 1],[1 1 1],3,3);
bz = sparse(3,1); bz(1) = 1;
x = ctmc_gmres(Z,bz);
ok = all(isfinite(x));
failures = failures + ~ok;
fprintf('zero diagonal finite=%d %s\n',ok,local_verdict(ok));

fprintf('\ntest_ctmc_gmres: %d failure(s)\n',failures);

function Q = local_mm1k(lam,mu,K)
% Generator of an M/M/1/K queue.
n = K+1;
I = []; J = []; V = [];
for i=1:n
    if i<n, I(end+1)=i; J(end+1)=i+1; V(end+1)=lam; end %#ok<AGROW>
    if i>1, I(end+1)=i; J(end+1)=i-1; V(end+1)=mu; end %#ok<AGROW>
end
Q = sparse(I,J,V,n,n);
Q = Q - spdiags(sum(Q,2),0,n,n);
end

function [A,b] = local_normalized_system(Q)
% Assembles the linear system ctmc_solve poses: the last column of Q is replaced
% by ones to carry the normalization, and the transposed system is solved
% against e_n.
n = size(Q,1);
Q(:,end) = 1;
A = Q';
b = sparse(n,1);
b(end) = 1;
end

function s = local_verdict(ok)
if ok
    s = 'PASS';
else
    s = 'FAIL';
end
end
