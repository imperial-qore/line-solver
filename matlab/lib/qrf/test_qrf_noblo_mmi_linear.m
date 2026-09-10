function test_qrf_noblo_mmi_linear
% Regression test for the QRF no-blocking MMI polytope (linear assembly).
%
% No LINE dependencies. The single-phase block below needs no LP solver either:
% it uses feasibility of a known-exact point as the oracle, which is what
% actually discriminates the polytope (an fmincon iterate does not). The
% MULTI-PHASE block that follows does run the solver, because the two defects it
% covers were both in the solve rather than in the polytope.
%
% The instance is the 2-station closed network with exponential service and
% load-dependent rates alpha(i,n), for which the exact stationary law is the
% birth-death chain
%
%     mu(2)*alpha(2,N-n) * g(n) = mu(1)*alpha(1,n+1) * g(n+1),   g(n)=P(n1=n).
%
% Three properties are asserted:
%   1. q is 5-D, [M M max(K) max(K) N+1], per qrboundsrsrd_skel.mod.
%   2. The exact point is feasible.
%   3. A point built from the WRONG rates is infeasible. This is the property
%      that failed before the 6-subscript-write / 4-subscript-read defect was
%      fixed: q then read as all zeros, every q-weighted row (THM1, THM30,
%      THM3, THM3f) was a vacuous 0=0 row, and the polytope did not depend on
%      mu or alpha at all.

N = 3; M = 2; rt = [0 1; 1 0];
tol = 1e-9;

cases = { ...
    struct('mu',[1.0 1.0],'alpha',[1 1 1]), ...
    struct('mu',[4.0 1.0],'alpha',[1 1 1]), ...
    struct('mu',[0.5 1.0],'alpha',[1 1 1]), ...
    struct('mu',[1.0 1.0],'alpha',[1 2 3]), ...
    struct('mu',[1.0 1.0],'alpha',[1 3 9])};

for t = 1:numel(cases)
    mu = cases{t}.mu;
    alpha = ones(M,N); alpha(1,:) = cases{t}.alpha;
    MAPs = {{-mu(1), mu(1)}; {-mu(2), mu(2)}};

    [~,~,~,~,lp] = qrf_noblo_mmi_linear(MAPs,N,rt,alpha,true);

    % 1. arity of q
    q = local_q(MAPs,N,rt,alpha);
    assert(isequal(size(q),[M M 1 1 N+1]), ...
        'case %d: q must be 5-D [M M K K N+1], got %s', t, mat2str(size(q)));

    % 2. the exact point is feasible
    g = local_exact(N,mu,alpha);
    x = local_point(N,M,g);
    [rEq,rIn] = local_residuals(lp,x);
    assert(rEq < tol, 'case %d: exact point violates an equality by %g', t, rEq);
    assert(rIn < tol, 'case %d: exact point violates an inequality by %g', t, rIn);

    % 3. a point built from the wrong rates is infeasible
    gBad = local_exact(N,[mu(1)*3 mu(2)],alpha);
    if max(abs(gBad-g)) > 1e-6
        xBad = local_point(N,M,gBad);
        rEqBad = local_residuals(lp,xBad);
        assert(rEqBad > 1e-3, ...
            ['case %d: a point built from the wrong service rates is feasible ' ...
             '(residual %g). The polytope does not depend on the rates.'], t, rEqBad);
    end
end

% ---- MULTI-PHASE (K > 1), the instances the block above cannot reach ----
%
% Until 2026-09-01 this file covered only the SINGLE-PHASE exponential instance,
% so neither of the two defects below was visible to it:
%
%   K = [2,2]  the inlined phase 1 was a bare minimum-norm quadprog with no
%              linprog fallback, and it stalls at an equality residual of
%              7.9e+00 on this polytope. The arm RAISED on a well-posed model.
%   K = [2,1]  sub_qrfvar filled the decision vector compactly while deltap2
%              addressed it with max(K)-padded strides, so every family indexed
%              the wrong columns and the arm returned UN = [0,0].
%
% Two oracles, neither needing a hand-computed CTMC. Visits are equal on this
% cycle, so U(i) = X*s(i) and the utilizations carry the DEMAND RATIO exactly;
% and the queue lengths must carry the population. The UN values additionally
% pin against native python, which returns the same numbers to the digits shown.
mpTol = 1e-6;
erl2 = @(mean) {[-2/mean, 2/mean; 0, -2/mean], [0, 0; 2/mean, 0]};
expo = @(mean) {-1/mean, 1/mean};
mp = { ...
    struct('name','K=[2,2]','MAPs',{{erl2(1/2); erl2(1/4)}},'UN',[0.75 0.375]), ...
    struct('name','K=[2,1]','MAPs',{{erl2(1/2); expo(1/4)}},'UN',[0.80 0.400])};

Nmp = 2; rtmp = [0 1; 1 0]; alphamp = ones(2,Nmp);
for t = 1:numel(mp)
    [UNmp,QNmp] = qrf_noblo_mmi_linear(mp{t}.MAPs,Nmp,rtmp,alphamp);
    assert(all(isfinite(UNmp)) && any(UNmp > 0), ...
        '%s: solver returned no utilization (UN = %s)', mp{t}.name, mat2str(UNmp));
    assert(abs(sum(QNmp) - Nmp) < mpTol, ...
        '%s: queue lengths carry %g jobs, not N = %d', mp{t}.name, sum(QNmp), Nmp);
    assert(abs(UNmp(1)/UNmp(2) - 2) < mpTol, ...
        '%s: demands are 2:1 but UN = %s', mp{t}.name, mat2str(UNmp,8));
    assert(max(abs(UNmp(:).' - mp{t}.UN)) < mpTol, ...
        '%s: UN = %s, expected %s', mp{t}.name, mat2str(UNmp,8), mat2str(mp{t}.UN));
end

fprintf('test_qrf_noblo_mmi_linear: %d single-phase + %d multi-phase cases passed\n', ...
    numel(cases), numel(mp));
end

function q = local_q(MAPs,N,rt,alpha)
% Rebuild q exactly as qrf_noblo_mmi_linear does, to assert its arity.
M = length(MAPs);
for i=1:M, K(i)=size(MAPs{i}{1},1); end
q = zeros(M,M,max(K),max(K),N+1);
for i = 1:M, for j = 1:M, for k = 1:K(i), for h = 1:K(i), for ni = 1:N
    mu = MAPs{i}{2}(k,h);
    if k==h, v = 0; else, v = MAPs{i}{1}(k,h); end
    if j ~= i
        q(i,j,k,h,1+ni) = rt(i,j)*mu*alpha(i,ni);
    else
        q(i,j,k,h,1+ni) = v*alpha(i,ni)+rt(i,i)*mu*alpha(i,ni);
    end
end, end, end, end, end
end

function g = local_exact(N,mu,alpha)
% g(1+n) = P(n1 = n) for the 2-station load-dependent birth-death chain.
g = zeros(1,N+1); g(1) = 1;
for n = 0:N-1
    g(1+n+1) = g(1+n) * (mu(2)*alpha(2,N-n)) / (mu(1)*alpha(1,n+1));
end
g = g / sum(g);
end

function x = local_point(N,M,g)
% Assemble the decision vector in the layout used by deltap2/deltae, for the
% single-phase two-station case: p2 is supported on n1+n2 = N.
MR = 1; Km = 1;
nx = M*(N+1)*Km*M*(N+1)*Km*MR + M*Km;
x = zeros(nx,1);
gi = @(i,n) g(1 + (i==1)*n + (i==2)*(N-n));   % P(n_i = n)
for j = 1:M, for nj = 0:N, for i = 1:M, for ni = 0:N
    if i == j
        if nj == ni, val = gi(i,ni); else, val = 0; end
    else
        if nj + ni == N, val = gi(i,ni); else, val = 0; end
    end
    if val ~= 0
        x(local_idx(N,M,j,1+nj,i,1+ni)) = val;
    end
end, end, end, end
% e(i,k) = P(n_i >= 1)
base = M*(N+1)*Km*M*(N+1)*Km*MR;
for i = 1:M
    x(base + i) = 1 - gi(i,0);
end
end

function idx = local_idx(N,M,j,nj,i,ni)
MR = 1; Km = 1; k = 1; h = 1; m = 1;
idx = (N+1)*Km*M*(N+1)*Km*MR*(j-1);
idx = idx + Km*M*(N+1)*Km*MR*(nj-1);
idx = idx + M*(N+1)*Km*MR*(k-1);
idx = idx + (N+1)*Km*MR*(i-1);
idx = idx + Km*MR*(ni-1);
idx = idx + MR*(h-1);
idx = idx + m;
end

function [rEq,rIn] = local_residuals(lp,x)
rEq = 0; rIn = 0;
if ~isempty(lp.Aeq)
    rEq = max(abs(full(lp.Aeq*x - lp.beq(:))));
end
if ~isempty(lp.A) && size(lp.A,1) > 0
    rIn = max([0; full(lp.A*x - lp.b(:))]);
end
end
