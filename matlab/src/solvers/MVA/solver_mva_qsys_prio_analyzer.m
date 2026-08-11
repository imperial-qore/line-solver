function [Q,U,R,T,C,X,lG,runtime,totiter,actualmethod] = solver_mva_qsys_prio_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME,ITER,ACTUALMETHOD] = SOLVER_MVA_QSYS_PRIO_ANALYZER(SN, OPTIONS)
%
% Exact non-preemptive priority (HOL) analyzer for a single open M/G/1 queue
% with Poisson per-class arrivals: dispatches to the Cobham formula
% (qsys_mg1_prio) instead of the AMVA preemptive shadow-server approximation,
% which underestimates the waiting time of every class.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0 = tic;
totiter = 0;
lG = 0;
actualmethod = 'mg1.prio'; % exact Cobham formula, reported in the solver banner
K = sn.nclasses;
M = sn.nstations;

Q = zeros(M,K); U = zeros(M,K);
R = zeros(M,K); T = zeros(M,K);
C = zeros(M,K); X = zeros(1,K);

source_ist = sn.nodeToStation(sn.nodetype == NodeType.Source);
queue_ist = sn.nodeToStation(sn.nodetype == NodeType.Queue);

% Order classes by priority (lowest value = highest priority, stable)
[~, order] = sort(sn.classprio(:)', 'ascend');

lambdaOrd = sn.rates(source_ist, order);
muOrd = sn.rates(queue_ist, order);
scvOrd = sn.scv(queue_ist, order);
scvOrd(~isfinite(scvOrd) | scvOrd <= 0) = 1;
csOrd = sqrt(scvOrd);

W = qsys_mg1_prio(lambdaOrd, muOrd, csOrd);

for j = 1:K
    r = order(j);
    lam = lambdaOrd(j);
    if lam <= 0 || ~isfinite(muOrd(j))
        continue
    end
    R(queue_ist, r) = W(j);
    C(queue_ist, r) = W(j);
    X(r) = lam;
    U(queue_ist, r) = lam / muOrd(j);
    T(queue_ist, r) = lam;
    Q(queue_ist, r) = lam * W(j);
    T(source_ist, r) = lam;
end

runtime = toc(T0);
end
