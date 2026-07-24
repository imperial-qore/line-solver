function [QN, UN, RN, TN] = solver_mam_ldqbd_avg(ld, piflat, levelOf)
% SOLVER_MAM_LDQBD_AVG Map an LD-QBD state distribution to mean metrics.
%
% [QN,UN,RN,TN] = SOLVER_MAM_LDQBD_AVG(LD, PIFLAT, LEVELOF) takes a probability
% vector PIFLAT over the flat LD-QBD state space (with per-state level LEVELOF
% from solver_mam_ldqbd_flatten) and returns the per-(station,class) mean queue
% length, utilization, response time and throughput for the single-class
% Delay/Queue (closed) or Source/Queue (open) model described by LD.
%
% Mirrors the steady-state metric formulas in solver_mam_ldqbd, applied to an
% arbitrary (e.g. transient-averaged) distribution rather than the stationary
% one. Used by the SolverENV state-vector analyzer's MAM backend.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Nlev = ld.Nlev;
M = ld.M;
qi = ld.queueIdx;
ri = ld.refIdx;
c = ld.nServers;

% Aggregate the flat distribution to per-level probabilities.
piflat = piflat(:)';
piflat(piflat < 0) = 0;
if sum(piflat) > 0
    piflat = piflat / sum(piflat);
end
pLevel = zeros(1, Nlev + 1);
for n = 0:Nlev
    pLevel(n+1) = sum(piflat(levelOf == n));
end

mean_queue = (0:Nlev) * pLevel(:);

% Utilization: P(busy) for a (load-dependent) single server, else average
% fraction of c servers in use.
if ld.hasLLD || c == 1
    util = 1 - pLevel(1);
else
    util = 0;
    for n = 1:Nlev
        util = util + (min(n, c) / c) * pLevel(n+1);
    end
end

QN = zeros(M, 1);
UN = zeros(M, 1);
RN = zeros(M, 1);
TN = zeros(M, 1);

if ld.isOpen
    X = ld.lambda_eff * (1 - pLevel(Nlev + 1));
    if X > 0, R_queue = mean_queue / X; else, R_queue = 0; end
    QN(ri) = 0;   UN(ri) = 0;    RN(ri) = 0;       TN(ri) = X;
    QN(qi) = mean_queue; UN(qi) = util; RN(qi) = R_queue; TN(qi) = X;
else
    mean_delay = ld.N - mean_queue;
    X = mean_delay * ld.lambda_eff;
    if X > 0, R_queue = mean_queue / X; else, R_queue = 0; end
    R_delay = 1 / ld.delayRate;
    QN(ri) = mean_delay; UN(ri) = mean_delay; RN(ri) = R_delay; TN(ri) = X;
    QN(qi) = mean_queue; UN(qi) = util / c;   RN(qi) = R_queue; TN(qi) = X;
end
end
