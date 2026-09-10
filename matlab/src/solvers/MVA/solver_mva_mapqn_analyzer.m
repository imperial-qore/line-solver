function [QN,UN,RN,TN,CN,XN,lG,runtime,lastiter,actualmethod] = solver_mva_mapqn_analyzer(sn, options)
% [QN,UN,RN,TN,CN,XN,LG,RUNTIME,LASTITER,ACTUALMETHOD] = SOLVER_MVA_MAPQN_ANALYZER(SN, OPTIONS)
%
% SolverMVA method 'amva.mapqn': the horizontal-cut mean value analysis of a
% closed multiclass network with one exponential delay station and one FCFS
% single-server queue whose class-r service is a MAP (mapqn_amva). The
% structural premise is checked by mva_mapqn_reason, the same predicate
% that offers and reports the method. Metrics: at the queue Q_r, U_r = X_r
% E[S_r] and R_r = Q_r / X_r; at the delay Q_r = U_r = X_r Z_r and R_r = Z_r;
% C_r = N_r / X_r. lG is NaN: the recursion carries no normalizing constant.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Tstart = tic;
reason = mva_mapqn_reason(sn);
if ~isempty(reason)
    line_error(mfilename, reason);
end
M = sn.nstations;
R = sn.nclasses;
isDelay = (sn.sched(:) == SchedStrategy.INF) | isinf(sn.nservers(:));
id = find(isDelay); iq = find(~isDelay);
N = round(sn.njobs(:)');
mu = ones(1, R); D0s = cell(1, R); D1s = cell(1, R);
for r = 1:R
    if N(r) > 0
        mu(r) = sn.rates(id, r);
        D0s{r} = sn.proc{iq}{r}{1};
        D1s{r} = sn.proc{iq}{r}{2};
    else
        D0s{r} = -1; D1s{r} = 1;   % absent class: inert single phase
    end
end
line_debug(options, 'amva.mapqn: %d classes, joint phase space of size %d', R, prod(cellfun(@(d) size(d, 1), D0s)));
[X, Qq, U] = mapqn_amva(mu, D0s, D1s, N);
QN = zeros(M, R); UN = zeros(M, R); RN = zeros(M, R); TN = zeros(M, R);
CN = zeros(1, R); XN = zeros(1, R);
for r = 1:R
    if N(r) <= 0 || X(r) <= 0
        continue
    end
    XN(r) = X(r);
    QN(iq, r) = Qq(r); UN(iq, r) = U(r); TN(iq, r) = X(r); RN(iq, r) = Qq(r) / X(r);
    QN(id, r) = X(r) / mu(r); UN(id, r) = X(r) / mu(r); TN(id, r) = X(r); RN(id, r) = 1 / mu(r);
    CN(r) = N(r) / X(r);
end
lG = NaN;
runtime = toc(Tstart);
lastiter = prod(N + 1);
actualmethod = 'amva.mapqn';
end
