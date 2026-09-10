function [StartN, PreemptN] = solver_ctmc_auxrates(DfiltAux, probSysState, wset, M, K)
% [STARTN, PREEMPTN] = SOLVER_CTMC_AUXRATES(DFILTAUX, PROBSYSSTATE, WSET, M, K)
%
% Reduce the derived START/PREEMPT filtrations to (station x class) rates:
% StartN(i,r) = pi * F_start{i,r} * e, the long-run number of class-r service
% starts per unit time at station i, and likewise for preemptions. Summing the
% row of the filtration and weighting by the stationary probability is the
% same reduction the departure rates use, so the two are directly comparable:
% StartN == TN + PreemptN at a lossless station.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

StartN = zeros(M,K);
PreemptN = zeros(M,K);
if isempty(DfiltAux) || ~isstruct(DfiltAux) || ~isfield(DfiltAux,'start')
    return
end
pi_row = probSysState(:)';
for i = 1:min(M, size(DfiltAux.start,1))
    for r = 1:min(K, size(DfiltAux.start,2))
        Fs = DfiltAux.start{i,r};
        if ~isempty(Fs) && nnz(Fs) > 0
            rowsum = full(sum(Fs,2));
            StartN(i,r) = pi_row * rowsum(wset);
        end
        Fp = DfiltAux.preempt{i,r};
        if ~isempty(Fp) && nnz(Fp) > 0
            rowsum = full(sum(Fp,2));
            PreemptN(i,r) = pi_row * rowsum(wset);
        end
    end
end
end
