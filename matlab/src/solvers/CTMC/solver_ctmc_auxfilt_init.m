function DfiltAux = solver_ctmc_auxfilt_init(sn, Q)
% DFILTAUX = SOLVER_CTMC_AUXFILT_INIT(SN, Q)
%
% Allocate the derived START/PREEMPT filtrations: one all-zero sparse matrix
% of the shape of Q per (station, class). DFILTAUX.start{i,r}(s,ns) will hold
% the rate at which the transition s -> ns starts a class-r service at station
% i, and DFILTAUX.preempt{i,r}(s,ns) the rate at which it pushes a class-r job
% in service back into the buffer there.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;
R = sn.nclasses;
DfiltAux.start = cell(M,R);
DfiltAux.preempt = cell(M,R);
Z = 0*Q;
for i = 1:M
    for r = 1:R
        DfiltAux.start{i,r} = Z;
        DfiltAux.preempt{i,r} = Z;
    end
end
end
