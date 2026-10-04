function tf = mam_chain_arrival_is_markovian(sn, c)
% MAM_CHAIN_ARRIVAL_IS_MARKOVIAN  True unless chain C is fed by an ME/RAP source.
%
% The system arrival MMAP of an open chain is assembled from the process of each
% of the chain's classes at its reference station, which for an open chain is the
% source. When any of those is a matrix-exponential (ME) or rational arrival
% process (RAP), the assembled (D0,D1) pair is legitimately non-Markovian: D0 may
% carry negative off-diagonal entries.
%
% Callers use this to decide whether MMAP_NORMALIZE may be applied. That routine
% clips negative entries to zero and re-derives the diagonal, which repairs
% numerical noise on a genuine MAP but on an ME or RAP substitutes a DIFFERENT,
% Markovian process with different autocorrelation. The distinction has to be
% made on the declared process type rather than on the sign pattern, because a
% MAP perturbed by roundoff also shows small negative entries and does need the
% repair.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

tf = true;
inchain = sn.inchain{c};
if isempty(inchain)
    return;
end
if ~isinf(sum(sn.njobs(inchain)))
    return;   % closed chain: a Poisson surrogate is used, always Markovian
end
ist = sn.refstat(inchain(1));
for k = inchain(:)'
    if sn.procid(ist, k) == ProcessType.ME || sn.procid(ist, k) == ProcessType.RAP
        tf = false;
        return;
    end
end
end
