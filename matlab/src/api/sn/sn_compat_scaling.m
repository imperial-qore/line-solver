function eta = sn_compat_scaling(compat, counts, rates, n)
% ETA = SN_COMPAT_SCALING(COMPAT, COUNTS, RATES, N)
%
% Rate scaling a compatibility declaration imposes on its station.
%
% This is what SolverLN carries onto the layer station, and it is NOT
% SN_COMPAT_RATE / SN_COMPAT_PEAK. The denominator is the rate the SAME
% population would obtain under FULL compatibility:
%
%   eta(n) = mu(n) / (peak * min(1, sum_j n(j) / S))
%
% so ETA isolates the effect of the compatibility GRAPH and nothing else. The
% denominator DAMPS BY OCCUPANCY RELATIVE TO THE SERVER COUNT, min(1, N/S),
% because that is precisely what the solver's own multiserver term contributes:
% it applies min(N,S) servers at the average server rate peak/S, so
%
%   min(N,S) * (peak/S) * eta(n) = mu(n)
%
% and the station clears the activated-server rate exactly, at every state.
%
% DAMPING BY min(1, N) INSTEAD -- which this did until 2026-08-28 -- leaves the
% effective law at min(N,S)/S * mu(n), which cancels the REDUNDANCY SPEED-UP the
% activated-server law exists to express: a pool of S servers facing one
% compatible job clears S, not 1, because every one of them works on it and the
% first to finish cancels the rest. Under the old normalization a fully
% compatible pool reduced to the plain multiserver, so the OI machinery did no
% work in the homogeneous case and LDES, which simulates mu(n) directly,
% disagreed with it by that factor.
%
% ETA is therefore ABOVE ONE at low occupancy, which is not a defect: it is the
% speed-up carried by servers that would otherwise be idle.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

total = sum(max(n(:)', 0));
if total <= 0
    eta = 1; % an empty station: nothing to scale
    return
end
nservers = sum(counts(:));
if nservers <= 0
    eta = 1;
    return
end
ref = sn_compat_peak(counts, rates) * min(1, total / nservers);
eta = sn_compat_rate(compat, counts, rates, n) / ref;
end
