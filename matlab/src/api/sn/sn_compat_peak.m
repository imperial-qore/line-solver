function peak = sn_compat_peak(counts, rates)
% PEAK = SN_COMPAT_PEAK(COUNTS, RATES)
%
% Rate a compatibility declaration clears with every pool active,
% sum_t counts(t)*rates(t).
%
% Utilization at a rate-scaled station is reported as U = T*S/peak, and the peak
% is a property of the DECLARATION rather than of a state, so it is computed
% once and handed to the solver beside the rate handle rather than recovered
% from SN_COMPAT_RATE at a guessed state.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if numel(rates) ~= numel(counts)
    line_error(mfilename, 'one rate per pool is required');
end
peak = sum(counts(:)' .* rates(:)');
end
