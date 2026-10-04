function [tau, bounded] = pas_saturation_cutoffs(muFun, K, cmax)
% [TAU, BOUNDED] = PAS_SATURATION_CUTOFFS(MUFUN, K, CMAX)
% Per-class saturation cutoffs of an order-independent (OI/PAS) rate mu(c).
%
% Port of the search in Python save_model (_saturation_cutoffs): the smallest
% per-class counts tau with mu(n) = mu(min(n, tau)) on every pair the search
% visits inside the buffer sum(n) <= CMAX (default Inf, no buffer). mu receives
% 1-based class tags, as in oi_rate_table (linemodel_save), and is compared as
% the table holds it (non-finite -> 0). BOUNDED is false when mu keeps changing
% past 256, or saturates no earlier than the buffer. Shared by linemodel_save
% (table extent and oiCutoffs) and ldes_require_finite_pas_buffers (the
% uncapped-buffer refusal, CMAX = Inf); the generic search is saturation_cutoffs.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(cmax)
    cmax = Inf;
end
[tau, bounded] = saturation_cutoffs(@(n) finite_or_zero(muFun(repelem(1:K, n))), K, cmax);
end

function v = finite_or_zero(v)
v = double(v);
v(~isfinite(v)) = 0;
end
