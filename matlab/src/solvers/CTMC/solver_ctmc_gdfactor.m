function gdFactor = solver_ctmc_gdfactor(sn, stateSpaceAggr, options)
% GDFACTOR = SOLVER_CTMC_GDFACTOR(SN, STATESPACEAGGR, OPTIONS)
%
% Tabulates the globally state-dependent rate scaling phi(n) declared through
% setGlobalDependence, one evaluation per CTMC state. Returns an
% (nstates x nstations*nclasses) matrix whose (s, i + (r-1)*M) entry is the
% scaling of class r at station i in state s, i.e. the column-major flattening
% of the (nstations x nclasses) matrix phi returns.
%
% The handle is evaluated ONCE per state and never per transition: phi may be
% expensive (a bandwidth-sharing allocation solves a convex program per call),
% and within a state it is a constant that multiplies every rate at that state.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;
K = sn.nclasses;
nstates = size(stateSpaceAggr,1);

maxEntries = 3e7;
if isfield(options,'config') && isfield(options.config,'gd_max_entries') && ~isempty(options.config.gd_max_entries)
    maxEntries = options.config.gd_max_entries;
end
if nstates*M*K > maxEntries
    line_error(mfilename, sprintf('The global dependence table would hold %d entries (%d states x %d stations x %d classes), above the budget %d. Lower options.cutoff or raise options.config.gd_max_entries.', nstates*M*K, nstates, M, K, maxEntries));
end

phi = sn.gdscaling;
gdFactor = zeros(nstates, M*K);
for s = 1:nstates
    % columns ((i-1)*K+1):i*K of stateSpaceAggr hold station i, so the reshape
    % is (K x M) and its transpose is the (M x K) matrix phi expects
    n = reshape(stateSpaceAggr(s,1:(M*K)), K, M)';
    v = phi(n);
    if isscalar(v)
        v = v * ones(M,K);
    elseif isequal(size(v),[M,1])
        v = repmat(v(:),1,K);
    elseif ~isequal(size(v),[M,K])
        line_error(mfilename, sprintf('The global dependence handle returned a %s array on a %dx%d population matrix; expected a scalar, an (%d x 1) column or an (%d x %d) matrix.', mat2str(size(v)), M, K, M, M, K));
    end
    if any(~isfinite(v(:))) || any(v(:) < 0)
        line_error(mfilename, sprintf('The global dependence handle returned a non-finite or negative scaling at state %d.', s));
    end
    gdFactor(s,:) = v(:)';
end
end
