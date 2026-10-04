function gdNow = solver_ssa_gdfactor(sn, cur_state)
% GDNOW = SOLVER_SSA_GDFACTOR(SN, CUR_STATE)
%
% Evaluates the globally state-dependent rate scaling phi(n) declared through
% setGlobalDependence on the CURRENT sample-path state, returning the
% (nstations x nclasses) matrix of scalings.
%
% SolverCTMC tabulates phi once per state of the enumerated space; a simulator
% instead has one state at a time, so the same factorization applies with the
% table collapsed to a single row. Within a state phi is a constant, so it
% multiplies the rate of every station service event at that state.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;
R = sn.nclasses;

nir = zeros(M,R);
for ind = 1:sn.nnodes
    if sn.isstation(ind)
        [~, nirow] = State.toMarginal(sn, ind, cur_state{sn.nodeToStateful(ind)});
        nir(sn.nodeToStation(ind), :) = nirow(:)';
    end
end

v = sn.gdscaling(nir);
if isscalar(v)
    v = v * ones(M,R);
elseif isequal(size(v),[M,1])
    v = repmat(v(:),1,R);
elseif ~isequal(size(v),[M,R])
    line_error(mfilename, sprintf('The global dependence handle returned a %s array on a %dx%d population matrix; expected a scalar, an (%d x 1) column or an (%d x %d) matrix.', mat2str(size(v)), M, R, M, M, R));
end
if any(~isfinite(v(:))) || any(v(:) < 0)
    line_error(mfilename, 'The global dependence handle returned a non-finite or negative scaling.');
end
gdNow = v;
end
