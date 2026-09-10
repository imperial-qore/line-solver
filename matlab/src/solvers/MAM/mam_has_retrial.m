function hasRetrial = mam_has_retrial(sn)
% HASRETRIAL = MAM_HAS_RETRIAL(SN)
%
% True when some station-class pair configures a retrial orbit
% (Queue.setRetrial / setOrbit), read off sn.retrialProc, which refreshStruct
% fills only for a configured delay that is not the Disabled placeholder.
% retrialProc is the ambiguity-free test: retrialType is 0 both for "none"
% and, in MATLAB, for an exponential delay (see CLAUDE.md).
%
% Lifted beside MAM_HAS_RENEGING_PATIENCE so that SolverMAM.supportsModelMethod
% can ask whether 'default' and 'dec.source' have an orbit to route to the
% retrial analyzer: QSYS_IS_RETRIAL answers about the SHAPE that analyzer
% needs, not about whether an orbit exists, and a model with an orbit outside
% that shape fell through to solver_mam_basic, which reads no sn.retrial*
% field and answered with the refused jobs simply lost.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

hasRetrial = isfield(sn, 'retrialProc') && ~isempty(sn.retrialProc) ...
    && any(~cellfun(@isempty, sn.retrialProc(:)));
end
