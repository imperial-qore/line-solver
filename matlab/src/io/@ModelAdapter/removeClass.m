function newmodel = removeClass(model, jobclass)
% NEWMODEL = REMOVECLASS(MODEL, JOBCLASS)
%
% Non-mutating variant of Network.removeClass: return a copy of MODEL in
% which JOBCLASS has been removed, leaving MODEL itself untouched. Use this
% when the original model must stay solvable (e.g. ablation studies, or a
% per-chain decomposition that peels one class at a time).
%
% The removal logic lives in @MNetwork/removeClass.m and is delegated to
% here rather than duplicated, so that the two entry points cannot drift.
% The lookup on the copy is by class name (see @MNetwork/removeClass.m), so
% passing the JOBCLASS object of the original model is correct.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

newmodel = model.copy();
newmodel.removeClass(jobclass);
end
