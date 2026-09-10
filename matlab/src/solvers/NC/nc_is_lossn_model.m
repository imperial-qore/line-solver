function [isLossn, hasShape] = nc_is_lossn_model(sn)
% [ISLOSSN, HASSHAPE] = NC_IS_LOSSN_MODEL(SN)
%
% Is this the loss network SOLVER_NC_LOSSN_ANALYZER solves, i.e. an OPEN model
% with a single finite capacity region whose only member is an infinite server,
% and whose admission rule DROPS every class?
%
% HASSHAPE reports the topology alone, without the DROP demand: a region of
% that shape under WAITQ (or any blocking rule) holds the arrival back instead
% of discarding it, which keeps the job in the region while it waits and is a
% queueing phenomenon the Erlang loss model has no state for. The two answers
% are returned separately because they have different remedies -- switching the
% rule to DROP makes the first solvable here, while the second needs a solver
% that carries the region as state.
%
% Factored out of @SolverNC/runAnalyzer, which decided the same question inline
% and could therefore not be asked it by the support gate; see
% NC_METHOD_REFUSAL.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

isLossn = false;
hasShape = false;

if ~isfield(sn,'nregions') || isempty(sn.nregions) || sn.nregions ~= 1
    return
end
if sn_has_closed_classes(sn)
    return
end
if ~isfield(sn,'region') || isempty(sn.region) || isempty(sn.region{1})
    return
end
regionMatrix = sn.region{1};
% Stations in the region: a non-negative per-class limit, or a non-negative
% aggregate limit in the last column. Same test runAnalyzer applied inline.
stationsInFCR = find(any(regionMatrix(:,1:end-1) >= 0, 2) | regionMatrix(:,end) >= 0);
if numel(stationsInFCR) ~= 1
    return
end
if ~isinf(sn.nservers(stationsInFCR(1)))
    return
end
hasShape = true;

if ~isfield(sn,'regionrule') || isempty(sn.regionrule)
    return
end
% ALL classes, not merely one: a region that discards one class and holds
% another back is a mixed system whose blocked class occupies the region while
% it waits, so the per-class loss probabilities the Erlang fixed point returns
% would not be the ones the model implies.
isLossn = all(sn.regionrule(1,:) == DropStrategy.DROP);
end
