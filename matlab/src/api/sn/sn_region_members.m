function mask = sn_region_members(sn, f, Rmat, memvec)
% MASK = SN_REGION_MEMBERS(SN, F, RMAT, MEMVEC) 1xM logical membership mask of
% finite capacity region F.
%
% Membership is read from sn.regionmembers{f}, which refreshRegions records
% directly from the region's node list. It cannot be derived from sn.region{f}:
% -1 there means "unbounded", which is indistinguishable from "not a member", so
% a region constrained only by regionlincon (or only by a memory budget) reads as
% empty and is silently ignored.
%
% RMAT and MEMVEC provide the legacy derivation, used only for an sn built before
% regionmembers existed (for instance one deserialised from an older model file).
% That derivation carries the ambiguity above and is not equivalent.

if isfield(sn,'regionmembers') && numel(sn.regionmembers) >= f && ~isempty(sn.regionmembers{f})
    mask = logical(sn.regionmembers{f}(:))';
    return
end
mask = (any(Rmat ~= -1, 2) | memvec(:) ~= -1)';
end
