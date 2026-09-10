function space = fromMarg(sn, ind, ntot, options)
% FROMMARG Generate the state space with a given TOTAL queue length
%
% @brief Creates the state space where a node holds NTOT jobs in total
% @param sn Network structure or Network object
% @param ind Node index
% @param ntot Total number of jobs at the node, all classes summed
% @param options Optional structure with configuration parameters
% @return space Generated state space with the requested total
%
% This is the class-summed counterpart of State.fromMarginal: where
% fromMarginal fixes how many jobs of EACH class the node holds, fromMarg
% fixes only how many jobs it holds ALTOGETHER, and returns the union of
% fromMarginal over every class split of NTOT the node can hold.
%
% A class that is disabled at the station has classcap 0 and is excluded from
% the split enumeration up front rather than after the fact. Asking
% fromMarginal for a job of such a class yields an EMPTY local space, and an
% empty factor is absorbed by the cartesian product instead of annihilating
% it, so the job would silently disappear; see the same trap documented at
% cpp/include/line/lang/qn/state.h.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin<4 %~exist('options','var')
    options.force = false;
end
if isa(sn,'Network')
    sn=sn.getStruct();
end

R = sn.nclasses;
space = [];

if ntot < 0
    return
end

if sn.isstation(ind) && ~isempty(sn.classcap)
    ist = sn.nodeToStation(ind);
    ccap = sn.classcap(ist,:);
else
    ccap = Inf*ones(1,R);
end

% Enumerate the class splits of NTOT that the node can hold. ntot=0 has the
% single empty split, which fromMarginal answers with the per-discipline empty
% state -- do not re-derive that width here.
if ntot == 0
    nset = zeros(1,R);
else
    nset = multichoose(R,ntot);
    nset = nset(all(nset <= repmat(ccap,size(nset,1),1),2),:);
end

% The buffer is RIGHT-aligned, so sub-spaces of different width must be padded
% on the left before they are stacked, exactly as fromMarginal does for the
% reply-block sub-spaces.
subspaces = cell(size(nset,1),1);
maxw = 0;
for j=1:size(nset,1)
    sj = State.fromMarginal(sn, ind, nset(j,:), options);
    if isempty(sj)
        continue
    end
    subspaces{j} = sj;
    maxw = max(maxw, size(sj,2));
end

for j=1:size(nset,1)
    sj = subspaces{j};
    if isempty(sj)
        continue
    end
    if size(sj,2) < maxw
        sj = [zeros(size(sj,1), maxw-size(sj,2)), sj]; %#ok<AGROW>
    end
    space = [space; sj]; %#ok<AGROW>
end

if isempty(space)
    return
end
space = unique(space,'rows'); % do not comment, required to sort empty state as first
space = space(end:-1:1,:); % so that states with jobs in phase 1 comes earlier
end
