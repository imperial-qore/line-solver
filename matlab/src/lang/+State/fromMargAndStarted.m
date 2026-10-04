function space = fromMargAndStarted(sn, ind, ntot, stot, options)
% FROMMARGANDSTARTED Generate the states with a given TOTAL queue length and
% a given TOTAL number of started jobs
%
% @brief Class-summed counterpart of State.fromMarginalAndStarted
% @param sn Network structure or Network object
% @param ind Node index
% @param ntot Total number of jobs at the node, all classes summed
% @param stot Total number of jobs that have started service
% @param options Optional structure with configuration parameters
% @return space Generated state space
%
% Where fromMarginalAndStarted takes one per-class vector N and one per-class
% vector S and builds ONE row, fromMargAndStarted takes only the two totals and
% returns the union of that row over every (N,S) pair consistent with them:
% sum(N)=NTOT, sum(S)=STOT, and S <= N elementwise.
%
% Classes disabled at the station are excluded from the enumeration through
% classcap, for the reason documented in State.fromMarg.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin<5 %~exist('options','var')
    options.force = true;
end
if isa(sn,'Network')
    sn=sn.getStruct();
end

R = sn.nclasses;
space = [];

if ntot < 0 || stot < 0 || stot > ntot
    return
end

if sn.isstation(ind) && ~isempty(sn.classcap)
    ist = sn.nodeToStation(ind);
    ccap = sn.classcap(ist,:);
else
    ccap = Inf*ones(1,R);
end

if ntot == 0
    nset = zeros(1,R);
else
    nset = multichoose(R,ntot);
    nset = nset(all(nset <= repmat(ccap,size(nset,1),1),2),:);
end

rows = {};
maxw = 0;
for j=1:size(nset,1)
    nj = nset(j,:);
    % s must be drawn from the jobs actually present, which is what
    % multichoosecon expresses; S=0 is its one uncovered base case.
    if stot == 0
        sset = zeros(1,R);
    else
        sset = multichoosecon(nj,stot);
    end
    for k=1:size(sset,1)
        sjk = State.fromMarginalAndStarted(sn, ind, nj, sset(k,:), options);
        if isempty(sjk)
            continue
        end
        rows{end+1} = sjk; %#ok<AGROW>
        maxw = max(maxw, size(sjk,2));
    end
end

for j=1:numel(rows)
    sj = rows{j};
    if size(sj,2) < maxw
        sj = [zeros(size(sj,1), maxw-size(sj,2)), sj]; %#ok<AGROW>
    end
    space = [space; sj]; %#ok<AGROW>
end

if isempty(space)
    return
end
space = unique(space,'rows');
space = space(end:-1:1,:);
end
