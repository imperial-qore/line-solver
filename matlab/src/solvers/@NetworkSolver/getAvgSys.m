function [CNchain,XNchain] = getAvgSys(self,R,T)
% [CNCHAIN,XNCHAIN] = GETAVGSYS(SELF,R,T)

% Return average system metrics at steady state
%
% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

sn = self.model.getStruct();
if nargin < 3
    R = self.getAvgRespTHandles;
    T = self.getAvgTputHandles;
end
[~,~,RN,TN] = self.getAvg([],[],R,T,[],[]);
RN(sn.nodeToStation(find(sn.nodetype == NodeType.ID_JOIN)), :) = 0;
for f=find(sn.nodetype == NodeType.ID_FORK)'
    % find join associated to fork f
    joinIdx = find(sn.fj(f,:));
    if isempty(joinIdx)
        % No join nodes for this fork, no synchronisation delay
        continue;
    end
    for c=1:sn.nchains
        inchain = sn.inchain{c};
        for r=inchain
            if RN(sn.nodeToStation(joinIdx), r) == 0
                % Find the parallel paths coming out of the fork
                [ri, stat, RN] = fj_paths(sn, self.model.getLinkedRoutingMatrix{r,r}, f, joinIdx, r, RN, 0, []);
                lambdai = 1./ri;
                d0 = 0;
                parallel_branches = length(ri);
                for pow=0:(parallel_branches - 1)
                    current_sum = sum(1./sum(nchoosek(lambdai, pow + 1),2));
                    d0 = d0 + (-1)^pow * current_sum;
                end
                RN(sn.nodeToStation(joinIdx), r) = d0;
                RN(stat, r) = 0;
            end
        end
    end
end

refstats = sn.refstat;
completes = true(1,sn.nclasses);
for r=1:sn.nclasses
    completes(r) = T{refstats(r),r}.class.completes;
end

%if any(isinf(sn.njobs')) % if the model has any open class
% TODO: this could be optimised by computing the statistics
% only for open chains

% compute chain visits
alpha = zeros(sn.nstations,sn.nclasses);
CNclass = zeros(1,sn.nclasses);
for c=1:sn.nchains
    inchain = sn.inchain{c};
    for r=inchain
        CNclass(r)=0;
        for i=1:sn.nstations
            if ~isempty(RN) && ~(isinf(sn.njobs(r)) && i==sn.refstat(r)) % not empty and not source
                CNclass(r) = CNclass(r) + sn.visits{c}(sn.stationToStateful(i),r)*RN(i,r)/sn.visits{c}(sn.stationToStateful(sn.refstat(r)),r);
            end
        end
    end
end

for c=1:sn.nchains
    inchain = sn.inchain{c};
    completingclasses = sn.chains(c,:) & completes;
    for i=1:sn.nstations
        if sn.refclass(c)>0
            for k=intersect(find(sn.refclass), inchain) % for all classes within the chain (a class belongs to a single chain, the reference station must be identical for all classes within a chain )
                alpha(i,k) = alpha(i,k) + sn.visits{c}(sn.stationToStateful(i),k)/sum(sn.visits{c}(sn.stationToStateful(sn.refstat(k)),completingclasses));
            end
        else
            for k=inchain % for all classes within the chain (a class belongs to a single chain, the reference station must be identical for all classes within a chain )
                alpha(i,k) = alpha(i,k) + sn.visits{c}(i,k)/sum(sn.visits{c}(sn.stationToStateful(sn.refstat(k)),completingclasses));
            end
        end
    end
end
alpha(~isfinite(alpha))=0;
%end

% compute average chain metrics
CNchain = zeros(1,sn.nchains);
XNchain = zeros(1,sn.nchains);
for c=1:sn.nchains
    inchain = sn.inchain{c};
    completingclasses = find(sn.chains(c,:) & completes);
    if ~isempty(TN)
        XNchain(c) = 0;
        % all classes in same chain must share the same refstation, so we use the first one
        ref = refstats(inchain(1));
        % we now compute the incoming system throughput to the
        % reference station from completing classes
        for i=1:sn.nstations
            for r=completingclasses(:)'
                if any(intersect(find(sn.refclass), inchain))
                    for s=intersect(find(sn.refclass), inchain)
                        if ~isnan(TN(i,r))
                            XNchain(c) = XNchain(c) + sn.rt((i-1)*sn.nclasses + r, (ref-1)*sn.nclasses + s )*TN(i,r);
                        end
                    end
                else
                    for s=inchain(:)'
                        if ~isnan(TN(i,r))
                            XNchain(c) = XNchain(c) + sn.rt((i-1)*sn.nclasses + r, (ref-1)*sn.nclasses + s )*TN(i,r);
                        end
                    end
                end
            end
        end
    end

    % If this is a closed chain we simply apply Little's law
    nJobsChain = sum(sn.njobs(find(sn.chains(c,:)))); %#ok<FNDSB>
    if isinf(nJobsChain)
        if length(inchain) ~= length(completingclasses)
            line_error(mfilename,'Edge-based chain definition not yet supported for open queueing networks.');
            %else
            % we use nan sum to disregard response at stations where
            % the class is not defined
            %    CNchain(c) = sumfinite(alpha(refstats(inchain(1)),inchain).*CNclass(inchain));
        end
    end
    CNchain(c) = sumfinite(alpha(refstats(inchain(1)),inchain).*CNclass(inchain));  
end
self.result.XN = XNchain;
self.result.CN = CNchain;
end

% Finds all the paths coming out of start, computes the total response time along them, and identifies all the stations along these paths
function [ri, stat, RN] = fj_paths(sn, P, start, endNode, r, RN, currentTime, stats)
    if start == endNode
        ri = currentTime;
        stat = stats;
        return
    end
    ri = [];
    stat = [];
    currentRn = 0;
    if ~isnan(sn.nodeToStation(start))
        currentRn = RN(sn.nodeToStation(start), r);
        stats(end + 1) = sn.nodeToStation(start);
    end
    for i=find(P(start, :))
        if sn.nodetype(i) == NodeType.ID_FORK & ~isempty(find(sn.fj(i,:)))
            % Encountered nested fork, compute the max along the parallel paths
            joinIdx = find(sn.fj(i,:));
            if RN(sn.nodeToStation(joinIdx), r) == 0
                [ri1, stat1, RN] = fj_paths(sn, P, i, joinIdx, r, RN, 0, []);
                lambdai = 1./ri1;
                d0 = 0;
                parallel_branches = length(ri1);
                for pow=0:(parallel_branches - 1)
                    current_sum = sum(1./sum(nchoosek(lambdai, pow + 1),2));
                    d0 = d0 + (-1)^pow * current_sum;
                end
                RN(sn.nodeToStation(joinIdx), r) = d0;
                RN(stat1, r) = 0;
            end
            [ri1, stat1, RN] = fj_paths(sn, P, joinIdx, endNode, r, RN, currentTime + currentRn, stats);
        else
            [ri1, stat1, RN] = fj_paths(sn, P, i, endNode, r, RN, currentTime + currentRn, stats);
        end
        ri = [ri, ri1];
        stat = [stat, stat1];
    end
end
