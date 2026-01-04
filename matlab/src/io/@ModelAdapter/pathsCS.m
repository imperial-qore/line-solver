% Finds all the paths coming out of start, computes the total response time along them, and identifies all the stations along these paths
function [ri, stat, RN] = pathsCS(sn, orignodes, P, curNode, endNode, curClass, RN, currentTime, stats)
    if curNode == endNode
        ri = currentTime;
        stat = stats;
        return
    end
    ri = [];
    stat = [];
    currentRn = 0;
    if ~isnan(sn.nodeToStation(curNode))
        currentRn = RN(sn.nodeToStation(curNode), curClass);
        stats(end + 1) = sn.nodeToStation(curNode);
    end
    for transition=find(P((curClass-1)*orignodes+curNode, :))
        nextClass = floor((transition - 1) / orignodes) + 1;
        nextNode = transition-(nextClass-1) * orignodes; % new class
        if sn.nodetype(nextNode) == NodeType.Fork & ~isempty(find(sn.fj(nextNode,:)))
            % Encountered nested fork, compute the max along the parallel paths
            joinIdx = find(sn.fj(nextNode,:));
            if RN(sn.nodeToStation(joinIdx), nextClass) == 0
                [ri1, stat1, RN] = pathsCS(sn, P, nextNode, joinIdx, nextClass, RN, 0, []);
                lambdai = 1./ri1;
                d0 = 0;
                parallel_branches = length(ri1);
                for pow=0:(parallel_branches - 1)
                    current_sum = sum(1./sum(nchoosek(lambdai, pow + 1),2));
                    d0 = d0 + (-1)^pow * current_sum;
                end
                RN(sn.nodeToStation(joinIdx), nextClass) = d0;
                RN(stat1, nextClass) = 0;
            end
            [ri1, stat1, RN] = pathsCS(sn, orignodes, P, joinIdx, endNode, nextClass, RN, currentTime + currentRn, stats);
        else
            [ri1, stat1, RN] = pathsCS(sn, orignodes, P, nextNode, endNode, nextClass, RN, currentTime + currentRn, stats);
        end
        ri = [ri, ri1];
        stat = [stat, stat1];
    end
end
