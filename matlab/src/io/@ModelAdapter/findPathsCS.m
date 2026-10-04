% Finds the response times along each path leading out of start until
% endNode taking into account class swithces
function ri = findPathsCS(sn, P, curNode, endNode, curClass, toMerge, QN, TN, currentTime, fjclassmap, fjforkmap, nonfjmodel, visited)
    if nargin < 13
        visited = [];
    end
    if curNode == endNode
        % Add the current response time of the parallel branch to the list, but subtract the synchronisation time
        qLen = sum(QN(sn.nodeToStation(curNode), toMerge));
        tput = sum(TN(sn.nodeToStation(curNode), toMerge));
        ri = currentTime - qLen / tput;
        return
    end
    ri = [];
    nfjsn = nonfjmodel.getStruct(false);
    orignodes = length(nfjsn.rtorig{1,1});
    here = (curClass-1)*orignodes+curNode;
    % a path that returns to a (node,class) it already holds is a routing loop,
    % such as the .Aux self-loop of a call with mean above 1, not a new branch
    if any(visited == here)
        return
    end
    visited(end+1) = here;
    for transition=find(P(here, :))
        curMerge = toMerge;
        nextClass = floor((transition - 1) / orignodes) + 1;
        nextNode = transition-(nextClass-1) * orignodes; % new class
        curMerge(1) = nextClass;
        qLen = 0;
        tput = 1;
        if ~isnan(sn.nodeToStation(nextNode))
            qLen = sum(QN(sn.nodeToStation(nextNode), curMerge));
            tput = sum(TN(sn.nodeToStation(nextNode), curMerge));
        end
        if sn.nodetype(nextNode) == NodeType.Fork
            joinIdx = find(sn.fj(nextNode,:));
            s = find(fjforkmap == nextNode & fjclassmap == nextClass);
            paths = ModelAdapter.findPathsCS(sn, P, nextNode, joinIdx, nextClass, [curMerge,s], QN, TN, 0, fjclassmap, fjforkmap, nonfjmodel, visited);
            % The inner join fires on the k-th branch completion, k = the branch
            % count on a standard join and the declared quorum on a PARTIAL one.
            kreq = sn_join_quorum(sn, joinIdx, nextClass, length(paths));
            d0 = fj_ordstat_exp(paths, kreq);
            for cls=[curMerge,s]
                nonfjmodel.nodes{joinIdx}.setService(nonfjmodel.classes{cls}, Exp.fitMean(max(d0 - mean(paths), 0)));
            end
            ri = [ri, ModelAdapter.findPathsCS(sn, P, joinIdx, endNode, nextClass, curMerge, QN, TN, currentTime + d0, fjclassmap, fjforkmap, nonfjmodel, visited)];
        else
            ri = [ri, ModelAdapter.findPathsCS(sn, P, nextNode, endNode, nextClass, curMerge, QN, TN, currentTime + qLen/tput, fjclassmap, fjforkmap, nonfjmodel, visited)];
        end
    end
end