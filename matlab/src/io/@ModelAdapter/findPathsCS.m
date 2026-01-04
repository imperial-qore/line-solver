% Finds the response times along each path leading out of start until
% endNode taking into account class swithces
function ri = findPathsCS(sn, P, curNode, endNode, curClass, toMerge, QN, TN, currentTime, fjclassmap, fjforkmap, nonfjmodel)
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
    for transition=find(P((curClass-1)*orignodes+curNode, :))
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
            paths = ModelAdapter.findPathsCS(sn, P, nextNode, joinIdx, nextClass, [curMerge,s], QN, TN, 0, fjclassmap, fjforkmap, nonfjmodel);
            lambdai = 1./paths;
            d0 = 0;
            parallel_branches = length(paths);
            for pow=0:(parallel_branches - 1)
                current_sum = sum(1./sum(nchoosek(lambdai, pow + 1),2));
                d0 = d0 + (-1)^pow * current_sum;
            end
            for cls=[curMerge,s]
                nonfjmodel.nodes{joinIdx}.setService(nonfjmodel.classes{cls}, Exp.fitMean((d0 - mean(paths))));
            end
            ri = [ri, ModelAdapter.findPathsCS(sn, P, joinIdx, endNode, nextClass, curMerge, QN, TN, currentTime + d0, fjclassmap, fjforkmap, nonfjmodel)];
        else
            ri = [ri, ModelAdapter.findPathsCS(sn, P, nextNode, endNode, nextClass, curMerge, QN, TN, currentTime + qLen/tput, fjclassmap, fjforkmap, nonfjmodel)];
        end
    end
end