% Finds the response times along each path leading out of start until endNode
function ri = findPaths(sn, P, start, endNode, r, toMerge, QN, TN, currentTime, fjclassmap, fjforkmap, nonfjmodel)
    if start == endNode
        % Add the current response time of the parallel branch to the list, but subtract the synchronisation time
        qLen = sum(QN(sn.nodeToStation(start), toMerge));
        tput = sum(TN(sn.nodeToStation(start), toMerge));
        ri = currentTime - qLen / tput;
        return
    end
    ri = [];
    for i=find(P(start, :))
        qLen = 0;
        tput = 1;
        if ~isnan(sn.nodeToStation(i))
            qLen = sum(QN(sn.nodeToStation(i), toMerge));
            tput = sum(TN(sn.nodeToStation(i), toMerge));
        end
        if sn.nodetype(i) == NodeType.Fork
            joinIdx = find(sn.fj(i,:));
            s = find(fjforkmap == i & fjclassmap == r);
            paths = ModelAdapter.findPaths(sn, P, i, joinIdx, r, [toMerge,s], QN, TN, 0, fjclassmap, fjforkmap, nonfjmodel);
            lambdai = 1./paths;
            d0 = 0;
            parallel_branches = length(paths);
            for pow=0:(parallel_branches - 1)
                current_sum = sum(1./sum(nchoosek(lambdai, pow + 1),2));
                d0 = d0 + (-1)^pow * current_sum;
            end
            for cls=[toMerge,s]
                nonfjmodel.nodes{joinIdx}.setService(nonfjmodel.classes{cls}, Exp.fitMean(d0 - mean(paths)));
            end
            ri = [ri, ModelAdapter.findPaths(sn, P, joinIdx, endNode, r, toMerge, QN, TN, currentTime + d0, fjclassmap, fjforkmap, nonfjmodel)];
        else
            ri = [ri, ModelAdapter.findPaths(sn, P, i, endNode, r, toMerge, QN, TN, currentTime + qLen/tput, fjclassmap, fjforkmap, nonfjmodel)];
        end
    end
end