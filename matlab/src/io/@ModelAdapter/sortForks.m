function [forks, parents] = sortForks(sn, fjforkmap, fjclassmap, nonfjmodel)
forks = zeros(length(sn.nodetype), max(fjclassmap));
forks(find(sn.nodetype == NodeType.Fork), :) = 1;
parents = zeros(1, length(sn.nodetype));
parents(find(sn.nodetype == NodeType.Fork)) = find(sn.nodetype == NodeType.Fork);
for f=find(sn.nodetype == NodeType.Fork)'
    joinIdx = find(sn.fj(f,:));
    forkauxclasses = find(fjforkmap==f);
    for s=forkauxclasses(:)'
        r = fjclassmap(s);
        discovered_forks = zeros(length(sn.nodetype), 1);
        discovered_forks(find(sn.nodetype == NodeType.Fork)) = 1;
        nested = nestedForks(f, joinIdx, nonfjmodel.getLinkedRoutingMatrix{r,r}, discovered_forks, sn);
        forks(:, r) = forks(:, r) & nested;
        parents(find(nested == 0)) = parents(f);
    end
end
end

function nested = nestedForks(startNode, endNode, conn, forks, sn)
if startNode == endNode
    nested = forks;
    return;
end

for i=find(conn(startNode, :))
    if sn.nodetype(i) == NodeType.Fork
        forks(i) = 0;
    end
    forks = forks & nestedForks(i, endNode, conn, forks, sn);
end
nested = forks;
end