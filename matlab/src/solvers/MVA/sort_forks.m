function [forks, parents] = sort_forks(sn, fjforkmap, fjclassmap, nonfjmodel)
forks = zeros(length(sn.nodetype), max(fjclassmap));
forks(find(sn.nodetype == NodeType.ID_FORK), :) = 1;
parents = zeros(1, length(sn.nodetype));
parents(find(sn.nodetype == NodeType.ID_FORK)) = find(sn.nodetype == NodeType.ID_FORK);
for f=find(sn.nodetype == NodeType.ID_FORK)'
    joinIdx = find(sn.fj(f,:));
    forkauxclasses = find(fjforkmap==f);
    for s=forkauxclasses(:)'
        r = fjclassmap(s);
        discovered_forks = zeros(length(sn.nodetype), 1);
        discovered_forks(find(sn.nodetype == NodeType.ID_FORK)) = 1;
        nested = nested_forks(f, joinIdx, nonfjmodel.getLinkedRoutingMatrix{r,r}, discovered_forks, sn);
        forks(:, r) = forks(:, r) & nested;
        parents(find(nested == 0)) = parents(f);
    end
end
end

function nested = nested_forks(startNode, endNode, conn, forks, sn)
if startNode == endNode
    nested = forks;
    return;
end

for i=find(conn(startNode, :))
    if sn.nodetype(i) == NodeType.ID_FORK
        forks(i) = 0;
    end
    forks = forks & nested_forks(i, endNode, conn, forks, sn);
end
nested = forks;
end