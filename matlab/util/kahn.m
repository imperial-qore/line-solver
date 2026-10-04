function [order, acyclic] = kahn(adj)
% Kahn's algorithm – topological sort + cycle flag
n        = size(adj,1);
indegree = sum(adj,1);
order    = zeros(1,n);   % pre-allocate
idx      = 1;
queue    = find(indegree==0);

while ~isempty(queue)
    i         = queue(1);      queue(1) = [];
    order(idx)= i;             idx      = idx + 1;

    for j = find(adj(i,:))
        indegree(j) = indegree(j) - 1;
        if indegree(j) == 0
            queue(end+1) = j;  %#ok<AGROW>
        end
    end
end

acyclic = (idx-1) == n;   % true  ⇔ no cycle
end
