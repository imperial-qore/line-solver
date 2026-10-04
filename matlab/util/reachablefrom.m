function V=reachablefrom(dag,node)
% determines all nodes that can be reached in a dag from the given node
children = find(dag(node,:));
V = [];
for child=children
    V = [V,child,reachablefrom(dag,child)];
end
end