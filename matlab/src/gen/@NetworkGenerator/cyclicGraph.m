function graph = cyclicGraph(numVertices)
% This is an alternative topologyFcn for NetworkGenerator. It takes in
% an integer number of vertices and returns a MATLAB digraph object
% representing a cyclic directed graph.

    graph = digraph;
    graph = addnode(graph, numVertices);
    if numVertices < 1
        error('cyclicGraph:nonPositiveVertices', 'numVertices must be positive');
    elseif numVertices == 1
        graph = addedge(graph, 1, 1);
    else
        for i = 1 : numVertices - 1
            graph = addedge(graph, i, i + 1);
        end
        graph = addedge(graph, numVertices, 1);
    end
end