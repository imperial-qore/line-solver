function graph = randGraph(numVertices)
% This function takes as input an integer number of vertices
% and returns a MATLAB digraph object representing a randomly
% generated, strongly connected graph

    if numVertices <= 0
        error('randGraph:positiveVerticesRequired', 'Positive vertices required');
    end
    
    if numVertices == 1
        graph = digraph;
        graph = addnode(graph, 1);
        graph = addedge(graph, 1, 1);
        return
    end

    % State variables needed for the DFS-based strongConnect function
    globalStartTime = 1;
    startTime = ones(1, numVertices) * -1;
    lowestLink = startTime;
    invStartTime = startTime;
     
    tree = randSpanningTree(numVertices);
    graph = strongConnect(tree, 1);
    graph = reordernodes(graph, randperm(numVertices)); % To generate arbitrary permutations
    
    % Nested function to avoid using global state variables
    function g = strongConnect(g, v)
        startTime(v) = globalStartTime; % Order of observation in DFS sequence
        invStartTime(startTime(v)) = v; % For inverse lookup when adding edges
        lowestLink(v) = startTime(v); % Initially can only reach itself
        globalStartTime = globalStartTime + 1;

        [~, vertexIDs] = outedges(g, v);
        for i = 1 : length(vertexIDs)
            w = vertexIDs(i);
            if startTime(w) == -1 % Vertex w has not been visited yet
                g = strongConnect(g, w);
                lowestLink(v) = min(lowestLink(v), lowestLink(w));
            else 
                lowestLink(v) = min(lowestLink(v), startTime(w));
            end
        end

        % If vertex v is the root of a SCC but not the entire graph's root
        if lowestLink(v) == startTime(v) && startTime(v) > 1
            descendantST = randi([startTime(v) globalStartTime - 1]);
            ancestorST = randi([1 startTime(v) - 1]);
            g = addedge(g, invStartTime(descendantST), invStartTime(ancestorST));
            lowestLink(v) = ancestorST;
        end 
    end
end

% Generates a random spanning tree with the specified number of vertices
function tree = randSpanningTree(numVertices)
    tree = digraph;
    tree = addnode(tree, numVertices);
    
    % Note: This always joins vertex 1 to 2, we will relabel in randGraph
    for i = 2 : numVertices
        tree = addedge(tree, randi(i - 1), i);
    end
end