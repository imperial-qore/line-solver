function mat = adj_to_mat(adj)
    edges = adj.getEdges();
    mat = [0];
    for i=1:length(edges)
        if ~isempty(edges{i})
            elem = edges{i};
            while ~isempty(elem)
                mat(i,elem.data(1)) = elem.data(2);
                elem = elem.next;
            end
        end
    end
end