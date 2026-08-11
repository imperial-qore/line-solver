function tf = is_dag(adj)
% Checks if the input adjacency matrix is a DAG using Kahn's algorithm
tol  = GlobalConstants.FineTol;
bin  = double(adj > tol);   % threshold & binarise
[~, tf] = kahn(bin);        % tf == acyclic flag
end

