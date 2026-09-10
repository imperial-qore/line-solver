%{ @file cache_gamma.m
 %  @brief Access factors of a multi-list cache whose lists form a general graph
 %
 %  @author LINE Development Team
%}

%{
 % @brief Access factors of a multi-list cache over a general access graph
 %
 % @details
 % Companion of cache_gamma_lp, which requires the access structure to be a
 % tree and walks the unique parent relation. Here the structure is only
 % required to be reachable: the path to list j is the breadth-first shortest
 % path in the access graph of item i, so a list with several parents is
 % admissible and the first shortest path found in node order is the one taken.
 % Along that path,
 %
 %   gamma(i,j) = (sum_v lambda(v,i,1)) prod_{edges (a,b)} sum_v lambda(v,i,a) R{v,i}(a,b)
 %
 % THREE DIVERGENCES FROM cache_gamma_lp, all of which change the number, so
 % the two are not substitutes:
 %   - the destination of column j is node j, not node 1+j: column 1 therefore
 %     has the trivial one-node path and carries no edge factor, where the tree
 %     version carries the miss-to-first-list edge;
 %   - the leading factor is the aggregate miss-node request rate, whereas the
 %     tree version starts the product at one;
 %   - each edge factor reads lambda at the source node a alone, whereas the
 %     tree version sums lambda(v,i,t) over every t <= a.
 % See _kb/09-ldes-and-cache.md. Use cache_gamma_lp for the access factors a
 % cache solver consumes.
 %
 % The adjacency is read from user 1 only; the per-user rates are then summed
 % along that one path, so a model whose users route an item differently is
 % analysed on the first user's graph. An unreachable list gives gamma(i,j)=0.
 %
 % @par Syntax:
 % @code
 % [gamma, u, n, h] = cache_gamma(lambda, R)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lambda<td>(u,n,h+1) request rate of user v for item i while at node t
 % <tr><td>R<td>(u,n) cell of (h+1)x(h+1) routing matrices
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>(n,h) access factors
 % <tr><td>u<td>Number of users
 % <tr><td>n<td>Number of items
 % <tr><td>h<td>Number of lists
 % </table>
%}
function [gamma,u,n,h]=cache_gamma(lambda,R)
u=size(lambda,1); % number of users
n=size(lambda,2); % number of items
h=size(lambda,3)-1; % number of lists

gamma=zeros(n,h);
for i = 1:n
    graph = R{1,i}; % the adjacency comes from the first user alone
    for j = 1:h
        Pij = bfs_path(graph, 1, j);
        if isempty(Pij)
            gamma(i,j) = 0;
        else
            g = 0;
            for v = 1:u
                g = g + lambda(v,i,1);
            end
            for li = 2:length(Pij)
                a = Pij(li-1);
                b = Pij(li);
                y = 0;
                for v = 1:u
                    y = y + lambda(v,i,a) * R{v,i}(a,b);
                end
                g = g * y;
            end
            gamma(i,j) = g;
        end
    end
end

end

function path = bfs_path(A, src, dst)
% breadth-first shortest path over the strictly positive entries of A
nn = size(A,1);
path = [];
if src > nn || dst > nn || src < 1 || dst < 1
    return
end
visited = false(1,nn);
parent = zeros(1,nn);
queue = src;
visited(src) = true;
while ~isempty(queue)
    current = queue(1);
    queue(1) = [];
    if current == dst
        path = dst;
        node = parent(dst);
        while node ~= 0
            path = [node, path]; %#ok<AGROW>
            node = parent(node);
        end
        return
    end
    for next = 1:nn
        if ~visited(next) && A(current,next) > 0
            visited(next) = true;
            parent(next) = current;
            queue(end+1) = next; %#ok<AGROW>
        end
    end
end
end
