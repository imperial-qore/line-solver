function [I, recurrent] = stronglyconncomp(A)
% [I, RECURRENT] = STRONGLYCONNCOMP(A)
%
% Strongly connected components of the directed graph with adjacency A, where a
% nonzero A(i,j) is an edge i -> j.
%
% I:         SCC index of each node, components numbered by decreasing size
%            (ties broken by discovery order)
% RECURRENT: true for the components with no edge leaving them, i.e. the bottom
%            SCCs (BSCCs). When A is the transition graph of a Markov chain,
%            these are exactly its recurrent classes and the remaining
%            components are transient.
%
% Tarjan's algorithm with an explicit depth-first stack. The recursion is
% unrolled deliberately: this is applied to CTMC/DTMC state spaces whose graphs
% routinely contain paths far longer than MATLAB's recursion limit, and a
% recursive formulation aborts on them.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

N = size(A, 1);
if N == 0
    I = zeros(1, 0);
    recurrent = false(1, 0);
    return
end

A = A ~= 0;
At = A.'; % column v of At lists the out-neighbours of v, cheap to slice

v_idx = zeros(1, N);   % discovery index, 0 while unvisited
v_low = zeros(1, N);   % lowlink
v_stk = false(1, N);   % currently on the Tarjan stack
comp = zeros(1, N);    % SCC index in discovery order
idx = 0;
ncomp = 0;

stk = zeros(1, N);     % Tarjan stack
nstk = 0;

frame_v = zeros(1, N); % DFS stack: node of each open frame
frame_k = zeros(1, N); % out-neighbours of that node already consumed
frame_n = cell(1, N);  % out-neighbour list, fetched once per frame

for root = 1:N
    if v_idx(root) ~= 0
        continue
    end

    idx = idx + 1;
    v_idx(root) = idx;
    v_low(root) = idx;
    nstk = nstk + 1; stk(nstk) = root; v_stk(root) = true;

    nf = 1;
    frame_v(nf) = root;
    frame_k(nf) = 0;
    frame_n{nf} = find(At(:, root))';

    while nf > 0
        v = frame_v(nf);
        nbrs = frame_n{nf};
        if frame_k(nf) < numel(nbrs)
            frame_k(nf) = frame_k(nf) + 1;
            w = nbrs(frame_k(nf));
            if v_idx(w) == 0
                idx = idx + 1;
                v_idx(w) = idx;
                v_low(w) = idx;
                nstk = nstk + 1; stk(nstk) = w; v_stk(w) = true;
                nf = nf + 1;
                frame_v(nf) = w;
                frame_k(nf) = 0;
                frame_n{nf} = find(At(:, w))';
            elseif v_stk(w)
                v_low(v) = min(v_low(v), v_idx(w));
            end
        else
            % v is exhausted: close its component if it is a root, then hand
            % its lowlink back to the parent frame
            if v_low(v) == v_idx(v)
                ncomp = ncomp + 1;
                while true
                    w = stk(nstk); nstk = nstk - 1;
                    v_stk(w) = false;
                    comp(w) = ncomp;
                    if w == v
                        break
                    end
                end
            end
            nf = nf - 1;
            if nf > 0
                p = frame_v(nf);
                v_low(p) = min(v_low(p), v_low(v));
            end
        end
    end
end

% Renumber by decreasing component size; sort is stable, so components of equal
% size keep their discovery order
counts = accumarray(comp(:), 1, [ncomp, 1]);
[~, order] = sort(counts, 'descend');
relabel = zeros(1, ncomp);
relabel(order) = 1:ncomp;
I = relabel(comp);

% A component is recurrent iff no edge leaves it
recurrent = true(1, ncomp);
[er, ec] = find(A);
Icol = I(:);
leaving = Icol(er) ~= Icol(ec);
recurrent(Icol(er(leaving))) = false;

end
