function [I, recurrent] = stronglyconncomp(A)
% find strongly connected components in directed graph A
N = size(A, 1);
idx = 0;
SCC = {};
stk = [];

v_idx = zeros(1, N);
v_low = zeros(1, N);
v_stk = false(1, N);

for i = 1:N
    if v_idx(i) == 0
        [v_idx, v_low, v_stk, SCC, stk, idx] = stronglyconncomp_aux(i, A, v_idx, v_low, v_stk, SCC, stk, idx);
    end
end

% Sort SCCs by size
[~, q] = sort(cellfun(@length, SCC), 'descend');
SCC = SCC(q);
I = zeros(1, N);
for j = 1:length(SCC)
    I(SCC{j}) = j;
end

recurrent=false(1, length(SCC));
for j = 1:length(SCC)
    scc = SCC{j};
    is_recurrent = true;

    % Check if any node in the SCC has outgoing edges to nodes outside the SCC
    for node = scc
        out_edges = find(A(node, :));
        if any(~ismember(out_edges, scc))
            is_recurrent = false;
            break;
        end
    end

    if is_recurrent
        recurrent(j) = true;
    end
end
end

function [v_idx, v_low, v_stk, SCC, stk, idx] = stronglyconncomp_aux(i, e, v_idx, v_low, v_stk, SCC, stk, idx)
idx = idx + 1;
v_idx(i) = idx;
v_low(i) = idx;
stk = [i, stk];
v_stk(i) = true;

out_edges = find(e(:, i));

for j = out_edges'
    if v_idx(j) == 0
        [v_idx, v_low, v_stk, SCC, stk, idx] = stronglyconncomp_aux(j, e, v_idx, v_low, v_stk, SCC, stk, idx);
        v_low(i) = min(v_low(i), v_low(j));
    elseif v_stk(j)
        v_low(i) = min(v_low(i), v_idx(j));
    end
end

if v_low(i) == v_idx(i)
    scc = stk(1:find(stk == i));
    stk = stk(length(scc) + 1:end);
    v_stk(scc) = false;
    SCC = [SCC, {scc}];
end
end
