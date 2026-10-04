function h = lqn_dep_layer_handle(f, cols, R, model)
% H = LQN_DEP_LAYER_HANDLE(F, COLS, R, MODEL)
%
% Lift a service-rate dependence handle declared on a LayeredNetwork server to
% the layer station that represents it. F maps the per-operand population vector
% of that server (task j of a Host, entry j of a Task) to a scalar scaling shared
% by every operand or to a per-operand vector. COLS{j} lists the layer classes
% through which operand j occupies the station and R is the number of classes in
% the layer.
%
% Solvers evaluate the handle in two different index spaces: CTMC and the exact
% recursions pass a per-class vector, while the AMVA and NC chain recursions pass
% a per-chain vector (solver_amvald works on Nchain throughout). The returned
% handle therefore reads NUMEL(N) to pick the space, aggregates the operand
% populations in it, and answers a vector of the SAME length, since the caller
% indexes the answer with the same index it passed in. An index that belongs to
% no operand keeps the neutral scaling 1.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

h = @(n) lqn_dep_expand(f, n, cols, R, model);
end

function v = lqn_dep_expand(f, n, cols, R, model)
% V = LQN_DEP_EXPAND(F, N, COLS, R, MODEL) evaluate F on operand populations, spread it back
L = numel(n);
if L ~= R
    cols = lqn_dep_chaincols(cols, model, L);
end
K = length(cols);
nop = zeros(1,K);
for j = 1:K
    if ~isempty(cols{j})
        nop(j) = sum(n(cols{j}));
    end
end
w = f(nop);
v = ones(1,L);
for j = 1:K
    if ~isempty(cols{j})
        v(cols{j}) = w(min(j,numel(w)));
    end
end
end

function ccols = lqn_dep_chaincols(cols, model, L)
% CCOLS = LQN_DEP_CHAINCOLS(COLS, MODEL, L) operand columns in the chain index space
sn = model.getStruct();
ccols = cell(size(cols));
for j = 1:length(cols)
    ch = [];
    for c = cols{j}
        ch = [ch, find(sn.chains(:,c))']; %#ok<AGROW>
    end
    ch = unique(ch);
    ccols{j} = ch(ch <= L);
end
end
