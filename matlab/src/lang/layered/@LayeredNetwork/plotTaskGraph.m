function plotTaskGraph(self, method)
% PLOTTASKGRAPH(SELF, METHOD)
%
% METHOD: nodes, names or ids

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin<2 %~exist('useNodes','var')
    method = 'nodes';
end

lqn = self.getStruct;
T = zeros(lqn.nhosts+lqn.ntasks);
for h=1:lqn.nhosts
    hidx = h;
    for tidx=lqn.tasksof{h}
        T(tidx,hidx) = 1;
    end
end
for t=1:lqn.ntasks
    tidx = lqn.tshift + t;
    [calling_idx, called_entries] = find(lqn.iscaller(:, lqn.entriesof{tidx})); %#ok<ASGLU>
    callers = intersect(lqn.tshift+(1:lqn.ntasks), unique(calling_idx)');
    T(callers,tidx) = 1;
end
figure;
switch method
    case 'nodes'
        plot(digraph(T),'Layout','layered','NodeLabel',{lqn.hashnames{1:(lqn.nhosts+lqn.ntasks)}});
    case 'names'
        plot(digraph(T),'Layout','layered','NodeLabel',{lqn.names{1:(lqn.nhosts+lqn.ntasks)}});
    case 'ids'
        plot(digraph(T),'Layout','layered');
end
title('Task graph');
end
