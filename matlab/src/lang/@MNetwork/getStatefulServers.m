function S = getStatefulServers(self)
% S = GETSTATEFULSERVERS()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

for i=1:getNumberOfStatefulNodes(self)    
    if isa(self.nodes{i}, 'Station')
        S(i,1) = self.nodes{i}.numberOfServers;
    elseif isa(self.nodes{i}, 'StatefulNode')
        if class(self.nodes{i}) == 'Transition'
            S(i,1) = self.nodes{i}.numberOfServers;
        else
            S(i,1) = 0;
        end
    end
end
end
