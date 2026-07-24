function nodes = resetNetwork(self, deleteCSNodes) % resets network topology
% NODES = RESETNETWORK(DELETECSNODES) % RESETS NETWORK TOPOLOGY

M = self.getNumberOfStations;
if nargin<2 %~exist('deleteNodes','var')
    deleteCSNodes = true;
end

% remove class switch nodes
if deleteCSNodes
    oldNodes = self.nodes;
    self.nodes = {};
    for notCS = find(~cellisa(oldNodes,'ClassSwitch'))'
        self.nodes{end+1,1} = oldNodes{notCS};
    end
end

for i = 1:M
    self.stations{i}.output.initDispatcherJobClasses(self.classes);
end

self.handles = {};
self.connections = [];
nodes = self.getNodes;
end
