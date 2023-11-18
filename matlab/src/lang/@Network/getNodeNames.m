function nodeNames = getNodeNames(self)
% NODENAMES = GETNODENAMES()

% The commented block causes issues with Logger nodes
% see e.g., getting_started_ex7

if self.hasStruct && isfield(self.sn,'nodenames')
    nodeNames = self.sn.nodenames;
else
    M = getNumberOfNodes(self);
    nodeNames = string([]); % string array
    for i=1:M
        nodeNames(i,1) = self.nodes{i}.name;
    end
end
end