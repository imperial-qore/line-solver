function nodeNames = getNodeNames(self)
% NODENAMES = GETNODENAMES()

% The commented block causes issues with Logger nodes
% see e.g., getting_started_ex7

if self.hasStruct && isfield(self.sn,'nodenames')
    nodeNames = self.sn.nodenames;
else
    M = getNumberOfNodes(self);
    nodeNames = string([]); % string array
    for ist=1:M
        nodeNames(ist,1) = self.nodes{ist}.name;
    end
end
end