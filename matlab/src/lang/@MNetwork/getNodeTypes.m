function nodeTypes = getNodeTypes(self)
% NODETYPES = GETNODETYPES()

nodeTypes = zeros(self.getNumberOfNodes,1);
for i=1:self.getNumberOfNodes
    switch class(self.nodes{i})
        case 'Cache'
            nodeTypes(i) = NodeType.Cache;
        case 'Logger'
            nodeTypes(i) = NodeType.Logger;
        case 'ClassSwitch'
            nodeTypes(i) = NodeType.ClassSwitch;
        case {'Queue','QueueingStation'}
            if isinf(self.nodes{i}.numberOfServers)
                nodeTypes(i) = NodeType.Delay;
            else                
                nodeTypes(i) = NodeType.Queue;
            end
        case 'Sink'
            nodeTypes(i) = NodeType.Sink;
        case 'Router'
            nodeTypes(i) = NodeType.Router;
        case {'Delay','DelayStation'}
            nodeTypes(i) = NodeType.Delay;
        case 'Fork'
            nodeTypes(i) = NodeType.Fork;
        case 'Join'
            nodeTypes(i) = NodeType.Join;
        case 'Source'
            nodeTypes(i) = NodeType.Source;
        case 'Place'
            nodeTypes(i) = NodeType.Place;
        case 'Transition'
            nodeTypes(i) = NodeType.Transition;
        otherwise
            line_error(mfilename,'Unknown node type.');
    end
end
end