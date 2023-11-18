function proc = getNodeHost(self,node)
% PROC = GETNODEHOST(SELF,NODE)
switch class(node)
    case {'Processor','Host'}
        proc = node;
    case 'Task'
        proc = node.parent;
    case {'Entry','Activity'}
        proc = node.parent.parent;
    otherwise
        line_error(mfilename,'Invalid input node');
end
end
