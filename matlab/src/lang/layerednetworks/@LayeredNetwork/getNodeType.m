function type = getNodeType(self,node)
% TYPE = GETNODETYPE(SELF,NODE)
switch class(node)
    case {'Processor','Host'}
        type = LayeredNetworkElement.HOST;
    case 'Task'
        type = LayeredNetworkElement.TASK;
    case 'Entry'
        type = LayeredNetworkElement.ENTRY;
    case 'Activity'
        type = LayeredNetworkElement.ACTIVITY;
    otherwise
        line_error(mfilename,'Invalid input node');
end
end
