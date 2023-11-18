function task = getNodeTask(self,node)
% PROC = GETNODETASK(SELF,NODE)
switch class(node)
    case 'Task'
        task = node;
    case {'Entry','Activity'}
        task = node.parent;
    otherwise
        line_error(mfilename,'Invalid input node');
end
end
