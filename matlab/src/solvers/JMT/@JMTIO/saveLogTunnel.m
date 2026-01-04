function [simDoc, section] = saveLogTunnel(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVELOGTUNNEL(SIMDOC, SECTION, NODEIDX)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.getStruct;
loggerNodesCP = {'java.lang.String','java.lang.String'};
for i=3:9 loggerNodesCP{i} = 'java.lang.Boolean'; end
loggerNodesCP{10} = 'java.lang.Integer';

loggerNodesNames = {'logfileName','logfilePath','logExecTimestamp', ...
    'logLoggerName','logTimeStamp','logJobID', ...
    'logJobClass','logTimeSameClass','logTimeAnyClass', ...
    'numClasses'};
numOfClasses = sn.nclasses;

% logger specific path does not work in JMT at the moment
if ~strcmpi(sn.nodeparam{ind}.filePath(end),filesep)
    currentNode.filePath = [sn.nodeparam{ind}.filePath, filesep];
end

loggerNodesValues = {sn.nodeparam{ind}.fileName, sn.nodeparam{ind}.filePath, ...
    sn.nodeparam{ind}.startTime,sn.nodeparam{ind}.loggerName, ...
    sn.nodeparam{ind}.timestamp,sn.nodeparam{ind}.jobID, ...
    sn.nodeparam{ind}.jobClass,sn.nodeparam{ind}.timeSameClass, ...
    sn.nodeparam{ind}.timeAnyClass,int2str(numOfClasses)};

for j=1:length(loggerNodesValues)
    loggerNode = simDoc.createElement('parameter');
    loggerNode.setAttribute('classPath', loggerNodesCP{j});
    loggerNode.setAttribute('name', loggerNodesNames{j});
    valueNode = simDoc.createElement('value');
    valueNode.appendChild(simDoc.createTextNode(loggerNodesValues{j}));
    loggerNode.appendChild(valueNode);
    section.appendChild(loggerNode);
end
end
