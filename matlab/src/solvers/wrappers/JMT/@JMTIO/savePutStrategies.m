function [simDoc, section] = savePutStrategies(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEPUTSTRATEGIES(SIMDOC, SECTION, ind)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
queuePutStrategyNode = simDoc.createElement('parameter');
queuePutStrategyNode.setAttribute('array', 'true');
queuePutStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategy');
queuePutStrategyNode.setAttribute('name', 'QueuePutStrategy');

sn = self.getStruct;
numOfClasses = sn.nclasses;
exportClasses = self.getExportableClasses();
i = sn.nodeToStation(ind);

for r=1:numOfClasses
    % Skip classes that should not be exported to JMT
    if ~exportClasses(r)
        continue;
    end

    refClassNode2 = simDoc.createElement('refClass');
    refClassNode2.appendChild(simDoc.createTextNode(sn.classnames{r}));
    queuePutStrategyNode.appendChild(refClassNode2);
    % Different to savePutStrategy.    
    switch sn.sched(i) %contrary to JMT, we assume identical across classes
        case SchedStrategy.SIRO
            subParameterNode2 = simDoc.createElement('subParameter');
            subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.RandStrategy');
            subParameterNode2.setAttribute('name', 'RandStrategy');
        case SchedStrategy.LCFS
            subParameterNode2 = simDoc.createElement('subParameter');
            subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.HeadStrategy');
            subParameterNode2.setAttribute('name', 'HeadStrategy');
        otherwise % treat as FCFS - this is required for PS
            subParameterNode2 = simDoc.createElement('subParameter');
            subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy');
            subParameterNode2.setAttribute('name', 'TailStrategy');
    end
    queuePutStrategyNode.appendChild(subParameterNode2);
    section.appendChild(queuePutStrategyNode);
end
end
