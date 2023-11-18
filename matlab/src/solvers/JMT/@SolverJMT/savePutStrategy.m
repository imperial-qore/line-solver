function [simDoc, section] = savePutStrategy(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEPUTSTRATEGY(SIMDOC, SECTION, CURRENTNODE)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.
queuePutStrategyNode = simDoc.createElement('parameter');
queuePutStrategyNode.setAttribute('array', 'true');
queuePutStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategy');
queuePutStrategyNode.setAttribute('name', 'QueuePutStrategy');

sn = self.getStruct;
numOfClasses = sn.nclasses;
for r=1:numOfClasses
    refClassNode2 = simDoc.createElement('refClass');
    refClassNode2.appendChild(simDoc.createTextNode(sn.classnames{r}));
    
    queuePutStrategyNode.appendChild(refClassNode2);
    %                switch currentNode.input.inputJobClasses{i}{2}
    
    if ~sn.isstation(ind) % if not a station treat as FCFS
        subParameterNode2 = simDoc.createElement('subParameter');
        subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy');
        subParameterNode2.setAttribute('name', 'TailStrategy');
    else % if a station
        switch sn.schedid(sn.nodeToStation(ind))
            case SchedStrategy.ID_SIRO
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.RandStrategy');
                subParameterNode2.setAttribute('name', 'RandStrategy');
            case SchedStrategy.ID_LJF
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.LJFStrategy');
                subParameterNode2.setAttribute('name', 'LJFStrategy');
            case SchedStrategy.ID_SJF
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.SJFStrategy');
                subParameterNode2.setAttribute('name', 'SJFStrategy');
            case SchedStrategy.ID_LEPT
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.LEPTStrategy');
                subParameterNode2.setAttribute('name', 'LEPTStrategy');
            case SchedStrategy.ID_SEPT
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.SEPTStrategy');
                subParameterNode2.setAttribute('name', 'SEPTStrategy');
            case SchedStrategy.ID_LCFS
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.HeadStrategy');
                subParameterNode2.setAttribute('name', 'HeadStrategy');
            case SchedStrategy.ID_LCFSPR
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.LCFSPRStrategy');
                subParameterNode2.setAttribute('name', 'LCFSPRStrategy');
            case SchedStrategy.ID_HOL
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategyPriority');
                subParameterNode2.setAttribute('name', 'TailStrategyPriority');
            otherwise % treat as FCFS - this is required for PS
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy');
                subParameterNode2.setAttribute('name', 'TailStrategy');
        end
    end
    queuePutStrategyNode.appendChild(subParameterNode2);
    section.appendChild(queuePutStrategyNode);
end
end
