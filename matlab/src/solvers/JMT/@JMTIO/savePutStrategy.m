function [simDoc, section] = savePutStrategy(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEPUTSTRATEGY(SIMDOC, SECTION, CURRENTNODE)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
queuePutStrategyNode = simDoc.createElement('parameter');
queuePutStrategyNode.setAttribute('array', 'true');
queuePutStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategy');
queuePutStrategyNode.setAttribute('name', 'QueuePutStrategy');

sn = self.getStruct;
numOfClasses = sn.nclasses;
exportClasses = self.getExportableClasses();
for r=1:numOfClasses
    % Skip classes that should not be exported to JMT
    if ~exportClasses(r)
        continue;
    end

    refClassNode2 = simDoc.createElement('refClass');
    refClassNode2.appendChild(simDoc.createTextNode(sn.classnames{r}));
    
    queuePutStrategyNode.appendChild(refClassNode2);    
    
    if ~sn.isstation(ind) % if not a station treat as FCFS
        subParameterNode2 = simDoc.createElement('subParameter');
        subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy');
        subParameterNode2.setAttribute('name', 'TailStrategy');
    else % if a station
        switch sn.sched(sn.nodeToStation(ind))
            case SchedStrategy.SIRO
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.RandStrategy');
                subParameterNode2.setAttribute('name', 'RandStrategy');
            case SchedStrategy.LJF
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.LJFStrategy');
                subParameterNode2.setAttribute('name', 'LJFStrategy');
            case SchedStrategy.SJF
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.SJFStrategy');
                subParameterNode2.setAttribute('name', 'SJFStrategy');
            case SchedStrategy.LEPT
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.LEPTStrategy');
                subParameterNode2.setAttribute('name', 'LEPTStrategy');
            case SchedStrategy.SEPT
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.SEPTStrategy');
                subParameterNode2.setAttribute('name', 'SEPTStrategy');
            case SchedStrategy.LCFS
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.HeadStrategy');
                subParameterNode2.setAttribute('name', 'HeadStrategy');
            case SchedStrategy.LCFSPRIO
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.HeadStrategyPriority');
                subParameterNode2.setAttribute('name', 'HeadStrategyPriority');
            case SchedStrategy.LCFSPR
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.LCFSPRStrategy');
                subParameterNode2.setAttribute('name', 'LCFSPRStrategy');
            case SchedStrategy.LCFSPI
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.LCFSPIStrategy');
                subParameterNode2.setAttribute('name', 'LCFSPIStrategy');
            case SchedStrategy.LCFSPRPRIO
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.LCFSPRStrategyPriority');
                subParameterNode2.setAttribute('name', 'LCFSPRStrategyPriority');
            case SchedStrategy.LCFSPIPRIO
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.LCFSPIStrategyPriority');
                subParameterNode2.setAttribute('name', 'LCFSPIStrategyPriority');
            case SchedStrategy.FCFSPR
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.FCFSPRStrategy');
                subParameterNode2.setAttribute('name', 'FCFSPRStrategy');
            case SchedStrategy.FCFSPI
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.FCFSPIStrategy');
                subParameterNode2.setAttribute('name', 'FCFSPIStrategy');
            case SchedStrategy.FCFSPRPRIO
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.FCFSPRStrategyPriority');
                subParameterNode2.setAttribute('name', 'FCFSPRStrategyPriority');
            case SchedStrategy.FCFSPIPRIO
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.FCFSPIStrategyPriority');
                subParameterNode2.setAttribute('name', 'FCFSPIStrategyPriority');
            case {SchedStrategy.HOL, SchedStrategy.FCFSPRIO}
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategyPriority');
                subParameterNode2.setAttribute('name', 'TailStrategyPriority');
            case SchedStrategy.EDD
                % Note: JMT does not natively support EDD yet. This generates XML for future compatibility.
                % For now, this will likely fall back to FCFS behavior in JMT execution.
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.EDDStrategy');
                subParameterNode2.setAttribute('name', 'EDDStrategy');
            case SchedStrategy.EDF
                % Note: JMT does not natively support EDF yet. This generates XML for future compatibility.
                % For now, this will likely fall back to FCFS behavior in JMT execution.
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.EDFStrategy');
                subParameterNode2.setAttribute('name', 'EDFStrategy');
            case SchedStrategy.SRPT
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.SRPTStrategy');
                subParameterNode2.setAttribute('name', 'SRPTStrategy');
            case SchedStrategy.SRPTPRIO
                subParameterNode2 = simDoc.createElement('subParameter');
                subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.SRPTStrategyPriority');
                subParameterNode2.setAttribute('name', 'SRPTStrategyPriority');
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
