function [simDoc, section] = saveGetStrategy(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEGETSTRATEGY(SIMDOC, SECTION, NODEIDX)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.getStruct;
ist = sn.nodeToStation(ind);
if sn.nodetype(ind) == NodeType.Queue && sn.sched(ist) == SchedStrategy.POLLING
    queueGetStrategyNode = simDoc.createElement('parameter');
    switch sn.nodeparam{ind}{1}.pollingType
        case PollingType.GATED
            queueGetStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.QueueGetStrategies.GatedPollingGetStrategy');
        case PollingType.EXHAUSTIVE
            queueGetStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.QueueGetStrategies.ExhaustivePollingGetStrategy');
        case PollingType.KLIMITED
            queueGetStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.QueueGetStrategies.LimitedPollingGetStrategy');
            pollingKNode = simDoc.createElement('subParameter');
            pollingKNode.setAttribute('classPath', 'java.lang.Integer');
            pollingKNode.setAttribute('name', 'pollingKValue');
            sn = self.getStruct;
            valueNode = simDoc.createElement('value');
            valueNode.appendChild(simDoc.createTextNode(int2str(sn.nodeparam{ind}{1}.pollingPar(1))));

            pollingKNode.appendChild(valueNode);
            queueGetStrategyNode.appendChild(pollingKNode);
    end
    queueGetStrategyNode.setAttribute('name', 'FCFSstrategy');
    section.appendChild(queueGetStrategyNode);
else
    % the get strategy is always fcfs for queues unless a polling queue
    queueGetStrategyNode = simDoc.createElement('parameter');
    queueGetStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy');
    queueGetStrategyNode.setAttribute('name', 'FCFSstrategy');
    section.appendChild(queueGetStrategyNode);
end
end