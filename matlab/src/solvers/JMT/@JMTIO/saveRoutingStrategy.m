function [simDoc, section] = saveRoutingStrategy(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEROUTINGSTRATEGY(SIMDOC, SECTION, NODEIDX)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
strategyNode = simDoc.createElement('parameter');
strategyNode.setAttribute('array', 'true');
strategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.RoutingStrategy');
strategyNode.setAttribute('name', 'RoutingStrategy');

sn = self.getStruct;
M = sn.nnodes;
K = sn.nclasses;
exportClasses = self.getExportableClasses();
i = ind;
% since the class switch node always outputs to a single node, it is faster
% to translate it to RAND. Also some problems with sn.rt value otherwise.
if sn.nodetype(i) == NodeType.ClassSwitch
    sn.routing(i,:) = RoutingStrategy.RAND;
end
for r=1:K
    % Skip classes that should not be exported to JMT
    if ~exportClasses(r)
        continue;
    end

    refClassNode = simDoc.createElement('refClass');
    refClassNode.appendChild(simDoc.createTextNode(sn.classnames{r}));
    strategyNode.appendChild(refClassNode);
    switch sn.routing(i,r)
        case RoutingStrategy.RAND
            concStratNode = simDoc.createElement('subParameter');
            concStratNode.setAttribute('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.RandomStrategy');
            concStratNode.setAttribute('name', 'Random');
        case RoutingStrategy.RROBIN
            concStratNode = simDoc.createElement('subParameter');
            concStratNode.setAttribute('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.RoundRobinStrategy');
            concStratNode.setAttribute('name', 'Round Robin');
        case RoutingStrategy.JSQ
            concStratNode = simDoc.createElement('subParameter');
            concStratNode.setAttribute('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.ShortestQueueLengthRoutingStrategy');
            concStratNode.setAttribute('name', 'Join the Shortest Queue (JSQ)');
        case RoutingStrategy.KCHOICES
            concStratNode = simDoc.createElement('subParameter');
            concStratNode.setAttribute('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.PowerOfKRoutingStrategy');
            concStratNode.setAttribute('name', 'Power of k');         
            concStratNode2 = simDoc.createElement('subParameter');
            concStratNode2.setAttribute('classPath', 'java.lang.Integer');
            concStratNode2.setAttribute('name', 'k');
            concStratNode2ValueNode = simDoc.createElement('value');
            concStratNode2ValueNode.appendChild(simDoc.createTextNode(sprintf('%d',sn.nodeparam{ind}{r}.k)));
            concStratNode2.appendChild(concStratNode2ValueNode);
            concStratNode3 = simDoc.createElement('subParameter');
            concStratNode3.setAttribute('classPath', 'java.lang.Boolean');
            concStratNode3.setAttribute('name', 'withMemory');
            concStratNode3ValueNode = simDoc.createElement('value');
            if sn.nodeparam{ind}{r}.withMemory
                concStratNode3ValueNode.appendChild(simDoc.createTextNode('true'));
            else
                concStratNode3ValueNode.appendChild(simDoc.createTextNode('false'));
            end
            concStratNode3.appendChild(concStratNode3ValueNode);
            concStratNode.appendChild(concStratNode2);
            concStratNode.appendChild(concStratNode3);
        case RoutingStrategy.WRROBIN
            concStratNode = simDoc.createElement('subParameter');
            concStratNode.setAttribute('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.WeightedRoundRobinStrategy');
            concStratNode.setAttribute('name', 'Weighted Round Robin');
            concStratNode2 = simDoc.createElement('subParameter');
            concStratNode2.setAttribute('array', 'true');
            concStratNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.WeightEntry');
            concStratNode2.setAttribute('name', 'WeightEntryArray');
            for j=find(sn.connmatrix(i,:)) % linked stations          
                    weight = sn.nodeparam{ind}{r}.weights(j);                    
                    concStratNode3 = simDoc.createElement('subParameter');
                    concStratNode3.setAttribute('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.WeightEntry');
                    concStratNode3.setAttribute('name', 'WeightEntry');
                    concStratNode4Station = simDoc.createElement('subParameter');
                    concStratNode4Station.setAttribute('classPath', 'java.lang.String');
                    concStratNode4Station.setAttribute('name', 'stationName');
                    concStratNode4StationValueNode = simDoc.createElement('value');
                    concStratNode4StationValueNode.appendChild(simDoc.createTextNode(sprintf('%s',sn.nodenames{j})));
                    concStratNode4Station.appendChild(concStratNode4StationValueNode);
                    concStratNode3.appendChild(concStratNode4Station);
                    concStratNode4Weight = simDoc.createElement('subParameter');
                    concStratNode4Weight.setAttribute('classPath', 'java.lang.Integer');
                    concStratNode4Weight.setAttribute('name', 'weight');
                    concStratNode4WeightValueNode = simDoc.createElement('value');
                    concStratNode4WeightValueNode.appendChild(simDoc.createTextNode(sprintf('%d',weight)));
                    concStratNode4Weight.appendChild(concStratNode4WeightValueNode);                    
                    concStratNode3.appendChild(concStratNode4Station);
                    concStratNode3.appendChild(concStratNode4Weight);
                    concStratNode2.appendChild(concStratNode3);
            end
            concStratNode.appendChild(concStratNode2);            
        case RoutingStrategy.PROB
            concStratNode = simDoc.createElement('subParameter');
            concStratNode.setAttribute('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.EmpiricalStrategy');
            concStratNode.setAttribute('name', RoutingStrategy.toText(RoutingStrategy.PROB));
            concStratNode2 = simDoc.createElement('subParameter');
            concStratNode2.setAttribute('array', 'true');
            concStratNode2.setAttribute('classPath', 'jmt.engine.random.EmpiricalEntry');
            concStratNode2.setAttribute('name', 'EmpiricalEntryArray');
            for j=find(sn.connmatrix(i,:)) % linked stations
                probRouting = sn.rtnodes((i-1)*K+r,(j-1)*K+r);
                if probRouting > 0
                    concStratNode3 = simDoc.createElement('subParameter');
                    concStratNode3.setAttribute('classPath', 'jmt.engine.random.EmpiricalEntry');
                    concStratNode3.setAttribute('name', 'EmpiricalEntry');
                    concStratNode4Station = simDoc.createElement('subParameter');
                    concStratNode4Station.setAttribute('classPath', 'java.lang.String');
                    concStratNode4Station.setAttribute('name', 'stationName');
                    concStratNode4StationValueNode = simDoc.createElement('value');
                    concStratNode4StationValueNode.appendChild(simDoc.createTextNode(sprintf('%s',sn.nodenames{j})));
                    concStratNode4Station.appendChild(concStratNode4StationValueNode);
                    concStratNode3.appendChild(concStratNode4Station);
                    concStratNode4Probability = simDoc.createElement('subParameter');
                    concStratNode4Probability.setAttribute('classPath', 'java.lang.Double');
                    concStratNode4Probability.setAttribute('name', 'probability');
                    concStratNode4ProbabilityValueNode = simDoc.createElement('value');
                    concStratNode4ProbabilityValueNode.appendChild(simDoc.createTextNode(sprintf('%12.12f',probRouting)));
                    concStratNode4Probability.appendChild(concStratNode4ProbabilityValueNode);                    
                    concStratNode3.appendChild(concStratNode4Station);
                    concStratNode3.appendChild(concStratNode4Probability);
                    concStratNode2.appendChild(concStratNode3);
                end
            end
            concStratNode.appendChild(concStratNode2);
        otherwise
            concStratNode = simDoc.createElement('subParameter');
            concStratNode.setAttribute('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.DisabledRoutingStrategy');
            concStratNode.setAttribute('name', 'Random');
    end
    strategyNode.appendChild(concStratNode);
    section.appendChild(strategyNode);
end
end
