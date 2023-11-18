function [simDoc, section] = saveDropStrategy(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEDROPSTRATEGY(SIMDOC, SECTION, NODEIDX)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

% @todo unfinished

% This is just creating an empty drop node
sn = self.getStruct;
numOfClasses = sn.nclasses;

schedStrategyNode = simDoc.createElement('parameter');
schedStrategyNode.setAttribute('array', 'true');
schedStrategyNode.setAttribute('classPath', 'java.lang.String');
schedStrategyNode.setAttribute('name', 'dropStrategies');
i = sn.nodeToStation(ind);
for r=1:numOfClasses
    refClassNode = simDoc.createElement('refClass');
    refClassNode.appendChild(simDoc.createTextNode(sn.classnames{r}));
    schedStrategyNode.appendChild(refClassNode);
    
    subParameterNode = simDoc.createElement('subParameter');
    subParameterNode.setAttribute('classPath', 'java.lang.String');
    subParameterNode.setAttribute('name', 'dropStrategy');
    
    valueNode2 = simDoc.createElement('value');
    if isnan(i) || sn.dropid(i,r)==0
        % JMT sets the field to 'drop' for nodes without a buffer
        valueNode2.appendChild(simDoc.createTextNode('drop'));
    else
        valueNode2.appendChild(simDoc.createTextNode(DropStrategy.toText(DropStrategy.fromId(sn.dropid(i,r)))));
    end
    
    subParameterNode.appendChild(valueNode2);
    schedStrategyNode.appendChild(subParameterNode);
    section.appendChild(schedStrategyNode);
end
end
