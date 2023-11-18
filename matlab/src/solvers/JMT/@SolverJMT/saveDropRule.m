function [simDoc, section] = saveDropRule(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEDROPRULE(SIMDOC, SECTION, NODEIDX)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

% @todo unfinished

sn = self.getStruct;
schedStrategyNode = simDoc.createElement('parameter');
schedStrategyNode.setAttribute('array', 'true');
schedStrategyNode.setAttribute('classPath', 'java.lang.String');
schedStrategyNode.setAttribute('name', 'dropRules');

numOfClasses = sn.nclasses;
i = sn.nodeToStation(ind);
for r=1:numOfClasses
    refClassNode = simDoc.createElement('refClass');
    refClassNode.appendChild(simDoc.createTextNode(sn.classnames{r}));
    schedStrategyNode.appendChild(refClassNode);
    
    subParameterNode = simDoc.createElement('subParameter');
    subParameterNode.setAttribute('classPath', 'java.lang.String');
    subParameterNode.setAttribute('name', 'dropRule');
    
    valueNode2 = simDoc.createElement('value');    
    valueNode2.appendChild(simDoc.createTextNode(DropStrategy.toText(DropStrategy.fromId(sn.dropid(i,r)))));     
    subParameterNode.appendChild(valueNode2);
    schedStrategyNode.appendChild(subParameterNode);
    section.appendChild(schedStrategyNode);
end
end
