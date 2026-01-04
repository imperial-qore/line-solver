function [simDoc, section] = savePreemptiveWeights(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEPREEMPTIVEWEIGHTS(SIMDOC, SECTION, NODEIDX)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
visitsNode = simDoc.createElement('parameter');
visitsNode.setAttribute('array', 'true');
visitsNode.setAttribute('classPath', 'java.lang.Double');
visitsNode.setAttribute('name', 'serviceWeights');

sn = self.getStruct;
numOfClasses = sn.nclasses;
exportClasses = self.getExportableClasses();
i = sn.nodeToStation(ind);
for r=1:numOfClasses
    % Skip classes that should not be exported to JMT
    if ~exportClasses(r)
        continue;
    end

    refClassNode = simDoc.createElement('refClass');
    refClassNode.appendChild(simDoc.createTextNode(sn.classnames{r}));
    visitsNode.appendChild(refClassNode);
    
    subParameterNode = simDoc.createElement('subParameter');
    subParameterNode.setAttribute('classPath', 'java.lang.Double');
    subParameterNode.setAttribute('name', 'serviceWeight');
    
    valueNode2 = simDoc.createElement('value');    
    %switch sn.sched(i)
        %case SchedStrategy.PS
        %    valueNode2.appendChild(simDoc.createTextNode(int2str(1)));
        %case {SchedStrategy.DPS, SchedStrategy.GPS}
            valueNode2.appendChild(simDoc.createTextNode(num2str(sn.schedparam(i,r))));
    %end
    
    subParameterNode.appendChild(valueNode2);
    visitsNode.appendChild(subParameterNode);
    section.appendChild(visitsNode);
end
end
