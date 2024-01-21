function [simDoc, section] = saveEnablingConditions(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEFIRINGOUTCOMES(SIMDOC, SECTION, NODEIDX)

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

enablingNode = simDoc.createElement('parameter');
enablingNode.setAttribute('array', 'true');
enablingNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix');
enablingNode.setAttribute('name', 'enablingConditions');

sn = self.getStruct;
numOfNodes = sn.nnodes;
numOfClasses = sn.nclasses;
numOfModes = sn.nodeparam{ind}.nmodes;
for m=1:(numOfModes)
    
    subEnablingConditionNode = simDoc.createElement('subParameter');
    subEnablingConditionNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix');
    subEnablingConditionNode.setAttribute('name', 'enablingCondition');
    
    subEnablingVectorsNode = simDoc.createElement('subParameter');
    subEnablingVectorsNode.setAttribute('array', 'true');
    subEnablingVectorsNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector');
    subEnablingVectorsNode.setAttribute('name', 'enablingVectors');
    
    for k=1:(numOfNodes)
        subEnablingVectorNode = simDoc.createElement('subParameter');
        subEnablingVectorNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector');
        subEnablingVectorNode.setAttribute('name', 'enablingVector');
        
        subStationNameNode = simDoc.createElement('subParameter');
        subStationNameNode.setAttribute('classPath', 'java.lang.String');
        subStationNameNode.setAttribute('name', 'stationName');
        
        placeNameValueNode = simDoc.createElement('value');
        placeNameValueNode.appendChild(simDoc.createTextNode(sn.nodenames{k}));
        subStationNameNode.appendChild(placeNameValueNode);
        
        subEnablingVectorNode.appendChild(subStationNameNode);
        
        subEnablingEntriesNode = simDoc.createElement('subParameter');
        subEnablingEntriesNode.setAttribute('array', 'true');
        subEnablingEntriesNode.setAttribute('classPath', 'java.lang.Integer');
        subEnablingEntriesNode.setAttribute('name', 'enablingEntries');
        
        exists = false;
        
        for r=1:(numOfClasses)
            refClassNode = simDoc.createElement('refClass');
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames{r}));
            subEnablingEntriesNode.appendChild(refClassNode);
            
            subParameterNode = simDoc.createElement('subParameter');
            subParameterNode.setAttribute('classPath', 'java.lang.Integer');
            subParameterNode.setAttribute('name', 'enablingEntry');
            
            valueNode2 = simDoc.createElement('value');
            
            if isinf(sn.nodeparam{ind}.enabling{m}(k,r))
                valueNode2.appendChild(simDoc.createTextNode(int2str(-1)));
                exists = true;
            elseif sn.nodeparam{ind}.enabling{m}(k,r) > 0
                valueNode2.appendChild(simDoc.createTextNode(int2str(sn.nodeparam{ind}.enabling{m}(k,r))));
                exists = true;
            elseif ~isinf(sn.nodeparam{ind}.inhibiting{m}(k,r)) && sn.nodeparam{ind}.inhibiting{m}(k,r) > 0
                valueNode2.appendChild(simDoc.createTextNode(int2str(0)));
                exists = true;
            end
            
            subParameterNode.appendChild(valueNode2);
            subEnablingEntriesNode.appendChild(subParameterNode);
            subEnablingVectorNode.appendChild(subEnablingEntriesNode);
        end
        if exists
            subEnablingVectorsNode.appendChild(subEnablingVectorNode);
        end
        
    end
    
    subEnablingConditionNode.appendChild(subEnablingVectorsNode);
    enablingNode.appendChild(subEnablingConditionNode);
end
section.appendChild(enablingNode);
end