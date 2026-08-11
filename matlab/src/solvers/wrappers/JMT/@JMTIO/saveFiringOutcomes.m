function [simDoc, section] = saveFiringOutcomes(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEFIRINGOUTCOMES(SIMDOC, SECTION, NODEIDX)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

firingOutcomesNode = simDoc.createElement('parameter');
firingOutcomesNode.setAttribute('array', 'true');
firingOutcomesNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix');
firingOutcomesNode.setAttribute('name', 'firingOutcomes');

sn = self.getStruct;
outputs = [find(sn.connmatrix(ind,:))];
connections = {sn.nodenames{outputs}};    
numOfOutput = length(connections);

numOfModes = sn.nodeparam{ind}.nmodes;
numOfClasses = length(self.model.classes);
for m=1:numOfModes
    
    subFiringOutcomeNode = simDoc.createElement('subParameter');
    subFiringOutcomeNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix');
    subFiringOutcomeNode.setAttribute('name', 'firingOutcome');
    
    subFiringVectorsNode = simDoc.createElement('subParameter');
    subFiringVectorsNode.setAttribute('array', 'true');
    subFiringVectorsNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector');
    subFiringVectorsNode.setAttribute('name', 'firingVectors');
    
    for k=1:(numOfOutput)
        subFiringVectorNode = simDoc.createElement('subParameter');
        subFiringVectorNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector');
        subFiringVectorNode.setAttribute('name', 'firingVector');
        
        subStationNameNode = simDoc.createElement('subParameter');
        subStationNameNode.setAttribute('classPath', 'java.lang.String');
        subStationNameNode.setAttribute('name', 'stationName');
        
        placeNameValueNode = simDoc.createElement('value');
        placeNameValueNode.appendChild(simDoc.createTextNode(connections(k)));
        subStationNameNode.appendChild(placeNameValueNode);
        
        subFiringVectorNode.appendChild(subStationNameNode);
        
        subFiringEntriesNode = simDoc.createElement('subParameter');
        subFiringEntriesNode.setAttribute('array', 'true');
        subFiringEntriesNode.setAttribute('classPath', 'java.lang.Integer');
        subFiringEntriesNode.setAttribute('name', 'firingEntries');
        
        for j=1:(numOfClasses)
            currentClass = self.model.classes{j,1};
            
            refClassNode = simDoc.createElement('refClass');
            refClassNode.appendChild(simDoc.createTextNode(currentClass.name));
            subFiringEntriesNode.appendChild(refClassNode);
            
            subParameterNode = simDoc.createElement('subParameter');
            subParameterNode.setAttribute('classPath', 'java.lang.Integer');
            subParameterNode.setAttribute('name', 'firingEntry');
            
            valueNode2 = simDoc.createElement('value');
            valueNode2.appendChild(simDoc.createTextNode(int2str(sn.nodeparam{ind}.firing{m}(outputs(k),j))));
            
            subParameterNode.appendChild(valueNode2);
            subFiringEntriesNode.appendChild(subParameterNode);
            subFiringVectorNode.appendChild(subFiringEntriesNode);
        end
        subFiringVectorsNode.appendChild(subFiringVectorNode);
    end
    
    subFiringOutcomeNode.appendChild(subFiringVectorsNode);
    firingOutcomesNode.appendChild(subFiringOutcomeNode);
end
section.appendChild(firingOutcomesNode);
end
