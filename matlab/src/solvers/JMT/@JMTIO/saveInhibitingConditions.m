
function [simDoc, section] = saveInhibitingConditions(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEINHIBITINGCONDITIONS(SIMDOC, SECTION, NODEIDX)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

inhibitingConditionsNode = simDoc.createElement('parameter');
inhibitingConditionsNode.setAttribute('array', 'true');
inhibitingConditionsNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix');
inhibitingConditionsNode.setAttribute('name', 'inhibitingConditions');

sn = self.getStruct;
inputs = [find(sn.connmatrix(:,ind))];
connections = {sn.nodenames{inputs}};
numOfInputs = length(connections);

numOfClasses = sn.nclasses;
numOfModes = sn.nodeparam{ind}.nmodes;
for m=1:numOfModes
    subInhibitingConditionNode = simDoc.createElement('subParameter');
    subInhibitingConditionNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix');
    subInhibitingConditionNode.setAttribute('name', 'inhibitingCondition');

    subInhibitingVectorsNode = simDoc.createElement('subParameter');
    subInhibitingVectorsNode.setAttribute('array', 'true');
    subInhibitingVectorsNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector');
    subInhibitingVectorsNode.setAttribute('name', 'inhibitingVectors');

    for k=1:numOfInputs
        subInhibitingVectorNode = simDoc.createElement('subParameter');
        subInhibitingVectorNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector');
        subInhibitingVectorNode.setAttribute('name', 'inhibitingVector');

        subStationNameNode = simDoc.createElement('subParameter');
        subStationNameNode.setAttribute('classPath', 'java.lang.String');
        subStationNameNode.setAttribute('name', 'stationName');

        placeNameValueNode = simDoc.createElement('value');
        placeNameValueNode.appendChild(simDoc.createTextNode(connections(k)));
        subStationNameNode.appendChild(placeNameValueNode);

        subInhibitingVectorNode.appendChild(subStationNameNode);

        subInhibitingEntriesNode = simDoc.createElement('subParameter');
        subInhibitingEntriesNode.setAttribute('array', 'true');
        subInhibitingEntriesNode.setAttribute('classPath', 'java.lang.Integer');
        subInhibitingEntriesNode.setAttribute('name', 'inhibitingEntries');

        for r=1:numOfClasses
            refClassNode = simDoc.createElement('refClass');
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames{r}));
            subInhibitingEntriesNode.appendChild(refClassNode);

            subParameterNode = simDoc.createElement('subParameter');
            subParameterNode.setAttribute('classPath', 'java.lang.Integer');
            subParameterNode.setAttribute('name', 'inhibitingEntry');

            valueNode2 = simDoc.createElement('value');

            if isinf(sn.nodeparam{ind}.inhibiting{m}(inputs(k),r))
                valueNode2.appendChild(simDoc.createTextNode(int2str(0)));
            else
                valueNode2.appendChild(simDoc.createTextNode(int2str(sn.nodeparam{ind}.inhibiting{m}(inputs(k),r))));
            end

            subParameterNode.appendChild(valueNode2);
            subInhibitingEntriesNode.appendChild(subParameterNode);
            subInhibitingVectorNode.appendChild(subInhibitingEntriesNode);
        end
        subInhibitingVectorsNode.appendChild(subInhibitingVectorNode);
    end

    subInhibitingConditionNode.appendChild(subInhibitingVectorsNode);
    inhibitingConditionsNode.appendChild(subInhibitingConditionNode);
end

section.appendChild(inhibitingConditionsNode);
end