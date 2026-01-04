function [simDoc, section] = saveEnablingConditions(self, simDoc, section, ind)
% SAVEENABLINGCONDITIONS Writes only relevant enabling vectors to XML.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

enablingNode = simDoc.createElement('parameter');
enablingNode.setAttribute('array', 'true');
enablingNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix');
enablingNode.setAttribute('name', 'enablingConditions');

sn = self.getStruct;
numOfNodes = sn.nnodes;
numOfClasses = sn.nclasses;
numOfModes = sn.nodeparam{ind}.nmodes;

for m = 1:numOfModes
    subEnablingConditionNode = simDoc.createElement('subParameter');
    subEnablingConditionNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix');
    subEnablingConditionNode.setAttribute('name', 'enablingCondition');

    subEnablingVectorsNode = simDoc.createElement('subParameter');
    subEnablingVectorsNode.setAttribute('array', 'true');
    subEnablingVectorsNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector');
    subEnablingVectorsNode.setAttribute('name', 'enablingVectors');

    for k = 1:numOfNodes
        hasRelevantEntry = false;

        % Check if there is any relevant enabling entry for this node
        for r = 1:numOfClasses
            en_val = sn.nodeparam{ind}.enabling{m}(k, r);
            in_val = sn.nodeparam{ind}.inhibiting{m}(k, r);
            if (isfinite(en_val) && en_val > 0) || (isfinite(in_val) && in_val > 0)
                hasRelevantEntry = true;
                break;
            end
        end

        if hasRelevantEntry
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

            for r = 1:numOfClasses
                refClassNode = simDoc.createElement('refClass');
                refClassNode.appendChild(simDoc.createTextNode(sn.classnames{r}));
                subEnablingEntriesNode.appendChild(refClassNode);

                subParameterNode = simDoc.createElement('subParameter');
                subParameterNode.setAttribute('classPath', 'java.lang.Integer');
                subParameterNode.setAttribute('name', 'enablingEntry');

                valueNode = simDoc.createElement('value');
                en_val = sn.nodeparam{ind}.enabling{m}(k, r);
                if isinf(en_val)
                    valueNode.appendChild(simDoc.createTextNode('-1'));
                else
                    valueNode.appendChild(simDoc.createTextNode(int2str(en_val)));
                end

                subParameterNode.appendChild(valueNode);
                subEnablingEntriesNode.appendChild(subParameterNode);
            end

            subEnablingVectorNode.appendChild(subEnablingEntriesNode);
            subEnablingVectorsNode.appendChild(subEnablingVectorNode);
        end
    end

    subEnablingConditionNode.appendChild(subEnablingVectorsNode);
    enablingNode.appendChild(subEnablingConditionNode);
end

section.appendChild(enablingNode);
end
