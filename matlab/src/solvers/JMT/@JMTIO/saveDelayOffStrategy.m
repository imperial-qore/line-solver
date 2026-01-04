function [simDoc, section] = saveDelayOffStrategy(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEDELAYOFFSTRATEGY(SIMDOC, SECTION, NODEIDX)
%
% Saves delay-off and setup time strategies for Queue nodes.
% This exports the delayOffTime and setUpTime parameters to JMT XML format
% when the queue has delay-off times enabled.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.getStruct;
K = sn.nclasses;
currentNode = self.model.getNodes{ind};

% Check if this node has delay-off enabled
if ~isa(currentNode, 'Queue') || isempty(currentNode.setupTime) || isempty(currentNode.delayoffTime)
    return;
end

% Check if any class has delay-off configured
hasDelayOff = false;
for r = 1:K
    if r <= length(currentNode.setupTime) && ~isempty(currentNode.setupTime{r})
        hasDelayOff = true;
        break;
    end
end

if ~hasDelayOff
    return;
end

% Export delayOffTime parameter
delayOffParamNode = simDoc.createElement('parameter');
delayOffParamNode.setAttribute('array', 'true');
delayOffParamNode.setAttribute('classPath', 'java.lang.Object');
delayOffParamNode.setAttribute('name', 'delayOffTime');

for r = 1:K
    refClassNode = simDoc.createElement('refClass');
    refClassNode.appendChild(simDoc.createTextNode(sn.classnames{r}));
    delayOffParamNode.appendChild(refClassNode);

    subParamNode = simDoc.createElement('subParameter');
    subParamNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy');
    subParamNode.setAttribute('name', 'delayOffTime');

    if r <= length(currentNode.delayoffTime) && ~isempty(currentNode.delayoffTime{r})
        dist = currentNode.delayoffTime{r};
        subParamNode = appendDistributionXml(simDoc, subParamNode, dist);
    else
        % Zero time (Immediate)
        subParamNode = appendZeroTimeXml(simDoc, subParamNode);
    end

    delayOffParamNode.appendChild(subParamNode);
end
section.appendChild(delayOffParamNode);

% Export setUpTime parameter
setupParamNode = simDoc.createElement('parameter');
setupParamNode.setAttribute('array', 'true');
setupParamNode.setAttribute('classPath', 'java.lang.Object');
setupParamNode.setAttribute('name', 'setUpTime');

for r = 1:K
    refClassNode = simDoc.createElement('refClass');
    refClassNode.appendChild(simDoc.createTextNode(sn.classnames{r}));
    setupParamNode.appendChild(refClassNode);

    subParamNode = simDoc.createElement('subParameter');
    subParamNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy');
    subParamNode.setAttribute('name', 'setUpTime');

    if r <= length(currentNode.setupTime) && ~isempty(currentNode.setupTime{r})
        dist = currentNode.setupTime{r};
        subParamNode = appendDistributionXml(simDoc, subParamNode, dist);
    else
        % Zero time (Immediate)
        subParamNode = appendZeroTimeXml(simDoc, subParamNode);
    end

    setupParamNode.appendChild(subParamNode);
end
section.appendChild(setupParamNode);

end

function parentNode = appendZeroTimeXml(simDoc, parentNode)
% Append zero/immediate service time XML

distributionNode = simDoc.createElement('subParameter');
distributionNode.setAttribute('classPath', 'jmt.engine.random.DeterministicDistr');
distributionNode.setAttribute('name', 'Deterministic');
parentNode.appendChild(distributionNode);

distrParNode = simDoc.createElement('subParameter');
distrParNode.setAttribute('classPath', 'jmt.engine.random.DeterministicDistrPar');
distrParNode.setAttribute('name', 'distrPar');

tNode = simDoc.createElement('subParameter');
tNode.setAttribute('classPath', 'java.lang.Double');
tNode.setAttribute('name', 't');
tValue = simDoc.createElement('value');
tValue.appendChild(simDoc.createTextNode('0.0'));
tNode.appendChild(tValue);
distrParNode.appendChild(tNode);
parentNode.appendChild(distrParNode);

end

function parentNode = appendDistributionXml(simDoc, parentNode, dist)
% Append distribution XML elements to a parent node for JMT export
% Supports Exp, Erlang, Det, and Immediate distributions

if isa(dist, 'Immediate')
    % Zero/Immediate service time
    parentNode = appendZeroTimeXml(simDoc, parentNode);

elseif isa(dist, 'Exp')
    distributionNode = simDoc.createElement('subParameter');
    distributionNode.setAttribute('classPath', 'jmt.engine.random.Exponential');
    distributionNode.setAttribute('name', 'Exponential');
    parentNode.appendChild(distributionNode);

    distrParNode = simDoc.createElement('subParameter');
    distrParNode.setAttribute('classPath', 'jmt.engine.random.ExponentialPar');
    distrParNode.setAttribute('name', 'distrPar');

    lambdaNode = simDoc.createElement('subParameter');
    lambdaNode.setAttribute('classPath', 'java.lang.Double');
    lambdaNode.setAttribute('name', 'lambda');
    lambdaValue = simDoc.createElement('value');
    lambdaValue.appendChild(simDoc.createTextNode(sprintf('%.12f', dist.getRate())));
    lambdaNode.appendChild(lambdaValue);
    distrParNode.appendChild(lambdaNode);
    parentNode.appendChild(distrParNode);

elseif isa(dist, 'Erlang')
    distributionNode = simDoc.createElement('subParameter');
    distributionNode.setAttribute('classPath', 'jmt.engine.random.Erlang');
    distributionNode.setAttribute('name', 'Erlang');
    parentNode.appendChild(distributionNode);

    distrParNode = simDoc.createElement('subParameter');
    distrParNode.setAttribute('classPath', 'jmt.engine.random.ErlangPar');
    distrParNode.setAttribute('name', 'distrPar');

    alphaNode = simDoc.createElement('subParameter');
    alphaNode.setAttribute('classPath', 'java.lang.Double');
    alphaNode.setAttribute('name', 'alpha');
    alphaValue = simDoc.createElement('value');
    meanVal = dist.getMean();
    phases = dist.getNumberOfPhases();
    alpha = phases / meanVal;
    alphaValue.appendChild(simDoc.createTextNode(sprintf('%.12f', alpha)));
    alphaNode.appendChild(alphaValue);
    distrParNode.appendChild(alphaNode);

    rNode = simDoc.createElement('subParameter');
    rNode.setAttribute('classPath', 'java.lang.Long');
    rNode.setAttribute('name', 'r');
    rValue = simDoc.createElement('value');
    rValue.appendChild(simDoc.createTextNode(sprintf('%d', phases)));
    rNode.appendChild(rValue);
    distrParNode.appendChild(rNode);
    parentNode.appendChild(distrParNode);

elseif isa(dist, 'Det')
    distributionNode = simDoc.createElement('subParameter');
    distributionNode.setAttribute('classPath', 'jmt.engine.random.DeterministicDistr');
    distributionNode.setAttribute('name', 'Deterministic');
    parentNode.appendChild(distributionNode);

    distrParNode = simDoc.createElement('subParameter');
    distrParNode.setAttribute('classPath', 'jmt.engine.random.DeterministicDistrPar');
    distrParNode.setAttribute('name', 'distrPar');

    tNode = simDoc.createElement('subParameter');
    tNode.setAttribute('classPath', 'java.lang.Double');
    tNode.setAttribute('name', 't');
    tValue = simDoc.createElement('value');
    tValue.appendChild(simDoc.createTextNode(sprintf('%.12f', dist.getMean())));
    tNode.appendChild(tValue);
    distrParNode.appendChild(tNode);
    parentNode.appendChild(distrParNode);

else
    % Default fallback: use mean as deterministic value
    distributionNode = simDoc.createElement('subParameter');
    distributionNode.setAttribute('classPath', 'jmt.engine.random.DeterministicDistr');
    distributionNode.setAttribute('name', 'Deterministic');
    parentNode.appendChild(distributionNode);

    distrParNode = simDoc.createElement('subParameter');
    distrParNode.setAttribute('classPath', 'jmt.engine.random.DeterministicDistrPar');
    distrParNode.setAttribute('name', 'distrPar');

    tNode = simDoc.createElement('subParameter');
    tNode.setAttribute('classPath', 'java.lang.Double');
    tNode.setAttribute('name', 't');
    tValue = simDoc.createElement('value');
    tValue.appendChild(simDoc.createTextNode(sprintf('%.12f', dist.getMean())));
    tNode.appendChild(tValue);
    distrParNode.appendChild(tNode);
    parentNode.appendChild(distrParNode);
end

end
