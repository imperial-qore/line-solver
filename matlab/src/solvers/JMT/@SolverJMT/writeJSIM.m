function outputFileName = writeJSIM(self, sn, outputFileName)
% FNAME = WRITEJSIM(SN, FNAME)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

if nargin<3 %~exist('outFileName','var')
    outputFileName = getJSIMTempPath(self);
end

[simXMLElem, simXMLDoc] = saveXMLHeader(self, self.model.getLogPath);
[simXMLElem, simXMLDoc] = saveClasses(self, simXMLElem, simXMLDoc);

numOfClasses = sn.nclasses;
numOfNodes = sn.nnodes;
for i=1:numOfNodes
    ind = i;
    currentNode = self.model.nodes{i,1};
    node = simXMLDoc.createElement('node');
    node.setAttribute('name', currentNode.name);
    
    nodeSections = getSections(currentNode);
    for j=1:length(nodeSections)
        xml_section = simXMLDoc.createElement('section');
        currentSection = nodeSections{1,j};
        if ~isempty(currentSection)            
            xml_section.setAttribute('className', currentSection.className);
            switch currentSection.className
                case 'Buffer'
                    xml_section.setAttribute('className', 'Queue'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveBufferCapacity(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveDropStrategy(self, simXMLDoc, xml_section, ind); % unfinished
                    [simXMLDoc, xml_section] = saveGetStrategy(self, simXMLDoc, xml_section);
                    [simXMLDoc, xml_section] = savePutStrategy(self, simXMLDoc, xml_section, ind);
                case {'Server','PreemptiveServer'}
                    [simXMLDoc, xml_section] = saveNumberOfServers(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveServerVisits(self, simXMLDoc, xml_section);
                    [simXMLDoc, xml_section] = saveServiceStrategy(self, simXMLDoc, xml_section, ind);
                case 'SharedServer'
                    xml_section.setAttribute('className', 'PSServer'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveNumberOfServers(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveServerVisits(self, simXMLDoc, xml_section);
                    [simXMLDoc, xml_section] = saveServiceStrategy(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = savePreemptiveStrategy(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = savePreemptiveWeights(self, simXMLDoc, xml_section, ind);
                case 'InfiniteServer'
                    xml_section.setAttribute('className', 'Delay'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveServiceStrategy(self, simXMLDoc, xml_section, ind);
                case 'RandomSource'
                    [simXMLDoc, xml_section] = saveArrivalStrategy(self, simXMLDoc, xml_section, ind);
                case 'Dispatcher'
                    xml_section.setAttribute('className', 'Router'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveRoutingStrategy(self, simXMLDoc, xml_section, ind);
                case 'StatelessClassSwitcher'
                    xml_section.setAttribute('className', 'ClassSwitch'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveClassSwitchStrategy(self, simXMLDoc, xml_section, ind);
                case 'LogTunnel'
                    [simXMLDoc, xml_section] = saveLogTunnel(self, simXMLDoc, xml_section, ind);
                case 'Joiner'
                    xml_section.setAttribute('className', 'Join'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveJoinStrategy(self, simXMLDoc, xml_section, ind);
                case 'Forker'
                    xml_section.setAttribute('className', 'Fork'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveForkStrategy(self, simXMLDoc, xml_section, ind);
                case 'Storage'
                    xml_section.setAttribute('className', 'Storage'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveTotalCapacity(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = savePlaceCapacities(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveDropRule(self, simXMLDoc, xml_section, ind); % unfinished
                    [simXMLDoc, xml_section] = saveGetStrategy(self, simXMLDoc, xml_section);
                    [simXMLDoc, xml_section] = savePutStrategies(self, simXMLDoc, xml_section, ind);
                case 'Enabling'
                    xml_section.setAttribute('className', 'Enabling'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveEnablingConditions(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveInhibitingConditions(self, simXMLDoc, xml_section, ind);
                case 'Firing'
                    xml_section.setAttribute('className', 'Firing'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveFiringOutcomes(self, simXMLDoc, xml_section, ind);
                case 'Timing'
                    xml_section.setAttribute('className', 'Timing'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveModeNames(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveNumbersOfServers(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveTimingStrategies(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveFiringPriorities(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveFiringWeights(self, simXMLDoc, xml_section, ind);
            end
            node.appendChild(xml_section);
        end
    end
    simXMLElem.appendChild(node);
end

[simXMLElem, simXMLDoc] = saveMetrics(self, simXMLElem, simXMLDoc);
[simXMLElem, simXMLDoc] = saveLinks(self, simXMLElem, simXMLDoc);
[simXMLElem, simXMLDoc] = saveRegions(self, simXMLElem, simXMLDoc);

hasReferenceNodes = false;
preloadNode = simXMLDoc.createElement('preload');
s0 = sn.state;
numOfStations = sn.nstations;
for i=1:numOfStations
    isReferenceNode = false;
    if sn.nodetype(sn.stationToNode(i))~=NodeType.Source && sn.nodetype(sn.stationToNode(i))~=NodeType.Join
        [~, nir] = State.toMarginal(sn, sn.stationToNode(i), s0{sn.stationToStateful(i)});
        stationPopulationsNode = simXMLDoc.createElement('stationPopulations');
        stationPopulationsNode.setAttribute('stationName', sn.nodenames{sn.stationToNode(i)});
        for r=1:numOfClasses
            classPopulationNode = simXMLDoc.createElement('classPopulation');
            if isinf(sn.njobs(r))
                % case 'open'
                isReferenceNode = true;
                classPopulationNode.setAttribute('population', sprintf('%d',round(nir(r))));
                classPopulationNode.setAttribute('refClass', sn.classnames{r});
                stationPopulationsNode.appendChild(classPopulationNode);
            else
                % case 'closed'
                isReferenceNode = true;
                classPopulationNode.setAttribute('population', sprintf('%d',round(nir(r))));
                classPopulationNode.setAttribute('refClass', sn.classnames{r});
                stationPopulationsNode.appendChild(classPopulationNode);
            end
        end
    end
    if isReferenceNode
        preloadNode.appendChild(stationPopulationsNode);
    end
    hasReferenceNodes = hasReferenceNodes + isReferenceNode;
end
if hasReferenceNodes
    simXMLElem.appendChild(preloadNode);
end

try
    xmlwrite(outputFileName, simXMLDoc);
catch ME
    ME
    ME.stack
    javaaddpath(which('xercesImpl-2.11.0.jar'));
    javaaddpath(which('xml-apis-2.11.0.jar'));
    pkg load io;
    xmlwrite(outputFileName, simXMLDoc);
end
end