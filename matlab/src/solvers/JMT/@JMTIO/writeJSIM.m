function outputFileName = writeJSIM(self, sn, outputFileName)
% FNAME = WRITEJSIM(SN, FNAME)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin<3 || isempty(outputFileName)
    % Generate default temp path
    filePath = lineTempName('jsim');
    outputFileName = [filePath, filesep, 'model.jsim'];
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
            sectionClassName = currentSection.className;
            % Override className for preemptive strategies - JMT requires PreemptiveServer
            if strcmp(sectionClassName, 'Server') && isa(currentNode, 'Queue')
                sched = currentNode.schedStrategy;
                if sched == SchedStrategy.SRPT || sched == SchedStrategy.SRPTPRIO || ...
                   sched == SchedStrategy.LCFSPR || sched == SchedStrategy.LCFSPRPRIO || ...
                   sched == SchedStrategy.LCFSPI || sched == SchedStrategy.LCFSPIPRIO || ...
                   sched == SchedStrategy.FCFSPR || sched == SchedStrategy.FCFSPRPRIO || ...
                   sched == SchedStrategy.FCFSPI || sched == SchedStrategy.FCFSPIPRIO
                    sectionClassName = 'PreemptiveServer';
                end
            end
            xml_section.setAttribute('className', sectionClassName);
            switch sectionClassName
                case 'Buffer'
                    xml_section.setAttribute('className', 'Queue'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveBufferCapacity(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveDropStrategy(self, simXMLDoc, xml_section, ind); % unfinished
                    [simXMLDoc, xml_section] = saveGetStrategy(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = savePutStrategy(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveImpatience(self, simXMLDoc, xml_section, ind);
                case 'Server'
                    [simXMLDoc, xml_section] = saveNumberOfServers(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveServerVisits(self, simXMLDoc, xml_section);
                    % Heterogeneous server configuration
                    [simXMLDoc, xml_section] = saveServerTypeNames(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveServersPerType(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveServerCompatibilities(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveHeteroSchedPolicy(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveServiceStrategy(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveDelayOffStrategy(self, simXMLDoc, xml_section, ind);
                    % Note: SwitchoverStrategy is only supported by PollingServer in JMT
                    % Regular Server class does not support switchover times
                    if ~isempty(currentNode.switchoverTime)
                        line_warning(mfilename, 'JMT does not support switchover times for non-polling queues. Switchover times will be ignored for node %s.', currentNode.name);
                    end
                case 'PreemptiveServer'
                    [simXMLDoc, xml_section] = saveNumberOfServers(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveServerVisits(self, simXMLDoc, xml_section);
                    [simXMLDoc, xml_section] = saveServiceStrategy(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveDelayOffStrategy(self, simXMLDoc, xml_section, ind);
                case 'SharedServer'
                    xml_section.setAttribute('className', 'PSServer'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveNumberOfServers(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveServerVisits(self, simXMLDoc, xml_section);
                    [simXMLDoc, xml_section] = saveServiceStrategy(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveDelayOffStrategy(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = savePreemptiveStrategy(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = savePreemptiveWeights(self, simXMLDoc, xml_section, ind);
                case 'InfiniteServer'
                    xml_section.setAttribute('className', 'Delay'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveServiceStrategy(self, simXMLDoc, xml_section, ind);
                case 'PollingServer'
                    % we assume identical polling type on each buffer
                    switch sn.nodeparam{ind}{1}.pollingType
                        case PollingType.GATED
                            xml_section.setAttribute('className', 'GatedPollingServer'); %overwrite with JMT class name
                        case PollingType.EXHAUSTIVE
                            xml_section.setAttribute('className', 'ExhaustivePollingServer'); %overwrite with JMT class name
                        case PollingType.KLIMITED
                            xml_section.setAttribute('className', 'LimitedPollingServer'); %overwrite with JMT class name
                    end
                    [simXMLDoc, xml_section] = saveNumberOfServers(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveServerVisits(self, simXMLDoc, xml_section);
                    [simXMLDoc, xml_section] = saveServiceStrategy(self, simXMLDoc, xml_section, ind);
                    [simXMLDoc, xml_section] = saveSwitchoverStrategy(self, simXMLDoc, xml_section, ind);
                case 'RandomSource'
                    [simXMLDoc, xml_section] = saveArrivalStrategy(self, simXMLDoc, xml_section, ind);
                case 'Dispatcher'
                    xml_section.setAttribute('className', 'Router'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveRoutingStrategy(self, simXMLDoc, xml_section, ind);
                case 'StatelessClassSwitcher'
                    xml_section.setAttribute('className', 'ClassSwitch'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveClassSwitchStrategy(self, simXMLDoc, xml_section, ind);
                case 'Cache'
                    % CacheClassSwitcher section has className='Cache'
                    [simXMLDoc, xml_section] = saveCacheStrategy(self, simXMLDoc, xml_section, ind);
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
                    [simXMLDoc, xml_section] = saveGetStrategy(self, simXMLDoc, xml_section, ind);
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
        % TODO: this assumes that the preloading comes from the first
        % state, but it could be that stateprior is put on another
        for r=1:numOfClasses
            % Skip closed classes with 0 customers, unless state has jobs via class switching
            if isfinite(sn.njobs(r)) && sn.njobs(r) == 0 && nir(1,r) == 0
                continue;
            end

            classPopulationNode = simXMLDoc.createElement('classPopulation');
            if isinf(sn.njobs(r))
                % case 'open'
                isReferenceNode = true;
                classPopulationNode.setAttribute('population', sprintf('%d',round(nir(1,r))));
                classPopulationNode.setAttribute('refClass', sn.classnames{r});
                stationPopulationsNode.appendChild(classPopulationNode);
            else
                % case 'closed'
                isReferenceNode = true;
                classPopulationNode.setAttribute('population', sprintf('%d',round(nir(1,r))));
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
    getReport(ME,'basic')
    javaaddpath(which('xercesImpl-2.11.0.jar'));
    javaaddpath(which('xml-apis-2.11.0.jar'));
    pkg load io;
    xmlwrite(outputFileName, simXMLDoc);
end
end