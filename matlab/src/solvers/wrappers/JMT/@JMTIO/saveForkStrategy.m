function [simDoc, section] = saveForkStrategy(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEFORKSTRATEGY(SIMDOC, SECTION, NODEIDX)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn=self.getStruct;

% Get exportable classes (handles cache classes and class-switching)
exportClasses = self.getExportableClasses();

jplNode = simDoc.createElement('parameter');
jplNode.setAttribute('classPath', 'java.lang.Integer');
jplNode.setAttribute('name', 'jobsPerLink');
valueNode = simDoc.createElement('value');
valueNode.appendChild(simDoc.createTextNode(int2str(sn.nodeparam{ind}.fanOut)));
jplNode.appendChild(valueNode);
section.appendChild(jplNode);

blockNode = simDoc.createElement('parameter');
blockNode.setAttribute('classPath', 'java.lang.Integer');
blockNode.setAttribute('name', 'block');
valueNode = simDoc.createElement('value');
valueNode.appendChild(simDoc.createTextNode(int2str(-1)));
blockNode.appendChild(valueNode);
section.appendChild(blockNode);

% isSimplifiedFork lets JMT ignore the branch list and send one job down every
% link. That is only the same model when every branch is certain and carries the
% same number of tasks, so a variable forking level switches it off and makes
% JMT read the per-branch entries emitted below.
fanOutLink = [];
fanOutProb = [];
fanOutDist = {};
if isfield(sn.nodeparam{ind},'fanOutLink')
    fanOutLink = sn.nodeparam{ind}.fanOutLink;
    fanOutProb = sn.nodeparam{ind}.fanOutProb;
    fanOutDist = sn.nodeparam{ind}.fanOutDist;
end
isSimplified = true;
if ~isempty(fanOutLink)
    taken = fanOutProb > 0;
    isSimplified = all(fanOutProb(taken) == 1) && ...
        all(fanOutLink(taken) == sn.nodeparam{ind}.fanOut) && ...
        all(cellfun(@isempty, fanOutDist(:)));
end

issimplNode = simDoc.createElement('parameter');
issimplNode.setAttribute('classPath', 'java.lang.Boolean');
issimplNode.setAttribute('name', 'isSimplifiedFork');
valueNode = simDoc.createElement('value');
if isSimplified
    valueNode.appendChild(simDoc.createTextNode('true'));
else
    valueNode.appendChild(simDoc.createTextNode('false'));
end
issimplNode.appendChild(valueNode);
section.appendChild(issimplNode);

strategyNode = simDoc.createElement('parameter');
strategyNode.setAttribute('array', 'true');
strategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ForkStrategy');
strategyNode.setAttribute('name', 'ForkStrategy');

i = ind;
numOfClasses = sn.nclasses;
for r=1:numOfClasses
    % Skip classes that should not be exported to JMT
    if ~exportClasses(r)
        continue;
    end

    refClassNode = simDoc.createElement('refClass');
    refClassNode.appendChild(simDoc.createTextNode(sn.classnames{r}));
    strategyNode.appendChild(refClassNode);
    
    classStratNode = simDoc.createElement('subParameter');
    classStratNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ForkStrategies.ProbabilitiesFork');
    classStratNode.setAttribute('name', 'Branch Probabilities');
    classStratNode2 = simDoc.createElement('subParameter');
    classStratNode2.setAttribute('array', 'true');
    classStratNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.ForkStrategies.OutPath');
    classStratNode2.setAttribute('name', 'EmpiricalEntryArray');
    switch sn.routing(i,r)
        case RoutingStrategy.PROB
            % One OutPathEntry per outgoing link. An earlier version built the
            % entry inside this loop but appended it outside, so only the last
            % connected node survived; that was invisible because
            % isSimplifiedFork makes JMT send one job down every link and ignore
            % the branch list, and it stops being invisible as soon as a branch
            % carries its own probability or its own jobs-per-link.
            for k=find(sn.connmatrix(i,:))
                classStratNode3 = simDoc.createElement('subParameter');
                classStratNode3.setAttribute('classPath', 'jmt.engine.NetStrategies.ForkStrategies.OutPath');
                classStratNode3.setAttribute('name', 'OutPathEntry');

                classStratNode4 = simDoc.createElement('subParameter');
                classStratNode4.setAttribute('classPath', 'jmt.engine.random.EmpiricalEntry');
                classStratNode4.setAttribute('name', 'outUnitProbability');
                classStratNode4Station = simDoc.createElement('subParameter');
                classStratNode4Station.setAttribute('classPath', 'java.lang.String');
                classStratNode4Station.setAttribute('name', 'stationName');
                classStratNode4StationValueNode = simDoc.createElement('value');
                classStratNode4StationValueNode.appendChild(simDoc.createTextNode(sprintf('%s',sn.nodenames{k})));
                classStratNode4Station.appendChild(classStratNode4StationValueNode);
                classStratNode4.appendChild(classStratNode4Station);
                % branch activation probability: JMT's outUnitProbability
                if isempty(fanOutProb)
                    branchp = 1.0;
                else
                    branchp = fanOutProb(k,r);
                end
                classStratNode4Probability = simDoc.createElement('subParameter');
                classStratNode4Probability.setAttribute('classPath', 'java.lang.Double');
                classStratNode4Probability.setAttribute('name', 'probability');
                classStratNode4ProbabilityValueNode = simDoc.createElement('value');
                classStratNode4ProbabilityValueNode.appendChild(simDoc.createTextNode(num2str(branchp,'%.15g')));
                classStratNode4Probability.appendChild(classStratNode4ProbabilityValueNode);
                classStratNode4.appendChild(classStratNode4Probability);
                classStratNode3.appendChild(classStratNode4);

                % JobsPerLinkDis is an EmpiricalEntry ARRAY: one entry per point
                % of the jobs-per-link distribution. A deterministic fork emits
                % the single degenerate entry it always did.
                if ~isempty(fanOutDist) && ~isempty(fanOutDist{k,r})
                    jplPoints = fanOutDist{k,r}.getParam(2).paramValue;
                    jplProbs = fanOutDist{k,r}.getParam(1).paramValue;
                    jplProbs = jplProbs / sum(jplProbs);
                elseif ~isempty(fanOutLink)
                    jplPoints = fanOutLink(k,r);
                    jplProbs = 1.0;
                else
                    jplPoints = sn.nodeparam{ind}.fanOut;
                    jplProbs = 1.0;
                end
                classStratNode4b = simDoc.createElement('subParameter');
                classStratNode4b.setAttribute('classPath', 'jmt.engine.random.EmpiricalEntry');
                classStratNode4b.setAttribute('array', 'true');
                classStratNode4b.setAttribute('name', 'JobsPerLinkDis');
                for e=1:length(jplPoints)
                    classStratNode5b = simDoc.createElement('subParameter');
                    classStratNode5b.setAttribute('classPath', 'jmt.engine.random.EmpiricalEntry');
                    classStratNode5b.setAttribute('name', 'EmpiricalEntry');
                    classStratNode5bStation = simDoc.createElement('subParameter');
                    classStratNode5bStation.setAttribute('classPath', 'java.lang.String');
                    classStratNode5bStation.setAttribute('name', 'numbers');
                    classStratNode5bStationValueNode = simDoc.createElement('value');
                    classStratNode5bStationValueNode.appendChild(simDoc.createTextNode(int2str(jplPoints(e))));
                    classStratNode5bStation.appendChild(classStratNode5bStationValueNode);
                    classStratNode5b.appendChild(classStratNode5bStation);
                    classStratNode5bProbability = simDoc.createElement('subParameter');
                    classStratNode5bProbability.setAttribute('classPath', 'java.lang.Double');
                    classStratNode5bProbability.setAttribute('name', 'probability');
                    classStratNode5bProbabilityValueNode = simDoc.createElement('value');
                    classStratNode5bProbabilityValueNode.appendChild(simDoc.createTextNode(num2str(jplProbs(e),'%.15g')));
                    classStratNode5bProbability.appendChild(classStratNode5bProbabilityValueNode);
                    classStratNode5b.appendChild(classStratNode5bProbability);
                    classStratNode4b.appendChild(classStratNode5b);
                end
                classStratNode3.appendChild(classStratNode4b);

                classStratNode2.appendChild(classStratNode3);
            end
    end
    classStratNode.appendChild(classStratNode2);
    strategyNode.appendChild(classStratNode);
end
section.appendChild(strategyNode);
end
