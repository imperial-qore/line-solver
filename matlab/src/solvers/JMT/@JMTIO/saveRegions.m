function [simElem, simDoc] = saveRegions(self, simElem, simDoc)
% [SIMELEM, SIMDOC] = SAVEREGIONS(SIMELEM, SIMDOC)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.getStruct;

% <blockingRegion name="FCRegion1" type="default">
% <regionNode nodeName="Queue 1"/>
% <regionNode nodeName="Queue 2"/>
% <globalConstraint maxJobs="2"/>
% <globalMemoryConstraint maxMemory="-1"/>
% <classConstraint jobClass="Class1" maxJobsPerClass="-1"/>
% <classMemoryConstraint jobClass="Class1" maxMemoryPerClass="-1"/>
% <dropRules dropThisClass="false" jobClass="Class1"/>
% <classWeight jobClass="Class1" weight="1"/>
% <classSize jobClass="Class1" size="1"/>
% </blockingRegion>

% First, create implicit FCR regions for LPS queues
% LPS uses FCR to limit concurrent jobs at PS station
lpsRegionIdx = length(self.model.regions);  % Start numbering after explicit regions
for i = 1:sn.nstations
    ist = i;
    if sn.sched(ist) == SchedStrategy.LPS
        lpsRegionIdx = lpsRegionIdx + 1;
        lpsLimit = sn.schedparam(ist, 1);  % LPS limit stored in first column
        ind = sn.stationToNode(ist);
        nodeName = sn.nodenames{ind};

        blockingRegion = simDoc.createElement('blockingRegion');
        blockingRegion.setAttribute('name', ['LPSRegion', num2str(lpsRegionIdx)]);
        blockingRegion.setAttribute('type', 'default');

        regionNode = simDoc.createElement('regionNode');
        regionNode.setAttribute('nodeName', nodeName);
        blockingRegion.appendChild(regionNode);

        globalConstraint = simDoc.createElement('globalConstraint');
        globalConstraint.setAttribute('maxJobs', num2str(lpsLimit));
        blockingRegion.appendChild(globalConstraint);

        globalMemoryConstraint = simDoc.createElement('globalMemoryConstraint');
        globalMemoryConstraint.setAttribute('maxMemory', '-1');
        blockingRegion.appendChild(globalMemoryConstraint);

        for c = 1:sn.nclasses
            % Skip closed classes with 0 customers (JMT cannot handle these)
            if isfinite(sn.njobs(c)) && sn.njobs(c) == 0
                continue;
            end

            % dropRules - LPS uses blocking (waitq), not drop
            dropRuleElem = simDoc.createElement('dropRules');
            dropRuleElem.setAttribute('jobClass', sn.classnames{c});
            dropRuleElem.setAttribute('dropThisClass', 'false');
            blockingRegion.appendChild(dropRuleElem);
        end

        simElem.appendChild(blockingRegion);
    end
end

% Now save explicit user-defined regions
% JMT XML schema requires elements in specific order:
% regionNode*, globalConstraint, globalMemoryConstraint,
% classConstraint*, classMemoryConstraint*, dropRules*, classWeight*, classSize*, classSoftDeadline*
for r=1:length(self.model.regions)
    blockingRegion = simDoc.createElement('blockingRegion');
    blockingRegion.setAttribute('name', ['FCRegion',num2str(r)]);
    blockingRegion.setAttribute('type', 'default');

    % 1. regionNode elements
    for i=1:length(self.model.regions{r}.nodes)
        regionNode = simDoc.createElement('regionNode');
        regionNode.setAttribute('nodeName', self.model.regions{r}.nodes{i}.getName);
        blockingRegion.appendChild(regionNode);
    end

    % 2. globalConstraint
    globalConstraint = simDoc.createElement('globalConstraint');
    globalConstraint.setAttribute('maxJobs', num2str(self.model.regions{r}.globalMaxJobs));
    blockingRegion.appendChild(globalConstraint);

    % 3. globalMemoryConstraint
    globalMemoryConstraint = simDoc.createElement('globalMemoryConstraint');
    globalMemoryConstraint.setAttribute('maxMemory', num2str(self.model.regions{r}.globalMaxMemory));
    blockingRegion.appendChild(globalMemoryConstraint);

    % 4. All classConstraint elements (for all classes)
    for c=1:sn.nclasses
        if isfinite(sn.njobs(c)) && sn.njobs(c) == 0
            continue;
        end
        if self.model.regions{r}.classMaxJobs(c) ~= Region.UNBOUNDED
            classConstraint = simDoc.createElement('classConstraint');
            classConstraint.setAttribute('jobClass', self.model.regions{r}.classes{c}.getName);
            classConstraint.setAttribute('maxJobsPerClass', num2str(self.model.regions{r}.classMaxJobs(c)));
            blockingRegion.appendChild(classConstraint);
        end
    end

    % 5. All classMemoryConstraint elements (for all classes)
    for c=1:sn.nclasses
        if isfinite(sn.njobs(c)) && sn.njobs(c) == 0
            continue;
        end
        if self.model.regions{r}.classMaxMemory(c) ~= Region.UNBOUNDED
            classMemoryConstraint = simDoc.createElement('classMemoryConstraint');
            classMemoryConstraint.setAttribute('jobClass', self.model.regions{r}.classes{c}.getName);
            classMemoryConstraint.setAttribute('maxMemoryPerClass', num2str(self.model.regions{r}.classMaxMemory(c)));
            blockingRegion.appendChild(classMemoryConstraint);
        end
    end

    % 6. All dropRules elements (for all classes)
    for c=1:sn.nclasses
        if isfinite(sn.njobs(c)) && sn.njobs(c) == 0
            continue;
        end
        % Always write dropRules element - JMT defaults to drop when not specified
        dropRuleElem = simDoc.createElement('dropRules');
        dropRuleElem.setAttribute('jobClass', self.model.regions{r}.classes{c}.getName);
        if self.model.regions{r}.dropRule(c) == DropStrategy.DROP
            dropRuleElem.setAttribute('dropThisClass', 'true');
        else
            dropRuleElem.setAttribute('dropThisClass', 'false');
        end
        blockingRegion.appendChild(dropRuleElem);
    end

    % 7. classWeight elements (currently disabled in JMT)
    % for c=1:sn.nclasses
    %     if isfinite(sn.njobs(c)) && sn.njobs(c) == 0
    %         continue;
    %     end
    %     if self.model.regions{r}.classWeight(c) ~= 1
    %         classWeightElem = simDoc.createElement('classWeight');
    %         classWeightElem.setAttribute('jobClass', self.model.regions{r}.classes{c}.getName);
    %         classWeightElem.setAttribute('weight', num2str(self.model.regions{r}.classWeight(c)));
    %         blockingRegion.appendChild(classWeightElem);
    %     end
    % end

    % 8. All classSize elements (for all classes)
    for c=1:sn.nclasses
        if isfinite(sn.njobs(c)) && sn.njobs(c) == 0
            continue;
        end
        if self.model.regions{r}.classSize(c) ~= 1
            classSizeElem = simDoc.createElement('classSize');
            classSizeElem.setAttribute('jobClass', self.model.regions{r}.classes{c}.getName);
            classSizeElem.setAttribute('size', num2str(self.model.regions{r}.classSize(c)));
            blockingRegion.appendChild(classSizeElem);
        end
    end

    simElem.appendChild(blockingRegion);
end
end
