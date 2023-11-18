function [simElem, simDoc] = saveRegions(self, simElem, simDoc)
% [SIMELEM, SIMDOC] = SAVEREGIONS(SIMELEM, SIMDOC)

% Copyright (c) 2012-2023, Imperial College London
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

for r=1:length(self.model.regions)
    blockingRegion = simDoc.createElement('blockingRegion');
    blockingRegion.setAttribute('name', ['FCRegion',num2str(r)]);
    blockingRegion.setAttribute('type', 'default');
    for i=1:length(self.model.regions{r}.nodes)
        regionNode = simDoc.createElement('regionNode');
        regionNode.setAttribute('nodeName', self.model.regions{r}.nodes{i}.getName);
        blockingRegion.appendChild(regionNode);
    end
    globalConstraint = simDoc.createElement('globalConstraint');
    globalConstraint.setAttribute('maxJobs', num2str(self.model.regions{r}.globalMaxJobs));
    blockingRegion.appendChild(globalConstraint);
    globalMemoryConstraint = simDoc.createElement('globalMemoryConstraint');
    globalMemoryConstraint.setAttribute('maxMemory', num2str(self.model.regions{r}.globalMaxMemory));
    blockingRegion.appendChild(globalMemoryConstraint);
    for c=1:sn.nclasses
        if self.model.regions{r}.classMaxJobs(c) ~= FiniteCapacityRegion.UNBOUNDED
            classConstraint = simDoc.createElement('classConstraint');
            classConstraint.setAttribute('jobClass', self.model.regions{r}.classes{c}.getName);
            classConstraint.setAttribute('maxJobsPerClass', num2str(self.model.regions{r}.classMaxJobs(c)));
            blockingRegion.appendChild(classConstraint);
        end
        if self.model.regions{r}.classMaxMemory(c) ~= FiniteCapacityRegion.UNBOUNDED
            classMemoryConstraint = simDoc.createElement('classMemoryConstraint');
            classMemoryConstraint.setAttribute('jobClass', self.model.regions{r}.classes{c}.getName);
            classMemoryConstraint.setAttribute('maxMemoryPerClass', num2str(self.model.regions{r}.classMaxMemory(c)));
            blockingRegion.appendChild(classMemoryConstraint);
        end
        if self.model.regions{r}.dropRule(c)
            dropRule = simDoc.createElement('dropRules');
            dropRule.setAttribute('jobClass', self.model.regions{r}.classes{c}.getName);
            dropRule.setAttribute('dropThisClass', 'true');
            blockingRegion.appendChild(dropRule);
        end
        %         % Weight disabled in JMT
        %         if self.model.regions{r}.classWeight(c) ~= 1
        %             classMemoryConstraint = simDoc.createElement('classWeight');
        %             classMemoryConstraint.setAttribute('jobClass', self.model.regions{r}.classes{c}.getName);
        %             classMemoryConstraint.setAttribute('weight', num2str(self.model.regions{r}.classWeight(c)));
        %             blockingRegion.appendChild(classMemoryConstraint);
        %         end
        if self.model.regions{r}.classSize(c) ~= 1
            classMemoryConstraint = simDoc.createElement('classSize');
            classMemoryConstraint.setAttribute('jobClass', self.model.regions{r}.classes{c}.getName);
            classMemoryConstraint.setAttribute('size', num2str(self.model.regions{r}.classSize(c)));
            blockingRegion.appendChild(classMemoryConstraint);
        end
    end
    simElem.appendChild(blockingRegion);
end
end
