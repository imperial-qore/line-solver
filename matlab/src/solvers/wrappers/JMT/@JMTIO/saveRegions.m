function [simElem, simDoc] = saveRegions(self, simElem, simDoc)
% [SIMELEM, SIMDOC] = SAVEREGIONS(SIMELEM, SIMDOC)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.getStruct;

% see _kb/06-solver-catalog.md (Wrappers: JMT per-class capacity export)
classCapCon = jmtClassCapCon(sn);
covered = false(sn.nstations,1); % stations already inside a region emitted below

% see _kb/06-solver-catalog.md (Wrappers: JMT per-class capacity export) for the blockingRegion XML shape

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

        % see _kb/06-solver-catalog.md (Wrappers: JMT per-class capacity export, LPS row)
        covered(ist) = true;
        for c = 1:sn.nclasses
            if isfinite(classCapCon(ist,c))
                jmtClassCapAssert(sn, ist, c); % rejects the non-loss cases first
                line_error(mfilename, sprintf(['Station %s has both LPS scheduling and a finite capacity %d for open class %s. ', ...
                    'JMT expresses both through a single blocking region, which admits only one drop rule per class, ', ...
                    'but LPS requires blocking while the open-class capacity requires dropping. Remove the per-class ', ...
                    'capacity or use a non-LPS scheduling strategy.'], nodeName, classCapCon(ist,c), sn.classnames{c}));
            end
        end

        for c = 1:sn.nclasses
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
% see _kb/06-solver-catalog.md (Wrappers: JMT blocking-region XML element order) for the required child-element order (steps 1-8 below)
for r=1:length(self.model.regions)
    blockingRegion = simDoc.createElement('blockingRegion');
    blockingRegion.setAttribute('name', ['FCRegion',num2str(r)]);
    blockingRegion.setAttribute('type', 'default');

    % see _kb/06-solver-catalog.md (Wrappers: JMT per-class capacity export, single/multi-node FCR rows)
    regionStations = [];
    for i=1:length(self.model.regions{r}.nodes)
        istr = sn.nodeToStation(self.model.regions{r}.nodes{i}.index);
        if istr > 0
            regionStations(end+1) = istr; %#ok<AGROW>
            covered(istr) = true;
        end
    end
    regionClassCap = Inf(1,sn.nclasses);
    for i=1:length(regionStations)
        istr = regionStations(i);
        for c=1:sn.nclasses
            if ~isfinite(classCapCon(istr,c))
                continue
            end
            if length(regionStations) > 1
                line_error(mfilename, sprintf(['Station %s carries a finite capacity %d for class %s and also belongs to ', ...
                    'the multi-station finite capacity region FCRegion%d. JMT constrains a blocking region as a whole and ', ...
                    'allows a node to belong to only one region, so a per-station class capacity cannot be expressed ', ...
                    'alongside it. Remove the per-class capacity or shrink the region to that single station.'], ...
                    sn.nodenames{sn.stationToNode(istr)}, classCapCon(istr,c), sn.classnames{c}, r));
            end
            jmtClassCapAssert(sn, istr, c);
            regionDrops = self.model.regions{r}.dropRule(c) == DropStrategy.DROP;
            if ~regionDrops
                line_error(mfilename, sprintf(['Station %s carries a finite capacity %d for class %s and also belongs to ', ...
                    'region FCRegion%d, whose drop rule for that class is "%s". A JMT blocking region admits a single drop ', ...
                    'rule per class, shared by all of its constraints, so the two cannot be merged: the per-class capacity ', ...
                    'of an open class is a loss constraint. Set the region drop rule for that class to DROP.'], ...
                    sn.nodenames{sn.stationToNode(istr)}, classCapCon(istr,c), sn.classnames{c}, r, ...
                    DropStrategy.toText(self.model.regions{r}.dropRule(c))));
            end
            regionClassCap(c) = classCapCon(istr,c);
        end
    end

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

    % 4. All classConstraint elements (for all classes; see _kb/06-solver-catalog.md
    % Wrappers: JMT blocking-region XML for why zero-population classes are included)
    for c=1:sn.nclasses
        % The two constraints apply to the same node set here, so the tighter
        % one subsumes the other and min() is exact.
        classMaxJobs = self.model.regions{r}.classMaxJobs(c);
        if classMaxJobs == Region.UNBOUNDED
            classMaxJobs = regionClassCap(c);
        else
            classMaxJobs = min(classMaxJobs, regionClassCap(c));
        end
        if isfinite(classMaxJobs) && classMaxJobs ~= Region.UNBOUNDED
            classConstraint = simDoc.createElement('classConstraint');
            classConstraint.setAttribute('jobClass', self.model.regions{r}.classes{c}.getName);
            classConstraint.setAttribute('maxJobsPerClass', num2str(classMaxJobs));
            blockingRegion.appendChild(classConstraint);
        end
    end

    % 5. All classMemoryConstraint elements (for all classes)
    for c=1:sn.nclasses
        if self.model.regions{r}.classMaxMemory(c) ~= Region.UNBOUNDED
            classMemoryConstraint = simDoc.createElement('classMemoryConstraint');
            classMemoryConstraint.setAttribute('jobClass', self.model.regions{r}.classes{c}.getName);
            classMemoryConstraint.setAttribute('maxMemoryPerClass', num2str(self.model.regions{r}.classMaxMemory(c)));
            blockingRegion.appendChild(classMemoryConstraint);
        end
    end

    % 6. All dropRules elements (for all classes)
    for c=1:sn.nclasses
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

    % 7. classWeight elements (drive the JMT 'FCR Capacity' measure, i.e. the
    % weighted occupation / Total Weight)
    for c=1:sn.nclasses
        if self.model.regions{r}.classWeight(c) ~= 1
            classWeightElem = simDoc.createElement('classWeight');
            classWeightElem.setAttribute('jobClass', self.model.regions{r}.classes{c}.getName);
            classWeightElem.setAttribute('weight', num2str(self.model.regions{r}.classWeight(c)));
            blockingRegion.appendChild(classWeightElem);
        end
    end

    % 8. All classSize elements (for all classes)
    for c=1:sn.nclasses
        if self.model.regions{r}.classSize(c) ~= 1
            classSizeElem = simDoc.createElement('classSize');
            classSizeElem.setAttribute('jobClass', self.model.regions{r}.classes{c}.getName);
            classSizeElem.setAttribute('size', num2str(self.model.regions{r}.classSize(c)));
            blockingRegion.appendChild(classSizeElem);
        end
    end

    simElem.appendChild(blockingRegion);
end

% see _kb/06-solver-catalog.md (Wrappers: JMT per-class capacity export) for the synthetic single-station ClassCapRegion export
classCapRegionIdx = 0;
for ist = 1:sn.nstations
    if covered(ist) || ~any(isfinite(classCapCon(ist,:)))
        continue
    end
    classCapRegionIdx = classCapRegionIdx + 1;
    ind = sn.stationToNode(ist);

    blockingRegion = simDoc.createElement('blockingRegion');
    blockingRegion.setAttribute('name', ['ClassCapRegion',num2str(classCapRegionIdx)]);
    blockingRegion.setAttribute('type', 'default');

    regionNode = simDoc.createElement('regionNode');
    regionNode.setAttribute('nodeName', sn.nodenames{ind});
    blockingRegion.appendChild(regionNode);

    globalConstraint = simDoc.createElement('globalConstraint');
    globalConstraint.setAttribute('maxJobs', '-1');
    blockingRegion.appendChild(globalConstraint);

    globalMemoryConstraint = simDoc.createElement('globalMemoryConstraint');
    globalMemoryConstraint.setAttribute('maxMemory', '-1');
    blockingRegion.appendChild(globalMemoryConstraint);

    % Validate every constrained class before emitting anything, so that an
    % inexpressible capacity errors out instead of half-writing a region.
    for c = 1:sn.nclasses
        if isfinite(classCapCon(ist,c))
            jmtClassCapAssert(sn, ist, c);
        end
    end

    for c = 1:sn.nclasses
        if isfinite(classCapCon(ist,c))
            classConstraint = simDoc.createElement('classConstraint');
            classConstraint.setAttribute('jobClass', sn.classnames{c});
            classConstraint.setAttribute('maxJobsPerClass', num2str(classCapCon(ist,c)));
            blockingRegion.appendChild(classConstraint);
        end
    end

    % see _kb/06-solver-catalog.md (Wrappers: JMT per-class capacity export) for the drop-rule semantics
    for c = 1:sn.nclasses
        if isfinite(classCapCon(ist,c))
            dropRuleElem = simDoc.createElement('dropRules');
            dropRuleElem.setAttribute('jobClass', sn.classnames{c});
            dropRuleElem.setAttribute('dropThisClass', 'true');
            blockingRegion.appendChild(dropRuleElem);
        end
    end

    simElem.appendChild(blockingRegion);
end
end

function classCapCon = jmtClassCapCon(sn)
% CLASSCAPCON = JMTCLASSCAPCON(SN)
%
% Per-class station capacities that must be exported as JMT blocking-region
% class constraints. classCapCon(i,r) is finite only when sn.classcap(i,r) is a
% genuine buffer limit that the rest of the exported model does not already
% imply; it is Inf otherwise, so that models without a per-class capacity emit
% byte-identical XML.
%
% refreshCapacity derives classcap(i,r) as
%   min(chain population of r, station.classCap(r), station.cap)
% so a value equal to the chain population or to the station total is not a
% per-class buffer at all: the closed population, respectively the Queue "size"
% written by saveBufferCapacity, already enforces it. Only a strictly tighter
% value carries information that JMT would otherwise lose.
classCapCon = Inf(sn.nstations, sn.nclasses);
if isempty(sn.classcap)
    return
end

% Population bound implied by the closed classes of each chain (Inf when open).
chainpop = Inf(1, sn.nclasses);
for c = 1:sn.nchains
    inchain = sn.inchain{c};
    chainpop(inchain) = sum(sn.njobs(inchain));
end

for ist = 1:sn.nstations
    ind = sn.stationToNode(ist);
    % A Source has no buffer, and a Place carries its per-class capacities in
    % the JMT Storage section, which does have a capacities array.
    if sn.nodetype(ind) == NodeType.Source || sn.nodetype(ind) == NodeType.Place
        continue
    end
    for r = 1:sn.nclasses
        % A class disabled at the station gets classcap 0 from refreshCapacity;
        % it is not routed there, and a zero constraint is not a buffer limit.
        if isnan(sn.rates(ist,r))
            continue
        end
        % Both Inf and intmax mean "unbounded": MATLAB leaves classcap at Inf,
        % while a struct marshalled from the JAR carries Integer.MAX_VALUE,
        % which jline.util.Utils.isInf also reads as infinite.
        if isfinite(sn.classcap(ist,r)) && sn.classcap(ist,r) < intmax && ...
                sn.classcap(ist,r) < min(sn.cap(ist), chainpop(r))
            classCapCon(ist,r) = sn.classcap(ist,r);
        end
    end
end
end

function jmtClassCapAssert(sn, ist, r)
% JMTCLASSCAPASSERT(SN, IST, R)
%
% Asserts that LINE's sn.classcap semantics for class r at station ist is
% expressible as a JMT blocking-region constraint, raising a descriptive error
% when it is not. Every constraint the caller emits is therefore a loss
% constraint, i.e. dropThisClass="true".
%
% LINE enforces classcap in State.afterEventStation by not enabling the
% arrival; sn.droprule is read there only to select BAS/BBS/RSRD blocking, so
% WAITQ and DROP behave identically for a capacity limit. The resulting
% reference semantics, verified against SolverCTMC, is
%   open class   -> the arrival is lost, upstream is unaffected
%   closed class -> the job stays at the upstream station and the population is
%                   conserved
% which is the same predicate the CTMC analyzer applies in canDropClass. Only
% the open case has a faithful blocking-region counterpart, dropThisClass=true:
% JMT discards the arrival at the region input station, exactly as LINE does.
%
% The closed case is rejected rather than mapped to dropThisClass=false. JMT
% does not hold a blocked job at the upstream station: it parks it in the
% region's synthetic input station, where it belongs to no station and leaves
% the upstream server free. Against SolverCTMC on a two-chain closed model that
% costs the population count (sum of queue lengths 4.79 instead of 5) and moves
% throughput by ~29%, so the region does not express this constraint at all.
if ~isinf(sn.njobs(r))
    line_error(mfilename, sprintf(['Station %s carries a finite capacity %d for closed class %s. LINE holds a blocked ', ...
        'closed job at its upstream station, whereas JMT can only express a per-class capacity as a blocking region, ', ...
        'which parks the job in the region input station instead, freeing the upstream server and losing it from the ', ...
        'population count. SolverJMT therefore cannot reproduce this model; use SolverCTMC, SolverSSA or SolverLDES, ', ...
        'or express the limit as the station capacity.'], sn.nodenames{sn.stationToNode(ist)}, ...
        sn.classcap(ist,r), sn.classnames{r}));
end
% The blocking strategies below have no blocking-region counterpart either: a
% region drops or defers the arrival, it cannot hold a job in the upstream
% server.
switch sn.droprule(ist,r)
    case {DropStrategy.BAS, DropStrategy.BBS, DropStrategy.RSRD, ...
          DropStrategy.RETRIAL, DropStrategy.RETRIAL_WITH_LIMIT}
        line_error(mfilename, sprintf(['Station %s applies drop strategy "%s" to class %s and also carries a finite ', ...
            'capacity %d for it. JMT exports a per-class capacity as a blocking region, which can only drop or defer an ', ...
            'arrival and cannot reproduce that strategy. Remove the per-class capacity or use the station capacity ', ...
            'instead, which is exported with its drop strategy.'], sn.nodenames{sn.stationToNode(ist)}, ...
            DropStrategy.toText(sn.droprule(ist,r)), sn.classnames{r}, sn.classcap(ist,r)));
end
end
