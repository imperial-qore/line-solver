function linemodel_save(model, filename)
% LINEMODEL_SAVE Save a LINE model to JSON.
%
%   LINEMODEL_SAVE(MODEL, FILENAME) saves the model to the specified JSON
%   file, conforming to the line-model.schema.json specification.
%
% Parameters:
%   model    - Network, LayeredNetwork, Workflow, or Environment object
%   filename - output file path (should end in .json)
%
% Example:
%   model = Network('M/M/1');
%   source = Source(model, 'Source');
%   queue  = Queue(model, 'Queue', SchedStrategy.FCFS);
%   sink   = Sink(model, 'Sink');
%   oclass = OpenClass(model, 'Class1');
%   source.setArrival(oclass, Exp(1.0));
%   queue.setService(oclass, Exp(2.0));
%   P = model.initRoutingMatrix();
%   P{1}(1,2) = 1; P{1}(2,3) = 1;
%   model.link(P);
%   linemodel_save(model, 'mm1.json');
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isa(model, 'LayeredNetwork')
    modelMap = layered2json(model);
elseif isa(model, 'Workflow')
    modelMap = workflow2json(model);
elseif isa(model, 'Environment')
    modelMap = environment2json(model);
else
    modelMap = network2json(model);
end

% Build the full document
sb = {};
sb{end+1} = '{';
sb{end+1} = '  "format": "line-model",';
sb{end+1} = '  "version": "1.0",';
sb{end+1} = ['  "model": ', encode_value(modelMap, 2)];
sb{end+1} = '}';
jsonStr = strjoin(sb, newline);

fid = fopen(filename, 'w');
if fid == -1
    error('linemodel_save:fileOpen', 'Cannot open file: %s', filename);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonStr);
end


% =========================================================================
%  Network serialization
% =========================================================================

function result = network2json(model)
% Convert a Network to a dictionary (the writer sorts keys, as containers.Map did)

result = configureDictionary('string','cell');
result{'type'} = 'Network';
result{'name'} = model.getName();

nodes = model.getNodes();
classes = model.getClasses();
K = length(classes);
M = length(nodes);

% A STATE ON A STRICT SUBSET OF THE STATEFUL NODES DOES NOT TRAVEL, because it
% is not an initialization here either: getState runs initDefault whenever
% hasInitState is false, so this side answers for the DEFAULT marking while a
% document naming the one node the caller moved made the reader combine that row
% with default markings for the rest. On Think -> Q1 with 2 jobs and Q1 alone set
% to 2, the C++ read (Think=2, Q1=2) and answered getProbSysAggr 0 for a joint
% state holding 4 of 2 jobs, against the 0.4 this side reports.
%
% A PAS PLACEMENT IS THE EXCEPTION, and it is initDefault.m's own: its `hasUser`
% branch keeps the user row for a pass-and-swap station while defaulting every
% other, the ordering being a required input there rather than a default.
stateTravels = true;
if ismethod(model, 'hasInitState')
    stateTravels = logical(model.hasInitState());
end

% --- Nodes ---
nodesJson = {};
for i = 1:M
    node = nodes{i};
    nodeStateTravels = stateTravels;
    if ~nodeStateTravels && isa(node, 'Station') && isprop(node, 'schedStrategy')
        nodeStateTravels = node.schedStrategy == SchedStrategy.PAS;
    end

    % Skip implicit ClassSwitch nodes (auto-created by link())
    if isa(node, 'ClassSwitch') && isprop(node, 'autoAdded') && node.autoAdded
        continue;
    end

    nj = configureDictionary('string','cell');
    nj{'name'} = node.name;
    nj{'type'} = node_type_str(node);

    % A Logger is defined by the file it writes: a Logger with no file name
    % exports as a LogTunnel that logs nowhere, so the solver runs and the
    % trace the user asked for is silently absent.
    if isa(node, 'Logger') && ~isempty(node.fileName)
        nj{'fileName'} = node.fileName;
    end

    % Place extends Station not Queue, so isa(node,'Queue') alone would drop its scheduling strategy
    if isa(node, 'Delay')
        nj{'scheduling'} = 'INF';
    elseif isa(node, 'Place')
        if node.isQueueing() && ~isempty(node.schedStrategy)
            nj{'scheduling'} = sched_id_to_str(node.schedStrategy);
        end
    elseif isa(node, 'Queue')
        sched = node.schedStrategy;
        if ~isempty(sched)
            nj{'scheduling'} = sched_id_to_str(sched);
        end
    end

    % Infinite server count deliberately not emitted -- see _kb/09-ldes-and-cache.md
    if isa(node, 'Station') && ~isa(node, 'Delay') && ...
            (isa(node, 'Queue') || (isa(node, 'Place') && node.isQueueing()))
        ns = node.numberOfServers;
        if isfinite(ns) && ns > 1
            nj{'servers'} = ns;
        end
    end

    % Station, not Queue: Place extends Station directly (see scheduling note above)
    if isa(node, 'Station')
        c = node.cap;
        if ~isempty(c) && isfinite(c) && c > 0
            nj{'buffer'} = c;
        end
    end

    % Per-class buffer capacity
    if isa(node, 'Station') && ~isempty(node.classCap)
        ccMap = configureDictionary('string','cell');
        for r = 1:K
            jc = classes{r};
            if r <= length(node.classCap) && isfinite(node.classCap(r))
                ccMap{jc.name} = node.classCap(r);
            end
        end
        if numEntries(ccMap) > 0
            nj{'classCap'} = ccMap;
        end
    end

    % Drop rules
    if isa(node, 'Station') && ~isempty(node.dropRule)
        drMap = configureDictionary('string','cell');
        for r = 1:K
            jc = classes{r};
            if r <= length(node.dropRule)
                dr = node.dropRule(r);
                drStr = droprule_to_str(dr);
                if ~isempty(drStr)
                    drMap{jc.name} = drStr;
                end
            end
        end
        if numEntries(drMap) > 0
            nj{'dropRule'} = drMap;
        end
    end

    % Load-dependent scaling
    if isa(node, 'Station') && ~isempty(node.lldScaling)
        ldMap = configureDictionary('string','cell');
        ldMap{'type'} = 'loadDependent';
        ldMap{'scaling'} = node.lldScaling(:)';
        nj{'loadDependence'} = ldMap;
    end

    % beta_{i,r}(n) handle materialized over the per-class box lattice (cannot cross JSON as a handle)
    if isa(node, 'Station') && ~isempty(node.lcdScaling)
        maxc = zeros(1, K);
        for r = 1:K
            if isa(classes{r}, 'ClosedClass') && isfinite(classes{r}.population)
                maxc(r) = round(classes{r}.population);
            else
                maxc(r) = 10;   % open-class saturation cutoff (beta clamped beyond)
            end
        end
        cdMap = configureDictionary('string','cell');
        cdMap{'type'} = 'classDependent';
        % num2cell keeps a single-class (1x1) vector from collapsing to a scalar.
        cdMap{'cutoffs'} = num2cell(double(maxc(:)'));
        cdMap{'scaling'} = cd_scaling_table(node.lcdScaling, maxc, K);
        % Declared peak rate scaling per class (Util = T*S/peak). Broadcast a
        % scalar to K entries so the reader restores a per-class vector.
        pk = node.lcdScalingPeak;
        if isscalar(pk)
            pk = repmat(pk, 1, K);
        end
        cdMap{'peak'} = num2cell(double(pk(:)'));
        nj{'classDependence'} = cdMap;
    end

    % eta_i(n) joint-dependence handle (non-product-form), materialized over the
    % same per-class box lattice. Wire key "jointDependence", twin of the
    % classDependence block above; matches the JAR/Python readers.
    if isa(node, 'Station') && ~isempty(node.ljdScaling)
        maxc = zeros(1, K);
        for r = 1:K
            if isa(classes{r}, 'ClosedClass') && isfinite(classes{r}.population)
                maxc(r) = round(classes{r}.population);
            else
                maxc(r) = 10;   % open-class saturation cutoff (eta clamped beyond)
            end
        end
        jdMap = configureDictionary('string','cell');
        jdMap{'type'} = 'jointDependent';
        jdMap{'cutoffs'} = num2cell(double(maxc(:)'));
        jdMap{'scaling'} = cd_scaling_table(node.ljdScaling, maxc, K);
        pk = node.ljdScalingPeak;
        if isscalar(pk)
            pk = repmat(pk, 1, K);
        end
        jdMap{'peak'} = num2cell(double(pk(:)'));
        nj{'jointDependence'} = jdMap;
    end

    % A queueing Place carries per-class service processes plus departure discipline; an ordinary Place has neither
    svc = configureDictionary('string','cell');
    depDisc = configureDictionary('string','cell');
    for r = 1:K
        jc = classes{r};
        dist = [];
        if isa(node, 'Source')
            try
                dist = node.getArrivalProcess(jc);
            catch
                dist = [];
            end
        elseif isa(node, 'Place')
            if node.isQueueing() && numel(node.serviceProcess) >= jc.index ...
                    && ~isempty(node.serviceProcess{jc.index})
                dist = node.serviceProcess{jc.index};
            end
        elseif isa(node, 'Queue') || isa(node, 'Delay')
            % An OI/PAS queue carries mu(c) below as oiServiceRate, but it ALSO
            % keeps the representative per-class Exp(mu([r])) that setService
            % derives, and the reader needs it: suppressing it left the station
            % with no service process for any class, and line-cli answered
            % pas_mmk with an EMPTY AvgTable. The native writer emits both.
            try
                dist = node.getService(jc);
            catch
                dist = [];
            end
        end
        if ~isempty(dist) && ~isa(dist, 'Disabled')
            dj = dist2json(dist);
            if ~isempty(dj)
                svc{jc.name} = dj;
                if isa(node, 'Place') && numel(node.departureDiscipline) >= jc.index
                    depDisc{jc.name} = depdisc_to_str(node.departureDiscipline(jc.index));
                end
            end
        end
    end
    if numEntries(svc) > 0
        nj{'service'} = svc;
    end
    if numEntries(depDisc) > 0
        nj{'departureDiscipline'} = depDisc;
    end

    % OI/PAS mu(c) serialized as a macrostate table keyed by per-class counts -- see _kb/09-ldes-and-cache.md
    if isa(node, 'Queue') && ~isempty(node.svcRateFun)
        maxc = zeros(1, K);
        for r = 1:K
            if isa(classes{r}, 'ClosedClass') && isfinite(classes{r}.population)
                maxc(r) = round(classes{r}.population);
            else
                maxc(r) = 10;   % open-class saturation cutoff (mu constant beyond)
            end
        end
        rateMap = oi_rate_table(node.svcRateFun, maxc, K);
        nj{'oiServiceRate'} = rateMap;
        % num2cell keeps jsonencode from collapsing a single-class (1x1) vector
        % to a scalar, which the JAR reader parses as an array.
        nj{'oiCutoffs'} = num2cell(double(maxc(:)'));
        if node.schedStrategy == SchedStrategy.PAS && ~isempty(node.swapGraph)
            sgRows = cell(1, size(node.swapGraph, 1));
            for sgi = 1:size(node.swapGraph, 1)
                sgRows{sgi} = num2cell(double(node.swapGraph(sgi, :)));
            end
            nj{'swapGraph'} = sgRows;
        end
    end

    % Batch arrivals: the batch-size law released at each arrival epoch, per
    % class. Separate from 'service' above, which only spaces the epochs.
    if isa(node, 'Source') && ~isempty(node.arrivalBatch)
        batchMap = configureDictionary('string','cell');
        for r = 1:K
            if numel(node.arrivalBatch) >= r && ~isempty(node.arrivalBatch{r})
                bj = dist2json(node.arrivalBatch{r});
                if ~isempty(bj)
                    batchMap{classes{r}.name} = bj;
                end
            end
        end
        if numEntries(batchMap) > 0
            nj{'arrivalBatch'} = batchMap;
        end
    end

    % Marked (MMAP) arrival binding: class names ordered by mark
    if isa(node, 'Source') && ~isempty(node.markedClasses)
        markedNames = cell(1, numel(node.markedClasses));
        for km = 1:numel(node.markedClasses)
            markedNames{km} = classes{node.markedClasses(km)}.name;
        end
        nj{'markedClasses'} = markedNames;
    end

    % ClassSwitch matrix
    if isa(node, 'ClassSwitch')
        csm = node.server.csMatrix;
        if ~isempty(csm)
            csDict = configureDictionary('string','cell');
            for ri = 1:K
                row = configureDictionary('string','cell');
                for ci = 1:K
                    if ri <= size(csm,1) && ci <= size(csm,2) && csm(ri,ci) ~= 0
                        row{classes{ci}.name} = csm(ri,ci);
                    end
                end
                if numEntries(row) > 0
                    csDict{classes{ri}.name} = row;
                end
            end
            if numEntries(csDict) > 0
                nj{'classSwitchMatrix'} = csDict;
            end
        end
    end

    % Cache config
    if isa(node, 'Cache')
        cc = configureDictionary('string','cell');
        cc{'items'} = node.items.nitems;
        ilc = node.itemLevelCap;
        % Replacement policy emitted verbatim; CLIMB is rewritten at solve time, never here -- see _kb/09-ldes-and-cache.md
        if isscalar(ilc)
            cc{'capacity'} = ilc;
        else
            cc{'capacity'} = ilc(:)';
        end
        cc{'replacement'} = repl_to_str(node.replacestrategy);
        if isprop(node,'admissionProb') && ~isempty(node.admissionProb)
            cc{'admissionProb'} = node.admissionProb;
        end
        % per-item storage costs and per-list cost caps (ton21cache Sec. IX)
        if isprop(node,'itemSize') && ~isempty(node.itemSize)
            cc{'itemSizes'} = node.itemSize(:)';
        end
        if isprop(node,'costCap') && ~isempty(node.costCap)
            if isprop(node,'costCapGlobal') && node.costCapGlobal
                cc{'costCaps'} = node.costCap(1);
            else
                cc{'costCaps'} = node.costCap(:)';
            end
        end

        % Hit/miss class mappings
        hc = full(node.server.hitClass);
        mc = full(node.server.missClass);
        if ~isempty(hc) && any(hc > 0)
            hitMap = configureDictionary('string','cell');
            for hi = 1:length(hc)
                if hc(hi) > 0 && hi <= K && hc(hi) <= K
                    hitMap{classes{hi}.name} = classes{hc(hi)}.name;
                end
            end
            if numEntries(hitMap) > 0
                cc{'hitClass'} = hitMap;
            end
        end
        if ~isempty(mc) && any(mc > 0)
            missMap = configureDictionary('string','cell');
            for mi = 1:length(mc)
                if mc(mi) > 0 && mi <= K && mc(mi) <= K
                    missMap{classes{mi}.name} = classes{mc(mi)}.name;
                end
            end
            if numEntries(missMap) > 0
                cc{'missClass'} = missMap;
            end
        end

        % Read popularity distributions (setRead). popularity is keyed
        % (itemSetIndex, class), so emit only this cache's own row: with more than
        % one cache the other rows belong to other item sets.
        if ~isempty(node.popularity)
            popMap = configureDictionary('string','cell');
            pi = node.items.index;
            if pi <= size(node.popularity, 1)
                for pj = 1:size(node.popularity, 2)
                    if ~isempty(node.popularity{pi, pj})
                        popDist = node.popularity{pi, pj};
                        dj = dist2json(popDist);
                        if ~isempty(dj) && pj <= K
                            popMap{classes{pj}.name} = dj;
                        end
                    end
                end
            end
            if numEntries(popMap) > 0
                cc{'popularity'} = popMap;
            end
        end

        % Per-item classes of a cache network (Cache.setItemReadClasses): the item
        % each such class reads, so the mapping survives the JSON round trip.
        if isprop(node,'itemOfClass') && ~isempty(node.itemOfClass)
            itemClsMap = configureDictionary('string','cell');
            for pj = 1:min(numel(node.itemOfClass), K)
                if node.itemOfClass(pj) > 0
                    itemClsMap{classes{pj}.name} = node.itemOfClass(pj);
                end
            end
            if numEntries(itemClsMap) > 0
                cc{'itemClass'} = itemClsMap;
            end
        end

        % Access-cost graph shared by all classes, or full per-class accessProb; default super-diagonal rebuilt on load
        if ~isempty(node.graph)
            gArr = cell(1, numel(node.graph));
            for gi = 1:numel(node.graph)
                gArr{gi} = full(node.graph{gi});
            end
            cc{'accessGraph'} = gArr;
        elseif ~isempty(node.accessProb)
            [Kap, Nap] = size(node.accessProb);
            apArr = cell(1, Kap);
            for k1 = 1:Kap
                rowArr = cell(1, Nap);
                for k2 = 1:Nap
                    if ~isempty(node.accessProb{k1, k2})
                        rowArr{k2} = full(node.accessProb{k1, k2});
                    else
                        rowArr{k2} = [];
                    end
                end
                apArr{k1} = rowArr;
            end
            cc{'accessProb'} = apArr;
        end

        % Initial cache state [class counts | contents | retrieval bitmap]
        cacheState = node.getState;
        if ~isempty(cacheState)
            cc{'initialState'} = num2cell(double(full(cacheState(1, :))));
        end

        nj{'cache'} = cc;

        % Flat cache fields for the Java LineModelIO reader, mirrored from cc so both forms agree -- see _kb/09-ldes-and-cache.md
        nj{'numItems'} = double(node.items.nitems);
        nj{'itemLevelCap'} = num2cell(double(ilc(:)'));
        nj{'replacementStrategy'} = cc{'replacement'};
        flatKeys = {'hitClass', 'missClass', 'popularity', 'accessGraph', ...
                    'accessProb', 'initialState', 'admissionProb', 'itemClass'};
        for fk = 1:numel(flatKeys)
            if isKey(cc, flatKeys{fk})
                nj{flatKeys{fk}} = cc{flatKeys{fk}};
            end
        end

        % Retrieval system flat block (setRetrievalSystem only) -- see _kb/09-ldes-and-cache.md
        if ~isempty(node.retrievalSystemCapacity) && node.retrievalSystemCapacity > 0
            byClass = configureDictionary('string','cell');
            nItemsR = node.items.nitems;
            rc = node.server.retrievalClasses;
            qKeys = keys(node.retrievalSystemQueueIndices);
            for kk = 1:numel(qKeys)
                key0 = qKeys(kk);                 % jobinClass.index - 1 (0-based)
                inIdx = double(key0) + 1;
                if inIdx < 1 || inIdx > K
                    continue;
                end
                entry = configureDictionary('string','cell');
                qidxs = node.retrievalSystemQueueIndices{key0};
                qnames = cell(1, numel(qidxs));
                for qi = 1:numel(qidxs)
                    qnames{qi} = nodes{qidxs(qi)}.name;
                end
                entry{'queues'} = qnames;
                itemsMap = configureDictionary('string','cell');
                for it = 1:nItemsR
                    if size(rc, 1) >= it && size(rc, 2) >= inIdx
                        rClassIdx = rc(it, inIdx);
                        if rClassIdx > 0 && rClassIdx <= K
                            itemsMap{num2str(it - 1)} = classes{rClassIdx}.name;
                        end
                    end
                end
                if numEntries(itemsMap) > 0
                    entry{'items'} = itemsMap;
                end
                byClass{classes{inIdx}.name} = entry;
            end
            if numEntries(byClass) > 0
                rsRoot = configureDictionary('string','cell');
                rsRoot{'capacity'} = double(node.retrievalSystemCapacity);
                rsRoot{'byClass'} = byClass;
                nj{'retrievalSystem'} = rsRoot;
            end
        end
    end

    % Fork tasksPerLink and the variable-forking-level overrides. Each override
    % list is emitted only when non-empty, so a plain fork's JSON is byte-for-
    % byte what it was.
    if isa(node, 'Fork')
        if ~isempty(node.output) && isprop(node.output, 'tasksPerLink') && node.output.tasksPerLink > 1
            nj{'tasksPerLink'} = node.output.tasksPerLink;
        end
        if ~isempty(node.output) && isprop(node.output, 'tasksPerLinkByDest') && ~isempty(node.output.tasksPerLinkByDest)
            ovs = {};
            for e = 1:length(node.output.tasksPerLinkByDest)
                ov = node.output.tasksPerLinkByDest(e);
                d = configureDictionary('string','cell');
                d{'dest'} = ov.dest;
                d{'class'} = double(ov.class);
                d{'value'} = double(ov.value);
                ovs{end+1} = d; %#ok<AGROW>
            end
            nj{'fanOutByDest'} = ovs;
        end
        if ~isempty(node.output) && isprop(node.output, 'tasksPerLinkDist') && ~isempty(node.output.tasksPerLinkDist)
            ovs = {};
            for e = 1:length(node.output.tasksPerLinkDist)
                ov = node.output.tasksPerLinkDist(e);
                d = configureDictionary('string','cell');
                d{'dest'} = ov.dest;
                d{'class'} = double(ov.class);
                % pmf and support, the two parameters a DiscreteSampler is
                d{'p'} = double(ov.dist.getParam(1).paramValue(:)');
                d{'x'} = double(ov.dist.getParam(2).paramValue(:)');
                ovs{end+1} = d; %#ok<AGROW>
            end
            nj{'fanOutDist'} = ovs;
        end
        if ~isempty(node.output) && isprop(node.output, 'branchProb') && ~isempty(node.output.branchProb)
            ovs = {};
            for e = 1:length(node.output.branchProb)
                ov = node.output.branchProb(e);
                d = configureDictionary('string','cell');
                d{'dest'} = ov.dest;
                d{'class'} = double(ov.class);
                d{'value'} = double(ov.value);
                ovs{end+1} = d; %#ok<AGROW>
            end
            nj{'fanOutProb'} = ovs;
        end
    end

    % Join paired fork and join strategy
    if isa(node, 'Join')
        if ~isempty(node.joinOf)
            nj{'forkNode'} = node.joinOf.name;
        end
        % Serialize per-class join strategy if non-default
        if ~isempty(node.input) && isprop(node.input, 'joinStrategy') && ~isempty(node.input.joinStrategy)
            for r = 1:K
                jc = classes{r};
                if r <= length(node.input.joinStrategy) && ~isempty(node.input.joinStrategy{r})
                    js = node.input.joinStrategy{r};
                    if js ~= JoinStrategy.STD
                        if js == JoinStrategy.PARTIAL
                            nj{'joinStrategy'} = 'PARTIAL';
                        end
                    end
                end
            end
        end
        if ~isempty(node.input) && isprop(node.input, 'joinRequired') && ~isempty(node.input.joinRequired)
            for r = 1:K
                jc = classes{r};
                if r <= length(node.input.joinRequired) && ~isempty(node.input.joinRequired{r})
                    jq = node.input.joinRequired{r};
                    if jq > 0
                        nj{'joinQuorum'} = jq;
                    end
                end
            end
        end
    end

    % DPS/GPS weights, including DPSPRIO/GPSPRIO (same schedStrategyPar); omitting them resets weights to 1 on reload
    if isa(node, 'Queue') && ~isa(node, 'Delay')
        sched = node.schedStrategy;
        if ~isempty(sched) && (sched == SchedStrategy.DPS || sched == SchedStrategy.GPS || ...
                sched == SchedStrategy.DPSPRIO || sched == SchedStrategy.GPSPRIO)
            sp = configureDictionary('string','cell');
            for r = 1:K
                jc = classes{r};
                try
                    w = node.schedStrategyPar(r);
                    if ~isempty(w) && isfinite(w) && w > 0
                        sp{jc.name} = w;
                    end
                catch
                end
            end
            if numEntries(sp) > 0
                nj{'schedParams'} = sp;
            end
        end
    end

    % Transition modes
    if isa(node, 'Transition')
        modesJson = {};
        nModes = node.getNumberOfModes();
        allNodes = model.getNodes();
        for mi = 1:nModes
            mj = configureDictionary('string','cell');
            if mi <= length(node.modeNames) && ~isempty(node.modeNames{mi})
                mj{'name'} = node.modeNames{mi};
            else
                mj{'name'} = sprintf('Mode%d', mi);
            end
            % Immediate mode's Exp(1) placeholder distribution is omitted, else a reader sees a timed mode at rate 1
            isImmediateMode = mi <= length(node.timingStrategies) && ...
                node.timingStrategies(mi) == TimingStrategy.IMMEDIATE;
            if ~isImmediateMode && mi <= length(node.distributions) && ...
                    ~isempty(node.distributions{mi})
                dj = dist2json(node.distributions{mi});
                if ~isempty(dj)
                    mj{'distribution'} = dj;
                end
            end
            % Timing strategy
            if mi <= length(node.timingStrategies)
                if node.timingStrategies(mi) == TimingStrategy.TIMED
                    mj{'timingStrategy'} = 'TIMED';
                else
                    mj{'timingStrategy'} = 'IMMEDIATE';
                end
            end
            % Number of servers
            if mi <= length(node.numberOfServers) && node.numberOfServers(mi) > 1
                mj{'numServers'} = node.numberOfServers(mi);
            end
            % Firing priority
            % Omit only when it equals the builder default of 1: an explicit 0 is a
            % legal JMT firing priority and has to survive the round trip (BUG-90).
            if mi <= length(node.firingPriorities) && node.firingPriorities(mi) ~= 1
                mj{'firingPriority'} = node.firingPriorities(mi);
            end
            % Firing weight
            if mi <= length(node.firingWeights) && node.firingWeights(mi) ~= 1.0
                mj{'firingWeight'} = node.firingWeights(mi);
            end
            % Marking-dependent firing-rate multiplier g_m(marking), materialized
            % over the enabling (place,class) box lattice (cannot cross JSON as a
            % handle). Timed modes only; empty handle => omitted (unit multiplier).
            if ~isImmediateMode && mi <= length(node.firingRateDependence) ...
                    && ~isempty(node.firingRateDependence{mi})
                g = node.firingRateDependence{mi};
                ecMat = node.enablingConditions{mi};
                slots = {}; slotIdx = []; caps = [];
                for ni = 1:size(ecMat, 1)
                    for ci = 1:size(ecMat, 2)
                        if ecMat(ni, ci) > 0
                            sm = configureDictionary('string','cell');
                            sm{'node'} = allNodes{ni}.name;
                            sm{'class'} = classes{ci}.name;
                            slots{end+1} = sm; %#ok<AGROW>
                            slotIdx(end+1,:) = [ni, ci]; %#ok<AGROW>
                            pcap = allNodes{ni}.cap;
                            if ~isfinite(pcap), pcap = 10; end % open-place saturation cutoff
                            caps(end+1) = round(pcap); %#ok<AGROW>
                        end
                    end
                end
                if ~isempty(slots)
                    frm = configureDictionary('string','cell');
                    frm{'slots'} = slots;
                    frm{'cutoffs'} = num2cell(double(caps(:)'));
                    nnodes_all = length(allNodes);
                    frm{'scaling'} = firingdep_scaling_table(g, slotIdx, caps, nnodes_all, K);
                    mj{'firingRateDependence'} = frm;
                end
            end
            % Enabling conditions
            if mi <= length(node.enablingConditions)
                ecMat = node.enablingConditions{mi};
                ecList = {};
                for ni = 1:size(ecMat, 1)
                    for ci = 1:size(ecMat, 2)
                        if ecMat(ni, ci) > 0
                            ec = configureDictionary('string','cell');
                            ec{'node'} = allNodes{ni}.name;
                            ec{'class'} = classes{ci}.name;
                            ec{'count'} = ecMat(ni, ci);
                            ecList{end+1} = ec; %#ok<AGROW>
                        end
                    end
                end
                if ~isempty(ecList)
                    mj{'enablingConditions'} = ecList;
                end
            end
            % Inhibiting conditions
            if mi <= length(node.inhibitingConditions)
                icMat = node.inhibitingConditions{mi};
                icList = {};
                for ni = 1:size(icMat, 1)
                    for ci = 1:size(icMat, 2)
                        if isfinite(icMat(ni, ci))
                            ic = configureDictionary('string','cell');
                            ic{'node'} = allNodes{ni}.name;
                            ic{'class'} = classes{ci}.name;
                            ic{'count'} = icMat(ni, ci);
                            icList{end+1} = ic; %#ok<AGROW>
                        end
                    end
                end
                if ~isempty(icList)
                    mj{'inhibitingConditions'} = icList;
                end
            end
            % Firing outcomes
            if mi <= length(node.firingOutcomes)
                foMat = node.firingOutcomes{mi};
                foList = {};
                for ni = 1:size(foMat, 1)
                    for ci = 1:size(foMat, 2)
                        if foMat(ni, ci) ~= 0
                            fo = configureDictionary('string','cell');
                            fo{'node'} = allNodes{ni}.name;
                            fo{'class'} = classes{ci}.name;
                            fo{'count'} = foMat(ni, ci);
                            foList{end+1} = fo; %#ok<AGROW>
                        end
                    end
                end
                if ~isempty(foList)
                    mj{'firingOutcomes'} = foList;
                end
            end
            modesJson{end+1} = mj; %#ok<AGROW>
        end
        if ~isempty(modesJson)
            nj{'modes'} = modesJson;
        end
    end

    % The initial state of every stateful node, wrapped in a cell so a single-class
    % row still encodes as a JSON array. Not only a Place's token counts: a reader
    % calls the model initialized only when EVERY stateful node carries a state, so
    % naming one node and staying silent about the rest reads as uninitialized and
    % the reader rebuilds the default, discarding the state space and prior below.
    %
    % A NEGATIVE ENTRY IS NOT A STATE. It is the "ignore this station" flag of
    % the getProb* family -- solver_nc_margaggr sets lPr(ist) = NaN for it, and
    % statepr_aggr_large writes [-1,-1,-1,-1] on the two stations whose marginal
    % it is not asking about. The C++ reader spells a declared state as a one-row
    % state space and the analyzers START THE CHAIN THERE, so sending the flag
    % moved the answer for the station that WAS asked about: getProbAggr came
    % back 0.2027766, the probability of that station being EMPTY, against the
    % 0.005511140 every other row reports. The flag is a query argument; it stays
    % on this side, where the query is formed.
    if isa(node, 'StatefulNode') && nodeStateTravels && ~isempty(node.state) && all(double(node.state(:)) >= 0)
        nj{'initialState'} = num2cell(double(node.state(:)'));
    end

    % Emitted in the node loop before getStruct() overwrites statePrior with the default -- see _kb/09-ldes-and-cache.md
    %
    % A trivial [1] prior over one row is NOT emitted here, and does not need to
    % be: `initialState` above already carries that row, and the C++ reader
    % spells it as exactly this pair. What this block is for is a prior over
    % SEVERAL rows, which `initialState` cannot express.
    if isa(node, 'StatefulNode') && nodeStateTravels && ~isempty(node.statePrior)
        prior = double(node.statePrior(:));
        trivialPrior = numel(prior) == 1 && abs(prior(1) - 1.0) < 1e-12;
        space = double(full(node.space));
        % The same "ignore this station" flag as above is not a state space either.
        if trivialPrior || (~isempty(space) && any(space(:) < 0))
            % nothing to emit; initialState carries the single row, and the flag
            % is a query argument rather than a state
        elseif isempty(space) || size(space, 1) ~= numel(prior)
            line_warning(mfilename, sprintf(['Node %s carries a state prior over %d states but a ' ...
                'state space of %d rows; the prior is not saved.'], ...
                node.getName(), numel(prior), size(space, 1)));
        else
            spaceRows = cell(1, size(space, 1));
            for si = 1:size(space, 1)
                spaceRows{si} = num2cell(space(si, :));
            end
            nj{'stateSpace'} = spaceRows;
            % num2cell keeps a single-state prior from collapsing to a bare scalar (reader requires an array)
            nj{'statePrior'} = num2cell(prior');
        end
    end

    nodesJson{end+1} = nj; %#ok<AGROW>
end
result{'nodes'} = nodesJson;

% --- Classes ---
classesJson = {};
for r = 1:K
    jc = classes{r};
    cj = configureDictionary('string','cell');
    cj{'name'} = jc.name;
    if isa(jc, 'OpenSignal')
        cj{'type'} = 'Signal';
        cj{'openOrClosed'} = 'Open';
        cj{'signalType'} = SignalType.toText(jc.signalType);
        if ~isempty(jc.targetJobClass)
            cj{'targetClass'} = jc.targetJobClass.name;
        end
        if ~isempty(jc.removalDistribution)
            cj{'removalDistribution'} = dist2json(jc.removalDistribution);
        end
        if ~isempty(jc.removalPolicy) && jc.removalPolicy ~= RemovalPolicy.RANDOM
            cj{'removalPolicy'} = RemovalPolicy.toText(jc.removalPolicy);
        end
    elseif isa(jc, 'ClosedSignal')
        cj{'type'} = 'Signal';
        cj{'openOrClosed'} = 'Closed';
        cj{'signalType'} = SignalType.toText(jc.signalType);
        if ~isempty(jc.refstat) && isprop(jc.refstat, 'name')
            cj{'refNode'} = jc.refstat.name;
        end
        if ~isempty(jc.targetJobClass)
            cj{'targetClass'} = jc.targetJobClass.name;
        end
        if ~isempty(jc.removalDistribution)
            cj{'removalDistribution'} = dist2json(jc.removalDistribution);
        end
        if ~isempty(jc.removalPolicy) && jc.removalPolicy ~= RemovalPolicy.RANDOM
            cj{'removalPolicy'} = RemovalPolicy.toText(jc.removalPolicy);
        end
    elseif isa(jc, 'Signal')
        % Bare Signal (not OpenSignal/ClosedSignal): 'openOrClosed' omitted deliberately -- see _kb/09-ldes-and-cache.md
        cj{'type'} = 'Signal';
        cj{'signalType'} = SignalType.toText(jc.signalType);
        if ~isempty(jc.targetJobClass)
            cj{'targetClass'} = jc.targetJobClass.name;
        end
        if ~isempty(jc.removalDistribution)
            cj{'removalDistribution'} = dist2json(jc.removalDistribution);
        end
        if ~isempty(jc.removalPolicy) && jc.removalPolicy ~= RemovalPolicy.RANDOM
            cj{'removalPolicy'} = RemovalPolicy.toText(jc.removalPolicy);
        end
    elseif isa(jc, 'SelfLoopingClass')
        % Must be tested before ClosedClass (which it subclasses), else it reloads as an ordinary closed class
        cj{'type'} = 'SelfLooping';
        cj{'population'} = jc.population;
        if ~isempty(jc.refstat) && isprop(jc.refstat, 'name')
            cj{'refNode'} = jc.refstat.name;
        end
    elseif isa(jc, 'ClosedClass')
        cj{'type'} = 'Closed';
        cj{'population'} = jc.population;
        if ~isempty(jc.refstat) && isprop(jc.refstat, 'name')
            cj{'refNode'} = jc.refstat.name;
        end
    elseif isa(jc, 'OpenClass')
        cj{'type'} = 'Open';
        % refstat is normally re-derived as Source on load; carry it only when setReferenceStation overrode it
        if ~isempty(jc.refstat) && isprop(jc.refstat, 'name') && ~isa(jc.refstat, 'Source')
            cj{'refNode'} = jc.refstat.name;
        end
    else
        cj{'type'} = 'Open';
    end
    if jc.priority ~= 0
        cj{'priority'} = jc.priority;
    end
    if isprop(jc, 'deadline') && isfinite(jc.deadline)
        cj{'deadline'} = jc.deadline;
    end
    if jc.isReferenceClass()
        cj{'isReferenceClass'} = true;
    end
    % Reply signal binding (sn.syncreply). Without it a REPLY signal class is
    % inert after a round-trip: nothing unblocks the servers waiting on it.
    if isprop(jc, 'replySignalClass') && ~isempty(jc.replySignalClass)
        cj{'replySignalClass'} = jc.replySignalClass.name;
    end
    % Spawn-on-completion binding (sn.classspawn): the class injected at the
    % same station whenever a job of this class completes service.
    if isprop(jc, 'spawnClass') && ~isempty(jc.spawnClass)
        cj{'spawnClass'} = jc.spawnClass.name;
    end
    % Class-level (global) patience, distinct from the node-scoped 'patience'
    % emitted per Queue. A node-scoped entry overrides this one on load.
    if isprop(jc, 'patience') && ~isempty(jc.patience) && ~isa(jc.patience, 'Disabled')
        cj{'patience'} = dist2json(jc.patience);
        if ~isempty(jc.impatienceType)
            cj{'impatienceType'} = ImpatienceType.toText(jc.impatienceType);
        end
    end
    classesJson{end+1} = cj; %#ok<AGROW>
end
result{'classes'} = classesJson;

% --- Routing ---
routingMap = configureDictionary('string','cell');
try
    sn = model.getStruct();
    % Prefer rtorig (original P matrix before ClassSwitch expansion)
    if ~isempty(sn) && isfield(sn, 'rtorig') && iscell(sn.rtorig) && ~isempty(sn.rtorig) && ~isempty(sn.rtorig{1,1})
        P_orig = sn.rtorig;
        M_orig = size(P_orig{1,1}, 1);
        % Only non-identity explicit ClassSwitch nodes need same-class collapsing -- see _kb/04-networkstruct.md
        nodes = model.getNodes();
        explicit_cs = false(1, M_orig);
        for ii = 1:min(M_orig, length(nodes))
            if isa(nodes{ii}, 'ClassSwitch') && ~nodes{ii}.autoAdded ...
                    && ~isIdentityClassSwitch(nodes{ii}.server.csMatrix, K)
                explicit_cs(ii) = true;
            end
        end
        % P_same(s,ii,jj) = sum_r P_orig{r,s}(ii,jj): same-class-only entries to avoid double-switching -- see _kb/04-networkstruct.md
        cs_same = zeros(K, M_orig, M_orig);
        for ii = 1:M_orig
            if explicit_cs(ii)
                for s = 1:K
                    for jj = 1:M_orig
                        total = 0;
                        for r = 1:K
                            Prs = P_orig{r,s};
                            if issparse(Prs); Prs = full(Prs); end
                            total = total + Prs(ii, jj);
                        end
                        cs_same(s, ii, jj) = total;
                    end
                end
            end
        end
        for r = 1:K
            for s = 1:K
                fromTo = configureDictionary('string','cell');
                Prs = P_orig{r,s};
                if issparse(Prs)
                    Prs = full(Prs);
                end
                for ii = 1:M_orig
                    if explicit_cs(ii)
                        % For explicit CS, use same-class routing only
                        if r == s
                            for jj = 1:M_orig
                                val = cs_same(s, ii, jj);
                                if val > 1e-14
                                    ni = sn.nodenames{ii};
                                    njn = sn.nodenames{jj};
                                    if ~isKey(fromTo, ni)
                                        fromTo{ni} = configureDictionary('string','cell');
                                    end
                                    dest = fromTo{ni};
                                    dest{njn} = val;
                                    fromTo{ni} = dest;
                                end
                            end
                        end
                        % Skip cross-class entries from explicit CS
                        continue;
                    end
                    for jj = 1:M_orig
                        val = Prs(ii, jj);
                        if val > 1e-14
                            ni = sn.nodenames{ii};
                            njn = sn.nodenames{jj};
                            if ~isKey(fromTo, ni)
                                fromTo{ni} = configureDictionary('string','cell');
                            end
                            dest = fromTo{ni};
                            dest{njn} = val;
                            fromTo{ni} = dest;
                        end
                    end
                end
                if numEntries(fromTo) > 0
                    key = char(sprintf('%s,%s', classes{r}.name, classes{s}.name));
                    routingMap{key} = fromTo;
                end
            end
        end
        % rtorig IS PRE-LINK, so a (node,class) pair the user left to a routing
        % STRATEGY carries no row at all: link() resolved it against the model's
        % connection matrix, and exporting rtorig alone handed the reader a class
        % that never leaves its station. line-cli answered test_gallery_erlerl1
        % (Class1 routed only Source->Queue by hand) with an EMPTY AvgTable.
        % Fill only the pairs rtorig left silent, same-class, from the resolved
        % rtnodes; every explicit row is left exactly as the user wrote it.
        % Node indices are shared with rtnodes only when link() added no
        % ClassSwitch, and a Sink absorbs, so both are skipped.
        % A RETRIEVAL QUEUE'S SILENT ROW IS NOT AN UNSPECIFIED ONE. link()
        % broadcasts the read class's edges over the retrieval-system queue set
        % onto the per-item retrieval classes and then CONSUMES them (zeroes the
        % read class's own rows), so those rows are silent by construction, not
        % because the user left them to a strategy. Refilling them from the
        % RAND-resolved rtnodes writes a read-class TEMPLATE back into the
        % document, and a reader that re-runs the same broadcast overwrites the
        % explicit per-item rows with it: on retrieval_default the exported
        % Queue_1 -> {Cache 0.5, Queue_2 0.5} clobbered item 1's Queue_1 ->
        % Cache 1.0, and the native engine answered MVA ArvR 0.58216 / ResidT
        % 0.55089 where MATLAB answers 0.57599 / 0.54262.
        % Keyed on (read class, queue) and not on the queue alone, because a
        % retrieval queue may also carry an unrelated class whose row link()
        % really did resolve from a strategy.
        isRetrievalTemplate = false(K, M_orig);
        for ii = 1:min(M_orig, numel(nodes))
            if ~isa(nodes{ii}, 'Cache') || ~isprop(nodes{ii}, 'retrievalSystemQueueIndices')
                continue;
            end
            rcMap = nodes{ii}.retrievalSystemQueueIndices;
            rcKeys = keys(rcMap);
            for kk = 1:numel(rcKeys)
                r = double(rcKeys(kk)) + 1;   % read class index (key is index-1)
                if r < 1 || r > K
                    continue;
                end
                Q = rcMap{rcKeys(kk)};
                Q = Q(Q >= 1 & Q <= M_orig);
                isRetrievalTemplate(r, Q) = true;
            end
        end
        if isfield(sn, 'rtnodes') && ~isempty(sn.rtnodes) && M_orig == sn.nnodes
            for r = 1:K
                for ii = 1:M_orig
                    if explicit_cs(ii) || isRetrievalTemplate(r, ii) ...
                            || (ii <= numel(nodes) && isa(nodes{ii}, 'Sink'))
                        continue;
                    end
                    hasRow = false;
                    for s = 1:K
                        Prs = P_orig{r,s};
                        if issparse(Prs); Prs = full(Prs); end
                        if any(Prs(ii,:) > 1e-14)
                            hasRow = true;
                            break;
                        end
                    end
                    if hasRow
                        continue;
                    end
                    key = char(sprintf('%s,%s', classes{r}.name, classes{r}.name));
                    if isKey(routingMap, key)
                        fromTo = routingMap{key};
                    else
                        fromTo = configureDictionary('string','cell');
                    end
                    added = false;
                    for jj = 1:sn.nnodes
                        val = sn.rtnodes((ii-1)*K+r, (jj-1)*K+r);
                        if val > 1e-14
                            ni = sn.nodenames{ii};
                            njn = sn.nodenames{jj};
                            if ~isKey(fromTo, ni)
                                fromTo{ni} = configureDictionary('string','cell');
                            end
                            dest = fromTo{ni};
                            dest{njn} = val;
                            fromTo{ni} = dest;
                            added = true;
                        end
                    end
                    if added
                        routingMap{key} = fromTo;
                    end
                end
            end
        end
    elseif ~isempty(sn) && isfield(sn, 'rtnodes') && ~isempty(sn.rtnodes)
        % Fallback to rtnodes if rtorig not available
        rt = sn.rtnodes;
        N = sn.nnodes;
        % rtnodes folds CS switching into cross-class entries; mirror the rtorig branch (same-class only) -- see _kb/04-networkstruct.md
        nodes = model.getNodes();
        explicit_cs = false(1, N);
        for ii = 1:min(N, length(nodes))
            if isa(nodes{ii}, 'ClassSwitch') && ~nodes{ii}.autoAdded ...
                    && ~isIdentityClassSwitch(nodes{ii}.server.csMatrix, K)
                explicit_cs(ii) = true;
            end
        end
        % cs_same(s,ii,jj) = sum_r rt((ii,r),(jj,s)), same formula as the rtorig branch above
        cs_same = zeros(K, N, N);
        for ii = 1:N
            if explicit_cs(ii)
                for s = 1:K
                    for jj = 1:N
                        total = 0;
                        for r = 1:K
                            total = total + rt((ii-1)*K+r, (jj-1)*K+s);
                        end
                        cs_same(s, ii, jj) = total;
                    end
                end
            end
        end
        for r = 1:K
            for s = 1:K
                fromTo = configureDictionary('string','cell');
                for ii = 1:N
                    if explicit_cs(ii)
                        % For explicit CS, use same-class routing only
                        if r == s
                            for jj = 1:N
                                val = cs_same(s, ii, jj);
                                if val > 1e-14
                                    ni = sn.nodenames{ii};
                                    njn = sn.nodenames{jj};
                                    if ~isKey(fromTo, ni)
                                        fromTo{ni} = configureDictionary('string','cell');
                                    end
                                    dest = fromTo{ni};
                                    dest{njn} = val;
                                    fromTo{ni} = dest;
                                end
                            end
                        end
                        % Skip cross-class entries from explicit CS
                        continue;
                    end
                    for jj = 1:N
                        val = rt((ii-1)*K+r, (jj-1)*K+s);
                        if val > 1e-14
                            ni = sn.nodenames{ii};
                            njn = sn.nodenames{jj};
                            if ~isKey(fromTo, ni)
                                fromTo{ni} = configureDictionary('string','cell');
                            end
                            dest = fromTo{ni};
                            dest{njn} = val;
                            fromTo{ni} = dest;
                        end
                    end
                end
                if numEntries(fromTo) > 0
                    key = char(sprintf('%s,%s', classes{r}.name, classes{s}.name));
                    routingMap{key} = fromTo;
                end
            end
        end
    end
catch
    % If struct not available, routing stays empty
end

routing = configureDictionary('string','cell');
routing{'type'} = 'matrix';
routing{'matrix'} = routingMap;
result{'routing'} = routing;

% --- Global (Whittle) dependence ---
% phi(n) reads the FULL population matrix, so unlike the per-station
% classDependence/jointDependence blocks it is materialized over the lattice of
% the whole network state. Only the (station,class) SLOTS a class can actually
% occupy vary: a class disabled at a station, and any station that holds no jobs
% (a Source), contribute no coordinate, which is what keeps the lattice finite.
if ismethod(model,'getGlobalDependence') && ~isempty(model.getGlobalDependence())
    result{'globalDependence'} = gd_block(model);
end

% --- Routing Strategies ---
try
    sn2 = model.getStruct();
    if ~isempty(sn2) && isfield(sn2, 'routing') && ~isempty(sn2.routing)
        routingStrategies = configureDictionary('string','cell');
        stratNames = configureDictionary('int32','cell');
        stratNames{int32(RoutingStrategy.RAND)} = 'RAND';
        stratNames{int32(RoutingStrategy.RROBIN)} = 'RROBIN';
        stratNames{int32(RoutingStrategy.WRROBIN)} = 'WRROBIN';
        stratNames{int32(RoutingStrategy.JSQ)} = 'JSQ';
        stratNames{int32(RoutingStrategy.SQ)} = 'SQ';
        stratNames{int32(RoutingStrategy.FIRING)} = 'FIRING';
        % DISABLED is DERIVED, not declared: refreshRouting marks a (node,class)
        % pair the class never visits, every loader skips it on read, and the
        % pairs it lands on include the AUTO-ADDED class-switch nodes that the
        % node loop above deliberately does not emit. Writing it therefore
        % produced a strategies entry naming a node absent from `nodes` -- a
        % dangling reference in every re-entrant and cache model.
        % Same skip as the node loop: an implicit ClassSwitch is not in `nodes`, so naming it here would be a dangling reference.
        allNodes = model.getNodes();
        for i = 1:sn2.nnodes
            if i <= length(allNodes) && isa(allNodes{i}, 'ClassSwitch') ...
                    && isprop(allNodes{i}, 'autoAdded') && allNodes{i}.autoAdded
                continue;
            end
            nodeStrats = configureDictionary('string','cell');
            for r = 1:K
                % RAND is DECLARED not derived: dropping it made JMT write Empirical where model wants Random. PROB/DISABLED absent from stratNames, isKey skips both.
                routVal = int32(sn2.routing(i, r));
                if isKey(stratNames, routVal)
                    nodeStrats{classes{r}.name} = stratNames{routVal};
                end
            end
            if numEntries(nodeStrats) > 0
                routingStrategies{sn2.nodenames{i}} = nodeStrats;
            end
        end
        if numEntries(routingStrategies) > 0
            result{'routingStrategies'} = routingStrategies;
        end

        % Krzesinski state-dependent routing. Carried by NODE NAME so the block
        % is language independent and a node reordering on either side cannot
        % shift a center. Branch index 1 is the complement M-V and is written as
        % an empty list, keeping the paper's own numbering.
        % See _kb/16-state-dependent-routing.md
        for i = 1:sn2.nnodes
            for r = 1:K
                if sn2.routing(i,r) ~= RoutingStrategy.SDR
                    continue
                end
                decl = nodes{i}.output.outputStrategy{r}{3}{1};
                sdrBlock = configureDictionary('string','cell');
                sdrBlock{'entry'} = sn2.nodenames{i};
                sdrBlock{'departure'} = decl.departure.getName();
                sdrBlock{'class'} = sn2.classnames{r};
                brNames = cell(1, numel(decl.branch));
                brNames{1} = {};
                for b = 2:numel(decl.branch)
                    bn = cell(1, numel(decl.branch{b}));
                    for q = 1:numel(decl.branch{b})
                        bn{q} = decl.branch{b}{q}.getName();
                    end
                    brNames{b} = bn;
                end
                sdrBlock{'branches'} = brNames;
                sdrBlock{'level'} = decl.level(:)';
                sdrBlock{'C'} = decl.C(:)';
                sdrBlock{'d'} = decl.d;
                result{'stateDepRouting'} = sdrBlock;
                break
            end
            if isKey(result,'stateDepRouting')
                break
            end
        end

        % WRROBIN weights indexed by NODE index i, not station index -- see _kb/09-ldes-and-cache.md
        routingWeights = configureDictionary('string','cell');
        for i = 1:sn2.nnodes
            nodeObj2 = nodes{i};
            nodeClassWeights = configureDictionary('string','cell');
            for r = 1:K
                if int32(sn2.routing(i, r)) == int32(RoutingStrategy.WRROBIN)
                    os = nodeObj2.output.outputStrategy;
                    if size(os,2) >= r
                        osEntry = os{1, r};
                        % osEntry = {className, stratName, forwardLinks}
                        if length(osEntry) >= 3
                            fwdLinks = osEntry{3};
                            destWeights = configureDictionary('string','cell');
                            for fi = 1:length(fwdLinks)
                                link = fwdLinks{fi};
                                % link = {destNode, weight}
                                if iscell(link) && length(link) >= 2 && isa(link{1}, 'Node')
                                    destWeights{link{1}.name} = link{2};
                                end
                            end
                            if numEntries(destWeights) > 0
                                nodeClassWeights{classes{r}.name} = destWeights;
                            end
                        end
                    end
                end
            end
            if numEntries(nodeClassWeights) > 0
                routingWeights{sn2.nodenames{i}} = nodeClassWeights;
            end
        end
        if numEntries(routingWeights) > 0
            result{'routingWeights'} = routingWeights;
        end

        % SQ(d) sampling width, from the nodeparam slot saveRoutingStrategy reads; without it loader rebuilds default d=2 and SQ(d) silently becomes SQ(2).
        routingParams = configureDictionary('string','cell');
        for i = 1:sn2.nnodes
            nodeClassParams = configureDictionary('string','cell');
            for r = 1:K
                if int32(sn2.routing(i, r)) == int32(RoutingStrategy.SQ)
                    if length(sn2.nodeparam) >= i && ~isempty(sn2.nodeparam{i}) && ...
                            length(sn2.nodeparam{i}) >= r && ~isempty(sn2.nodeparam{i}{r}) && ...
                            isfield(sn2.nodeparam{i}{r}, 'd')
                        dEntry = configureDictionary('string','cell');
                        dEntry{'d'} = sn2.nodeparam{i}{r}.d;
                        nodeClassParams{classes{r}.name} = dEntry;
                    end
                end
            end
            if numEntries(nodeClassParams) > 0
                routingParams{sn2.nodenames{i}} = nodeClassParams;
            end
        end
        if numEntries(routingParams) > 0
            result{'routingParams'} = routingParams;
        end
    end
catch
end

% --- Setup / Delay-Off, Polling Type and Switchover Times ---
nodesCellTmp = result{'nodes'};
for i = 1:M
    nodeObj = nodes{i};
    if ~isa(nodeObj, 'Queue') || isa(nodeObj, 'Delay')
        continue;
    end
    njIdx = 0;
    for nj_idx = 1:length(nodesCellTmp)
        if strcmp(nodesCellTmp{nj_idx}('name'), nodeObj.name)
            njIdx = nj_idx;
            break;
        end
    end
    if njIdx == 0
        continue;
    end
    nj = nodesCellTmp{njIdx};

    % Setup and delay-off emitted as a pair: setDelayOff requires both on reload
    setupMap = configureDictionary('string','cell');
    delayOffMap = configureDictionary('string','cell');
    for r = 1:K
        if r <= length(nodeObj.setupTime) && r <= length(nodeObj.delayoffTime)
            suDist = nodeObj.setupTime{1,r};
            doffDist = nodeObj.delayoffTime{1,r};
            if ~isempty(suDist) && ~isempty(doffDist) && ...
                    ~isa(suDist, 'Disabled') && ~isa(doffDist, 'Disabled')
                setupMap{classes{r}.name} = dist2json(suDist);
                delayOffMap{classes{r}.name} = dist2json(doffDist);
            end
        end
    end
    if numEntries(setupMap) > 0
        nj{'setupTime'} = setupMap;
        nj{'delayOffTime'} = delayOffMap;
    end

    % Server breakdown/repair. The degraded down-server service is written per
    % class, which flattens the class-independent form setBreakdown also
    % accepts: both rebuild the same sn.downServiceRates row.
    if ~isempty(nodeObj.breakdownFailure) && ~isempty(nodeObj.breakdownRepair)
        bdMap = configureDictionary('string','cell');
        bdMap{'failure'} = dist2json(nodeObj.breakdownFailure);
        bdMap{'repair'} = dist2json(nodeObj.breakdownRepair);
        downMap = configureDictionary('string','cell');
        dsvcList = nodeObj.breakdownDownService;
        for r = 1:K
            dsvc = [];
            if length(dsvcList) == 1
                dsvc = dsvcList{1};
            elseif length(dsvcList) >= r
                dsvc = dsvcList{r};
            end
            if ~isempty(dsvc) && ~isa(dsvc, 'Disabled')
                downMap{classes{r}.name} = dist2json(dsvc);
            end
        end
        if numEntries(downMap) > 0
            bdMap{'downService'} = downMap;
        end
        nj{'breakdown'} = bdMap;
    end

    isPolling = SchedStrategy.toId(nodeObj.schedStrategy) == SchedStrategy.POLLING;

    % Polling type written by name: ids agree with Java but Python assigns via auto()
    if isPolling && ~isempty(nodeObj.pollingType)
        ptId = PollingType.toId(nodeObj.pollingType{1,1});
        nj{'pollingType'} = PollingType.toName(ptId);
        if ptId == PollingType.KLIMITED && ~isempty(nodeObj.pollingPar)
            nj{'pollingPar'} = nodeObj.pollingPar;
        end
    end

    % Under POLLING, switchover is indexed by departing class alone (no "to" field); otherwise a KxK (from,to) cell
    soTimes = {};
    if ~isempty(nodeObj.switchoverTime)
        [soRows, soCols] = size(nodeObj.switchoverTime);
        if isPolling
            for r = 1:min(K, soCols)
                dist = nodeObj.switchoverTime{1,r};
                if ~isempty(dist) && ~isa(dist, 'Disabled')
                    so = configureDictionary('string','cell');
                    so{'from'} = classes{r}.name;
                    so{'distribution'} = dist2json(dist);
                    soTimes{end+1} = so;
                end
            end
        else
            for r = 1:min(K, soRows)
                for s = 1:min(K, soCols)
                    dist = nodeObj.switchoverTime{r,s};
                    if ~isempty(dist) && ~isa(dist, 'Disabled')
                        so = configureDictionary('string','cell');
                        so{'from'} = classes{r}.name;
                        so{'to'} = classes{s}.name;
                        so{'distribution'} = dist2json(dist);
                        soTimes{end+1} = so;
                    end
                end
            end
        end
    end
    if ~isempty(soTimes)
        nj{'switchoverTimes'} = soTimes;
    end
    nodesCellTmp{njIdx} = nj;
end
result{'nodes'} = nodesCellTmp;

% --- Heterogeneous Server Types ---
try
    nodesCellTmp = result{'nodes'};
    for i = 1:M
        nodeObj = nodes{i};
        if isa(nodeObj, 'Queue') && nodeObj.isHeterogeneous()
            stArr = {};
            for ti = 1:length(nodeObj.serverTypes)
                st = nodeObj.serverTypes{ti};
                stj = configureDictionary('string','cell');
                stj{'name'} = st.name;
                stj{'count'} = st.numOfServers;
                % Compatible classes
                ccNames = {};
                for cci = 1:length(st.compatibleClasses)
                    ccNames{end+1} = st.compatibleClasses{cci}.name; %#ok<AGROW>
                end
                if ~isempty(ccNames)
                    stj{'compatibleClasses'} = ccNames;
                end
                % Per-class service distributions
                svcMap = configureDictionary('string','cell');
                for r = 1:K
                    jc = classes{r};
                    dist = nodeObj.getHeteroService(jc, st);
                    if ~isempty(dist) && ~isa(dist, 'Disabled')
                        svcMap{jc.name} = dist2json(dist);
                    end
                end
                if numEntries(svcMap) > 0
                    stj{'service'} = svcMap;
                end
                stArr{end+1} = stj; %#ok<AGROW>
            end
            if ~isempty(stArr)
                for nj_idx = 1:length(nodesCellTmp)
                    nj = nodesCellTmp{nj_idx};
                    if strcmp(nj{'name'}, nodeObj.name)
                        nj{'serverTypes'} = stArr;
                        % Scheduling policy
                        policy = nodeObj.getHeteroSchedPolicy();
                        if ~isempty(policy) && policy ~= HeteroSchedPolicy.ORDER
                            nj{'heteroSchedPolicy'} = HeteroSchedPolicy.toText(policy);
                        end
                        nodesCellTmp{nj_idx} = nj;
                        break;
                    end
                end
            end
        end
    end
    result{'nodes'} = nodesCellTmp;
catch
end

% --- Balking, Retrial, Patience, Orbit Impatience, Immediate Feedback ---
% Deliberately no try/catch: a bare catch here used to silently drop this entire block
nodesCellTmp = result{'nodes'};
nodeByName = configureDictionary('string','cell');
for i = 1:M
    nodeByName{nodes{i}.name} = nodes{i};
end
for nj_idx = 1:length(nodesCellTmp)
    nj = nodesCellTmp{nj_idx};
    nodeName = nj{'name'};
    nodeObj = nodeByName{nodeName};
    % Immediate feedback is a per-class node property on any Station.
    if isa(nodeObj, 'Queue')
        ifMap = configureDictionary('string','cell');
        for r = 1:K
            if nodeObj.hasImmediateFeedback(classes{r})
                ifMap{classes{r}.name} = true;
            end
        end
        if numEntries(ifMap) > 0
            nj{'immediateFeedback'} = ifMap;
        end
    end
    if isa(nodeObj, 'Queue')
        % Balking
        balkJson = configureDictionary('string','cell');
        for r = 1:K
            jc = classes{r};
            if nodeObj.hasBalking(jc)
                [strategy, thresholds] = nodeObj.getBalking(jc);
                bjc = configureDictionary('string','cell');
                switch strategy
                    case BalkingStrategy.QUEUE_LENGTH, bjc{'strategy'} = 'QUEUE_LENGTH';
                    case BalkingStrategy.EXPECTED_WAIT, bjc{'strategy'} = 'EXPECTED_WAIT';
                    case BalkingStrategy.COMBINED, bjc{'strategy'} = 'COMBINED';
                end
                thArr = {};
                for ti = 1:length(thresholds)
                    th = thresholds{ti};
                    tjson = configureDictionary('string','cell');
                    tjson{'minJobs'} = th{1};
                    if isinf(th{2})
                        tjson{'maxJobs'} = -1;
                    else
                        tjson{'maxJobs'} = th{2};
                    end
                    tjson{'probability'} = th{3};
                    thArr{end+1} = tjson;
                end
                bjc{'thresholds'} = thArr;
                balkJson{jc.name} = bjc;
            end
        end
        if numEntries(balkJson) > 0
            nj{'balking'} = balkJson;
        end
        % Retrial
        retrialJson = configureDictionary('string','cell');
        for r = 1:K
            jc = classes{r};
            if nodeObj.hasRetrial(jc)
                [delayDist, maxAttempts] = nodeObj.getRetrial(jc);
                rjc = configureDictionary('string','cell');
                rjc{'delay'} = dist2json(delayDist);
                rjc{'maxAttempts'} = maxAttempts;
                retrialJson{jc.name} = rjc;
            end
        end
        if numEntries(retrialJson) > 0
            nj{'retrial'} = retrialJson;
        end
        % Patience
        patienceJson = configureDictionary('string','cell');
        for r = 1:K
            jc = classes{r};
            patDist = nodeObj.getPatience(jc);
            if ~isempty(patDist) && ~isa(patDist, 'Disabled')
                pjc = configureDictionary('string','cell');
                pjc{'distribution'} = dist2json(patDist);
                impType = nodeObj.getImpatienceType(jc);
                if ~isempty(impType)
                    pjc{'impatienceType'} = ImpatienceType.toText(impType);
                end
                patienceJson{jc.name} = pjc;
            end
        end
        if numEntries(patienceJson) > 0
            nj{'patience'} = patienceJson;
        end
        % Orbit impatience (abandonment from the retrial orbit), distinct from
        % the queue patience above.
        orbitJson = configureDictionary('string','cell');
        for r = 1:K
            jc = classes{r};
            orbDist = nodeObj.getOrbitImpatience(jc);
            if ~isempty(orbDist) && ~isa(orbDist, 'Disabled')
                orbitJson{jc.name} = dist2json(orbDist);
            end
        end
        if numEntries(orbitJson) > 0
            nj{'orbitImpatience'} = orbitJson;
        end
        % Batch rejection probability (retrial queues), per class
        brpJson = configureDictionary('string','cell');
        for r = 1:K
            jc = classes{r};
            brp = nodeObj.getBatchRejectProbability(jc);
            if ~isempty(brp) && brp > 0
                brpJson{jc.name} = brp;
            end
        end
        if numEntries(brpJson) > 0
            nj{'batchRejectProb'} = brpJson;
        end
        % Job parallelism: servers seized at once by a job, per class
        parJson = configureDictionary('string','cell');
        for r = 1:K
            jc = classes{r};
            npar = nodeObj.getServerParallelism(jc);
            if npar > 1
                parJson{jc.name} = npar;
            end
        end
        if numEntries(parJson) > 0
            nj{'serverParallelism'} = parJson;
        end
    end
    nodesCellTmp{nj_idx} = nj;
end
result{'nodes'} = nodesCellTmp;

% --- Finite Capacity Regions ---
try
    regions = model.regions;
    if ~isempty(regions)
        fcrArray = {};
        for ri = 1:length(regions)
            reg = regions{ri};
            rj = configureDictionary('string','cell');
            rj{'name'} = reg.name;
            % Stations with per-class details
            stationsJson = {};
            for ni = 1:length(reg.nodes)
                sj = configureDictionary('string','cell');
                sj{'node'} = reg.nodes{ni}.name;
                % Per-class classCap
                if isprop(reg, 'classMaxJobs') && ~isempty(reg.classMaxJobs)
                    ccMap = configureDictionary('string','cell');
                    for r = 1:K
                        jc = classes{r};
                        if r <= length(reg.classMaxJobs) && isfinite(reg.classMaxJobs(r))
                            ccMap{jc.name} = reg.classMaxJobs(r);
                        end
                    end
                    if numEntries(ccMap) > 0
                        sj{'classCap'} = ccMap;
                    end
                end
                % Per-class classWeight
                if isprop(reg, 'classWeight') && ~isempty(reg.classWeight)
                    cwMap = configureDictionary('string','cell');
                    for r = 1:K
                        jc = classes{r};
                        if r <= length(reg.classWeight) && reg.classWeight(r) ~= 1
                            cwMap{jc.name} = reg.classWeight(r);
                        end
                    end
                    if numEntries(cwMap) > 0
                        sj{'classWeight'} = cwMap;
                    end
                end
                % Per-class classSize
                if isprop(reg, 'classSize') && ~isempty(reg.classSize)
                    csMap = configureDictionary('string','cell');
                    for r = 1:K
                        jc = classes{r};
                        if r <= length(reg.classSize) && reg.classSize(r) ~= 1
                            csMap{jc.name} = reg.classSize(r);
                        end
                    end
                    if numEntries(csMap) > 0
                        sj{'classSize'} = csMap;
                    end
                end
                stationsJson{end+1} = sj; %#ok<AGROW>
            end
            rj{'stations'} = stationsJson;
            if isprop(reg, 'globalMaxJobs') && isfinite(reg.globalMaxJobs)
                rj{'globalMaxJobs'} = reg.globalMaxJobs;
            end
            if isprop(reg, 'globalMaxMemory') && isfinite(reg.globalMaxMemory)
                rj{'globalMaxMemory'} = reg.globalMaxMemory;
            end
            % Per-class classMaxJobs at region level
            if isprop(reg, 'classMaxJobs') && ~isempty(reg.classMaxJobs)
                cmjMap = configureDictionary('string','cell');
                for r = 1:K
                    jc = classes{r};
                    if r <= length(reg.classMaxJobs) && isfinite(reg.classMaxJobs(r))
                        cmjMap{jc.name} = reg.classMaxJobs(r);
                    end
                end
                if numEntries(cmjMap) > 0
                    rj{'classMaxJobs'} = cmjMap;
                end
            end
            % Region classMaxMemory folds into the equivalent job cap on read -- see _kb/09-ldes-and-cache.md
            if isprop(reg, 'classMaxMemory') && ~isempty(reg.classMaxMemory)
                cmmMap = configureDictionary('string','cell');
                for r = 1:K
                    jc = classes{r};
                    if r <= length(reg.classMaxMemory) && isfinite(reg.classMaxMemory(r)) ...
                            && reg.classMaxMemory(r) >= 0
                        cmmMap{jc.name} = reg.classMaxMemory(r);
                    end
                end
                if numEntries(cmmMap) > 0
                    rj{'classMaxMemory'} = cmmMap;
                end
            end
            % Drop rule
            if isprop(reg, 'dropRule') && ~isempty(reg.dropRule)
                drMap = configureDictionary('string','cell');
                for r = 1:K
                    jc = classes{r};
                    if r <= length(reg.dropRule)
                        drStr = droprule_to_str(reg.dropRule(r));
                        if ~isempty(drStr)
                            drMap{jc.name} = drStr;
                        end
                    end
                end
                if numEntries(drMap) > 0
                    rj{'dropRule'} = drMap;
                end
            end
            % Linear constraints (A * x <= b) if present
            if ismethod(reg, 'hasLinearConstraints') && reg.hasLinearConstraints()
                [A, b] = reg.getLinearConstraints();
                if ~isempty(A) && ~isempty(b)
                    % Serialize A row-by-row as cell of arrays for JSON compat
                    Acell = cell(1, size(A,1));
                    for ri = 1:size(A,1)
                        Acell{ri} = A(ri,:);
                    end
                    rj{'constraintA'} = Acell;
                    rj{'constraintB'} = b(:)';
                end
            end
            fcrArray{end+1} = rj; %#ok<AGROW>
        end
        if ~isempty(fcrArray)
            result{'finiteCapacityRegions'} = fcrArray;
        end
    end
catch
end

% --- Rewards ---
rewardsJson = rewards2json(model);
if ~isempty(rewardsJson)
    result{'rewards'} = rewardsJson;
end
end


% =========================================================================
%  Reward serialization
% =========================================================================

function rewardsJson = rewards2json(model)
% Serialize the model's reward definitions in the declarative form
%   {name, type, node, class}
% Only rewards created through a Reward.* template carry the structural
% metadata needed to reproduce them. A reward defined from a bare function
% handle (or via Reward.custom) is not reproducible from JSON: warn and omit
% it rather than emit a reward that would be wrong on reload.
rewardsJson = {};
sn = model.sn;
if isempty(sn) || ~isfield(sn, 'reward') || isempty(sn.reward)
    return;
end
% Emitted in name order to match Python/JAR (JAR's HashMap has no insertion order)
rewardNames = cell(1, length(sn.reward));
for i = 1:length(sn.reward)
    rewardNames{i} = sn.reward{i}.name;
end
[~, order] = sort(rewardNames);
for oi = 1:length(order)
    rw = sn.reward{order(oi)};
    descriptor = [];
    if isfield(rw, 'descriptor')
        descriptor = rw.descriptor;
    end
    if isempty(descriptor) || ~isa(descriptor, 'RewardDescriptor')
        line_warning(mfilename, sprintf(['Reward "%s" is defined by a bare function handle and cannot be ' ...
            'serialized to JSON; it is omitted from the saved model. Use a Reward.* template ' ...
            '(Reward.queueLength/utilization/blocking) for a serializable reward.'], rw.name));
        continue;
    end
    if strcmp(descriptor.kind, 'Custom')
        line_warning(mfilename, sprintf(['Reward "%s" is a custom reward wrapping an arbitrary function and ' ...
            'cannot be serialized to JSON; it is omitted from the saved model.'], rw.name));
        continue;
    end
    rj = configureDictionary('string','cell');
    rj{'name'} = rw.name;
    rj{'type'} = descriptor.kind;
    if isempty(descriptor.node)
        line_warning(mfilename, sprintf(['Reward "%s" of type %s has no associated node and cannot be ' ...
            'serialized to JSON; it is omitted from the saved model.'], rw.name, descriptor.kind));
        continue;
    end
    rj{'node'} = descriptor.node.name;
    if ~isempty(descriptor.jobclass)
        rj{'class'} = descriptor.jobclass.name;
    end
    rewardsJson{end+1} = rj; %#ok<AGROW>
end
end


% =========================================================================
%  LayeredNetwork serialization
% =========================================================================

function result = layered2json(model)
result = configureDictionary('string','cell');
result{'type'} = 'LayeredNetwork';
result{'name'} = model.getName();

% --- Processors ---
procsJson = {};
hosts = model.hosts;
for i = 1:length(hosts)
    h = hosts{i};
    pj = configureDictionary('string','cell');
    pj{'name'} = h.name;
    mult = h.multiplicity;
    if ~isfinite(mult)
        pj{'multiplicity'} = inf_multiplicity();
    elseif mult > 1
        pj{'multiplicity'} = mult;
    end
    schedStr = h.scheduling;
    % INF is emitted like any other discipline, as the task path already does:
    % suppressing it left a reader to default, and PS is the common default
    if ~isempty(schedStr)
        pj{'scheduling'} = upper(schedStr);
    end
    q = h.quantum;
    if q > 0 && q ~= 0.001
        pj{'quantum'} = q;
    end
    sf = h.speedFactor;
    if sf ~= 1.0
        pj{'speedFactor'} = sf;
    end
    repl = h.replication;
    if repl > 1
        pj{'replication'} = repl;
    end
    % Admission constraints: columns of a host constraint are its tasks
    colNames = {};
    for ci = 1:length(h.tasks)
        colNames{end+1} = h.tasks(ci).name; %#ok<AGROW>
    end
    rowsJson = lincon2json(h, colNames);
    if ~isempty(rowsJson)
        pj{'admissionConstraints'} = rowsJson;
    end
    procsJson{end+1} = pj; %#ok<AGROW>
end
result{'hosts'} = procsJson;

% --- Tasks ---
tasksJson = {};
tasksList = model.tasks;
for i = 1:length(tasksList)
    t = tasksList{i};
    tj = configureDictionary('string','cell');
    tj{'name'} = t.name;
    if ~isempty(t.parent)
        tj{'host'} = t.parent.name;
    end
    mult = t.multiplicity;
    if ~isfinite(mult)
        tj{'multiplicity'} = inf_multiplicity();
    elseif mult > 1
        tj{'multiplicity'} = mult;
    end
    schedStr = t.scheduling;
    if ~isempty(schedStr)
        tj{'scheduling'} = upper(schedStr);
    end
    % Think time
    ttMean = t.thinkTimeMean;
    if ~isempty(ttMean) && ttMean > GlobalConstants.FineTol
        if ~isempty(t.thinkTime) && isa(t.thinkTime, 'Distribution')
            tj{'thinkTime'} = dist2json(t.thinkTime);
        else
            params = configureDictionary('string','cell');
            params{'lambda'} = 1.0 / ttMean;
            dj = configureDictionary('string','cell');
            dj{'type'} = 'Exp';
            dj{'params'} = params;
            tj{'thinkTime'} = dj;
        end
    end
    % Fan in
    if ~isempty(t.fanInSource) && ischar(t.fanInSource) && ~isempty(t.fanInSource)
        fi = configureDictionary('string','cell');
        fi{t.fanInSource} = t.fanInValue;
        tj{'fanIn'} = fi;
    end
    % Fan out
    if ~isempty(t.fanOutDest)
        fo = configureDictionary('string','cell');
        for fi_idx = 1:length(t.fanOutDest)
            fo{t.fanOutDest{fi_idx}} = t.fanOutValue(fi_idx);
        end
        tj{'fanOut'} = fo;
    end
    repl = t.replication;
    if repl > 1
        tj{'replication'} = repl;
    end
    % Admission constraints: columns of a task constraint are its entries
    colNames = {};
    for ci = 1:length(t.entries)
        colNames{end+1} = t.entries(ci).name; %#ok<AGROW>
    end
    rowsJson = lincon2json(t, colNames);
    if ~isempty(rowsJson)
        tj{'admissionConstraints'} = rowsJson;
    end
    % SetupTask detection
    if isa(t, 'SetupTask')
        tj{'taskType'} = 'SetupTask';
    end
    % Setup time / delay-off time (on any Task)
    if ~isempty(t.setupTime) && isa(t.setupTime, 'Distribution')
        stMean = t.setupTimeMean;
        if stMean > GlobalConstants.FineTol
            tj{'setupTime'} = dist2json(t.setupTime);
        end
    end
    if ~isempty(t.delayOffTime) && isa(t.delayOffTime, 'Distribution')
        dotMean = t.delayOffTimeMean;
        if dotMean > GlobalConstants.FineTol
            tj{'delayOffTime'} = dist2json(t.delayOffTime);
        end
    end
    % CacheTask detection
    if isa(t, 'CacheTask')
        tj{'taskType'} = 'CacheTask';
        tj{'totalItems'} = t.items;
        tj{'cacheCapacity'} = t.itemLevelCap;
        rs = t.replacestrategy;
        rsNameMap = dictionary([ReplacementStrategy.RR, ReplacementStrategy.FIFO, ...
            ReplacementStrategy.SFIFO, ReplacementStrategy.LRU], ...
            {'RR', 'FIFO', 'SFIFO', 'LRU'});
        if isKey(rsNameMap, rs)
            tj{'replacementStrategy'} = rsNameMap{rs};
        else
            tj{'replacementStrategy'} = 'FIFO';
        end
    end
    tasksJson{end+1} = tj; %#ok<AGROW>
end
result{'tasks'} = tasksJson;

% --- Entries ---
entriesJson = {};
entriesList = model.entries;
for i = 1:length(entriesList)
    e = entriesList{i};
    ej = configureDictionary('string','cell');
    ej{'name'} = e.name;
    if ~isempty(e.parent)
        ej{'task'} = e.parent.name;
    end
    % Entry arrival distribution
    if ~isempty(e.arrival) && isa(e.arrival, 'Distribution')
        ej{'arrival'} = dist2json(e.arrival);
    end
    % ItemEntry detection
    if isa(e, 'ItemEntry')
        ej{'entryType'} = 'ItemEntry';
        ej{'totalItems'} = e.cardinality;
        if ~isempty(e.popularity)
            if isa(e.popularity, 'Distribution')
                ej{'accessProb'} = dist2json(e.popularity);
            end
        end
    end
    entriesJson{end+1} = ej; %#ok<AGROW>
end
result{'entries'} = entriesJson;

% --- Build reply map: activityName -> entryName ---
replyMap = configureDictionary('string','cell');
for i = 1:length(entriesList)
    e = entriesList{i};
    if ~isempty(e.replyActivity)
        for j = 1:length(e.replyActivity)
            replyMap{e.replyActivity{j}} = e.name;
        end
    end
end

% --- Activities ---
actsJson = {};
actsList = model.activities;
for i = 1:length(actsList)
    a = actsList{i};
    aj = configureDictionary('string','cell');
    aj{'name'} = a.name;
    if ~isempty(a.parent)
        if isa(a.parent, 'Task') || isa(a.parent, 'Entry')
            aj{'task'} = a.parent.name;
        elseif ischar(a.parent) || isstring(a.parent)
            aj{'task'} = char(a.parent);
        elseif ischar(a.parentName) && ~isempty(a.parentName)
            aj{'task'} = a.parentName;
        end
    elseif ~isempty(a.parentName) && ischar(a.parentName)
        aj{'task'} = a.parentName;
    end
    % Host demand
    if ~isempty(a.hostDemand) && isa(a.hostDemand, 'Distribution')
        if ~isa(a.hostDemand, 'Immediate')
            aj{'hostDemand'} = dist2json(a.hostDemand);
        end
    elseif ~isempty(a.hostDemandMean) && a.hostDemandMean > GlobalConstants.FineTol
        params = configureDictionary('string','cell');
        params{'lambda'} = 1.0 / a.hostDemandMean;
        dj = configureDictionary('string','cell');
        dj{'type'} = 'Exp';
        dj{'params'} = params;
        aj{'hostDemand'} = dj;
    end
    % Bound to entry
    if ~isempty(a.boundToEntry)
        aj{'boundToEntry'} = a.boundToEntry;
    end
    % Replies to entry
    if isKey(replyMap, a.name)
        aj{'repliesTo'} = replyMap{a.name};
    end
    % Synch calls
    if ~isempty(a.syncCallDests)
        synchCalls = {};
        for j = 1:length(a.syncCallDests)
            sc = configureDictionary('string','cell');
            sc{'dest'} = a.syncCallDests{j};
            if j <= length(a.syncCallMeans) && a.syncCallMeans(j) ~= 1.0
                sc{'mean'} = a.syncCallMeans(j);
            end
            synchCalls{end+1} = sc; %#ok<AGROW>
        end
        aj{'synchCalls'} = synchCalls;
    end
    % Asynch calls
    if ~isempty(a.asyncCallDests)
        asynchCalls = {};
        for j = 1:length(a.asyncCallDests)
            ac = configureDictionary('string','cell');
            ac{'dest'} = a.asyncCallDests{j};
            if j <= length(a.asyncCallMeans) && a.asyncCallMeans(j) ~= 1.0
                ac{'mean'} = a.asyncCallMeans(j);
            end
            asynchCalls{end+1} = ac; %#ok<AGROW>
        end
        aj{'asynchCalls'} = asynchCalls;
    end
    actsJson{end+1} = aj; %#ok<AGROW>
end
result{'activities'} = actsJson;

% --- Precedences ---
precsJson = {};
for i = 1:length(tasksList)
    t = tasksList{i};
    precs = t.precedences;
    if isempty(precs), continue; end
    for j = 1:length(precs)
        p = precs(j);
        pj = configureDictionary('string','cell');
        pj{'task'} = t.name;

        preType  = p.preType;
        postType = p.postType;

        % Determine JSON precedence type and collect activity names
        if preType == ActivityPrecedenceType.PRE_SEQ && postType == ActivityPrecedenceType.POST_SEQ
            pj{'type'} = 'Serial';
            pj{'activities'} = [p.preActs, p.postActs];
        elseif preType == ActivityPrecedenceType.PRE_SEQ && postType == ActivityPrecedenceType.POST_AND
            pj{'type'} = 'AndFork';
            pj{'activities'} = [p.preActs, p.postActs];
        elseif preType == ActivityPrecedenceType.PRE_AND && postType == ActivityPrecedenceType.POST_SEQ
            pj{'type'} = 'AndJoin';
            pj{'activities'} = [p.preActs, p.postActs];
        elseif preType == ActivityPrecedenceType.PRE_SEQ && postType == ActivityPrecedenceType.POST_OR
            pj{'type'} = 'OrFork';
            pj{'activities'} = [p.preActs, p.postActs];
            if ~isempty(p.postParams)
                pj{'probabilities'} = p.postParams(:)';
            end
        elseif preType == ActivityPrecedenceType.PRE_OR && postType == ActivityPrecedenceType.POST_SEQ
            pj{'type'} = 'OrJoin';
            pj{'activities'} = [p.preActs, p.postActs];
        elseif postType == ActivityPrecedenceType.POST_LOOP
            pj{'type'} = 'Loop';
            % For Loop, preActs is the trigger, postActs is the loop body
            pj{'activities'} = p.postActs;
            if ~isempty(p.preActs)
                pj{'preActivity'} = p.preActs{1};
            end
            if ~isempty(p.postParams)
                pj{'loopCount'} = p.postParams(1);
            end
        elseif preType == ActivityPrecedenceType.PRE_SEQ && postType == ActivityPrecedenceType.POST_CACHE
            pj{'type'} = 'CacheAccess';
            pj{'activities'} = [p.preActs, p.postActs];
        else
            continue;
        end
        precsJson{end+1} = pj; %#ok<AGROW>
    end
end
if ~isempty(precsJson)
    result{'precedences'} = precsJson;
end
end


% =========================================================================
%  Workflow serialization
% =========================================================================

function result = workflow2json(model)
% Convert a Workflow to a dictionary for JSON output.
result = configureDictionary('string','cell');
result{'type'} = 'Workflow';
result{'name'} = model.getName();

% --- Activities ---
actsJson = {};
acts = model.activities;
for i = 1:length(acts)
    act = acts{i};
    aj = configureDictionary('string','cell');
    aj{'name'} = act.name;
    if ~isempty(act.hostDemand) && isa(act.hostDemand, 'Distribution')
        dj = dist2json(act.hostDemand);
        if ~isempty(dj)
            aj{'hostDemand'} = dj;
        end
    end
    actsJson{end+1} = aj; %#ok<AGROW>
end
result{'activities'} = actsJson;

% --- Precedences ---
precsJson = {};
precs = model.precedences;
for i = 1:length(precs)
    p = precs(i);
    pj = configureDictionary('string','cell');

    % preActs
    preActsJson = {};
    for a = 1:length(p.preActs)
        preActsJson{end+1} = p.preActs{a}; %#ok<AGROW>
    end
    pj{'preActs'} = preActsJson;

    % postActs
    postActsJson = {};
    for a = 1:length(p.postActs)
        postActsJson{end+1} = p.postActs{a}; %#ok<AGROW>
    end
    pj{'postActs'} = postActsJson;

    % preType / postType - convert numeric IDs to JAR-compatible strings
    pj{'preType'} = prectype_to_str(p.preType);
    pj{'postType'} = prectype_to_str(p.postType);

    % preParams
    if ~isempty(p.preParams)
        pj{'preParams'} = p.preParams(:)';
    end

    % postParams
    if ~isempty(p.postParams)
        pj{'postParams'} = p.postParams(:)';
    end

    precsJson{end+1} = pj; %#ok<AGROW>
end
result{'precedences'} = precsJson;
end


% =========================================================================
%  Environment serialization
% =========================================================================

function result = environment2json(model)
% Convert an Environment to a dictionary for JSON output.
result = configureDictionary('string','cell');
result{'type'} = 'Environment';
result{'name'} = model.getName();

E = height(model.envGraph.Nodes);
result{'numStages'} = E;

% --- Stages ---
stagesJson = {};
for e = 1:E
    sj = configureDictionary('string','cell');
    sj{'name'} = model.envGraph.Nodes.Name{e};
    % The stage TYPE, which the reader already looks for ("type") and which no
    % writer emitted: an Environment round-tripped through JSON came back with
    % every stage type blanked, so getStageTable() and any consumer keying off
    % UP/DOWN read a different environment from the one that was saved.
    if ~isempty(model.envGraph.Nodes.Type{e})
        sj{'type'} = char(model.envGraph.Nodes.Type{e});
    end
    % Serialize the stage's Network model
    if e <= length(model.ensemble) && ~isempty(model.ensemble{e})
        sj{'model'} = network2json(model.ensemble{e});
    end
    stagesJson{end+1} = sj; %#ok<AGROW>
end
result{'stages'} = stagesJson;

% --- Transitions ---
transJson = {};
for e = 1:E
    for h = 1:E
        if ~isempty(model.env) && e <= size(model.env, 1) && h <= size(model.env, 2) ...
                && ~isempty(model.env{e,h}) && ~isa(model.env{e,h}, 'Disabled')
            tj = configureDictionary('string','cell');
            tj{'from'} = e - 1;  % Convert to 0-indexed for JAR compatibility
            tj{'to'} = h - 1;    % Convert to 0-indexed for JAR compatibility
            dj = dist2json(model.env{e,h});
            if ~isempty(dj)
                tj{'distribution'} = dj;
                transJson{end+1} = tj; %#ok<AGROW>
            end
        end
    end
end
result{'transitions'} = transJson;

% --- Node failures ---
% Declarative record of addNodeBreakdown/addNodeRepair; carries the queue-length reset policies (function handles), otherwise unrecoverable
nfJson = {};
for i = 1:length(model.nodeFailures)
    nf = model.nodeFailures{i};
    nj = configureDictionary('string','cell');
    nj{'node'} = nf.node;
    bj = dist2json(nf.breakdown);
    if isempty(bj)
        line_warning(mfilename, sprintf(['Node failure on "%s" has a breakdown distribution that cannot be ' ...
            'serialized; the nodeFailures entry is omitted.'], nf.node));
        continue;
    end
    nj{'breakdownRate'} = bj;
    if ~isempty(nf.repair)
        rj = dist2json(nf.repair);
        if isempty(rj)
            line_warning(mfilename, sprintf(['Node failure on "%s" has a repair distribution that cannot be ' ...
                'serialized; the nodeFailures entry is omitted.'], nf.node));
            continue;
        end
        nj{'repairRate'} = rj;
    end
    dj = dist2json(nf.downService);
    if isempty(dj)
        line_warning(mfilename, sprintf(['Node failure on "%s" has a down-service distribution that cannot ' ...
            'be serialized; the nodeFailures entry is omitted.'], nf.node));
        continue;
    end
    nj{'downService'} = dj;
    if strcmp(nf.breakdownResetPolicy, 'custom')
        line_warning(mfilename, sprintf(['Node failure on "%s" uses a custom breakdown reset function, which ' ...
            'cannot be serialized to JSON; the saved model falls back to the ''keep'' policy on reload.'], nf.node));
    else
        nj{'breakdownResetPolicy'} = nf.breakdownResetPolicy;
    end
    if ~isempty(nf.repairResetPolicy)
        if strcmp(nf.repairResetPolicy, 'custom')
            line_warning(mfilename, sprintf(['Node failure on "%s" uses a custom repair reset function, which ' ...
                'cannot be serialized to JSON; the saved model falls back to the ''keep'' policy on reload.'], nf.node));
        else
            nj{'repairResetPolicy'} = nf.repairResetPolicy;
        end
    end
    nfJson{end+1} = nj; %#ok<AGROW>
end
if ~isempty(nfJson)
    result{'nodeFailures'} = nfJson;
end
end


% =========================================================================
%  Admission constraint serialization
% =========================================================================

function rowsJson = lincon2json(elem, colNames)
% ROWSJSON = LINCON2JSON(ELEM, COLNAMES) admission constraint rows of a Task or
% Host on the wire.
%
% Both declaration forms are normalised to the named form, so the wire is
% order-independent: a positional setConstraint(A,b) matrix is resolved against
% COLNAMES (the element's entries, or its tasks) at write time. COLNAMES must be
% in the same declaration order the positional columns assume.
rowsJson = {};
if isempty(elem.linConA) && isempty(elem.linConRows)
    return
end
A = elem.linConA;
for r = 1:size(A,1)
    nz = find(A(r,:) ~= 0);
    if isempty(nz)
        continue
    end
    if any(nz > length(colNames))
        line_error(mfilename,'Admission constraint on %s references column %d but the element has only %d operands.', elem.name, max(nz), length(colNames));
    end
    rj = configureDictionary('string','cell');
    rj{'operands'} = colNames(nz);
    rj{'coeffs'} = A(r,nz);
    rj{'cap'} = elem.linConB(r);
    rowsJson{end+1} = rj; %#ok<AGROW>
end
for r = 1:length(elem.linConRows)
    namedRow = elem.linConRows{r};
    rj = configureDictionary('string','cell');
    rj{'operands'} = namedRow.names;
    rj{'coeffs'} = namedRow.coeffs;
    rj{'cap'} = namedRow.cap;
    rowsJson{end+1} = rj; %#ok<AGROW>
end
end

% =========================================================================
%  Multiplicity serialization
% =========================================================================

function v = inf_multiplicity()
% Wire sentinel for infinite host/task multiplicity: Java's Integer.MAX_VALUE,
% which the JAR uses as its infinite-multiplicity marker. Written as a literal
% rather than taken from GlobalConstants.MaxInt, which is settable at runtime:
% a sentinel that varies per session would not survive a round-trip between two
% differently configured readers.
v = 2147483647;
end

% =========================================================================
%  Distribution serialization
% =========================================================================

function d = dist2json(dist)
% Convert a Distribution to a dictionary for JSON output.
if isempty(dist)
    d = [];
    return;
end
d = configureDictionary('string','cell');
cn = builtin('class', dist);
switch cn
    case 'Disabled'
        d{'type'} = 'Disabled';
    case 'Immediate'
        d{'type'} = 'Immediate';
    case 'Exp'
        d{'type'} = 'Exp';
        params = configureDictionary('string','cell');
        params{'lambda'} = dist.getParam(1).paramValue;
        d{'params'} = params;
    case 'Det'
        d{'type'} = 'Det';
        params = configureDictionary('string','cell');
        params{'value'} = dist.getParam(1).paramValue;
        d{'params'} = params;
    case 'Erlang'
        d{'type'} = 'Erlang';
        params = configureDictionary('string','cell');
        params{'lambda'} = dist.getParam(1).paramValue;
        params{'k'} = dist.getParam(2).paramValue;
        d{'params'} = params;
    case 'HyperExp'
        % param2/param3 both store the FULL n-phase rate vector -- see _kb/09-ldes-and-cache.md
        d{'type'} = 'HyperExp';
        params = configureDictionary('string','cell');
        p = dist.getParam(1).paramValue;
        l1 = dist.getParam(2).paramValue;
        l2 = dist.getParam(3).paramValue;
        if isscalar(p)
            params{'p'} = [p, 1-p];
            params{'lambda'} = [l1, l2];
        else
            params{'p'} = p(:)';
            if isscalar(l1) && isscalar(l2)
                params{'lambda'} = [l1, l2];
            else
                % n-phase: param2 already holds all n rates (param3 duplicates it)
                params{'lambda'} = l1(:)';
            end
        end
        d{'params'} = params;
    case 'Gamma'
        d{'type'} = 'Gamma';
        params = configureDictionary('string','cell');
        params{'alpha'} = dist.getParam(1).paramValue;
        params{'beta'} = dist.getParam(2).paramValue;
        d{'params'} = params;
    case 'Lognormal'
        d{'type'} = 'Lognormal';
        params = configureDictionary('string','cell');
        params{'mu'} = dist.getParam(1).paramValue;
        params{'sigma'} = dist.getParam(2).paramValue;
        d{'params'} = params;
    case 'Uniform'
        d{'type'} = 'Uniform';
        params = configureDictionary('string','cell');
        params{'a'} = dist.getParam(1).paramValue;
        params{'b'} = dist.getParam(2).paramValue;
        d{'params'} = params;
    case 'Zipf'
        d{'type'} = 'Zipf';
        params = configureDictionary('string','cell');
        params{'s'} = dist.getParam(3).paramValue;
        params{'n'} = dist.getParam(4).paramValue;
        d{'params'} = params;
    case 'Pareto'
        d{'type'} = 'Pareto';
        params = configureDictionary('string','cell');
        params{'alpha'} = dist.getParam(1).paramValue;
        params{'scale'} = dist.getParam(2).paramValue;
        d{'params'} = params;
    case 'Weibull'
        d{'type'} = 'Weibull';
        params = configureDictionary('string','cell');
        params{'alpha'} = dist.getParam(1).paramValue;
        params{'beta'} = dist.getParam(2).paramValue;
        d{'params'} = params;
    case 'Normal'
        d{'type'} = 'Normal';
        params = configureDictionary('string','cell');
        params{'mu'} = dist.getParam(1).paramValue;
        params{'sigma'} = dist.getParam(2).paramValue;
        d{'params'} = params;
    case 'Geometric'
        d{'type'} = 'Geometric';
        params = configureDictionary('string','cell');
        params{'p'} = dist.getParam(1).paramValue;
        d{'params'} = params;
    case 'Binomial'
        d{'type'} = 'Binomial';
        params = configureDictionary('string','cell');
        params{'n'} = dist.getParam(1).paramValue;
        params{'p'} = dist.getParam(2).paramValue;
        d{'params'} = params;
    case 'Poisson'
        d{'type'} = 'Poisson';
        params = configureDictionary('string','cell');
        params{'lambda'} = dist.getParam(1).paramValue;
        d{'params'} = params;
    case 'Bernoulli'
        d{'type'} = 'Bernoulli';
        params = configureDictionary('string','cell');
        params{'p'} = dist.getParam(1).paramValue;
        d{'params'} = params;
    case 'DiscreteUniform'
        d{'type'} = 'DiscreteUniform';
        params = configureDictionary('string','cell');
        params{'min'} = dist.getParam(1).paramValue;
        params{'max'} = dist.getParam(2).paramValue;
        d{'params'} = params;
    case {'Coxian', 'Cox2'}
        d{'type'} = 'Coxian';
        params = configureDictionary('string','cell');
        params{'mu'} = dist.getMu()';
        params{'phi'} = dist.getPhi()';
        d{'params'} = params;
    case {'PH', 'APH'}
        % Keep the concrete class: an APH written back as a generic PH is a
        % lossy downgrade, since solver feature sets admit APH but not PH.
        d{'type'} = cn;
        ph = configureDictionary('string','cell');
        alpha = dist.getInitProb();
        T = dist.getSubgenerator();
        if isvector(alpha)
            ph{'alpha'} = num2cell(double(full(alpha(:)')));
        else
            ph{'alpha'} = mat2json(alpha);
        end
        ph{'T'} = mat2json(T);
        d{'ph'} = ph;
    case 'MAP'
        d{'type'} = 'MAP';
        mapSpec = configureDictionary('string','cell');
        mapSpec{'D0'} = mat2json(dist.getParam(1).paramValue);
        mapSpec{'D1'} = mat2json(dist.getParam(2).paramValue);
        d{'map'} = mapSpec;
    case 'DMAP'
        % Discrete-time MAP. Distinct from MAP on the wire: D0+D1 is stochastic,
        % not an infinitesimal generator, so a reader must not rebuild it as MAP.
        d{'type'} = 'DMAP';
        params = configureDictionary('string','cell');
        params{'D0'} = mat2json(dist.getParam(1).paramValue);
        params{'D1'} = mat2json(dist.getParam(2).paramValue);
        d{'params'} = params;
    case {'ME', 'CME'}
        % CME goes on the wire as its (alpha, A) ME representation; the subclass tag is not preserved
        d{'type'} = 'ME';
        params = configureDictionary('string','cell');
        alphaME = dist.getParam(1).paramValue;
        params{'alpha'} = num2cell(double(full(alphaME(:)')));
        params{'A'} = mat2json(dist.getParam(2).paramValue);
        d{'params'} = params;
    case 'RAP'
        d{'type'} = 'RAP';
        params = configureDictionary('string','cell');
        params{'H0'} = mat2json(dist.getParam(1).paramValue);
        params{'H1'} = mat2json(dist.getParam(2).paramValue);
        d{'params'} = params;
    case 'BMAP'
        % Must never be emitted through the MMAP branch (mark index is a batch size, not a class) -- see _kb/09-ldes-and-cache.md
        d{'type'} = 'BMAP';
        params = configureDictionary('string','cell');
        Kb = dist.getNumberOfTypes;
        dArr = cell(1, Kb + 1);
        dArr{1} = mat2json(dist.getParam(1).paramValue);      % D0
        for kb = 1:Kb
            dArr{1+kb} = mat2json(dist.getParam(2+kb).paramValue);   % Dk, batch size k
        end
        params{'D'} = dArr;
        d{'params'} = params;
    case 'MMDP2'
        d{'type'} = 'MMDP2';
        params = configureDictionary('string','cell');
        params{'r0'} = dist.getParam(1).paramValue;
        params{'r1'} = dist.getParam(2).paramValue;
        params{'sigma0'} = dist.getParam(3).paramValue;
        params{'sigma1'} = dist.getParam(4).paramValue;
        d{'params'} = params;
    case 'MarkedMMPP'
        % Extends MarkovModulated, not MarkedMAP, so needs its own branch -- see _kb/09-ldes-and-cache.md
        d{'type'} = 'MarkedMMPP';
        params = configureDictionary('string','cell');
        Km = dist.getNumberOfTypes;
        dArr = cell(1, Km + 1);
        dArr{1} = mat2json(dist.getParam(1).paramValue);      % D0
        for km = 1:Km
            dArr{1+km} = mat2json(dist.getParam(2+km).paramValue);   % D1k
        end
        params{'D'} = dArr;
        params{'K'} = Km;
        d{'params'} = params;
    case 'EmpiricalCDF'
        % data is [F, x] rows (cdf value, support point), as assembled by the
        % two-argument ctor; emit the two columns separately.
        d{'type'} = 'EmpiricalCDF';
        params = configureDictionary('string','cell');
        ecdf = dist.data;
        params{'F'} = ecdf(:, 1)';
        params{'x'} = ecdf(:, 2)';
        d{'params'} = params;
    case 'Expolynomial'
        % Nested object, mirroring the Python writer (the density is an expolynomial
        % expression string, not a numeric parameter).
        d{'type'} = 'Expolynomial';
        ep = configureDictionary('string','cell');
        ep{'density'} = dist.getParam(1).paramValue;
        ep{'eft'} = dist.getParam(2).paramValue;
        lft = dist.getParam(3).paramValue;
        if isfinite(lft)
            ep{'lft'} = lft;
        else
            ep{'lft'} = 'Inf';
        end
        d{'expolynomial'} = ep;
    case 'MarkedMAP'
        % Marked MAP: {D0, per-mark D1k}; the aggregate D1 is rebuilt on load
        d{'type'} = 'MMAP';
        mmapSpec = configureDictionary('string','cell');
        mmapSpec{'D0'} = mat2json(dist.getParam(1).paramValue);
        Kmarks = dist.getNumberOfTypes;
        d1k = cell(1, Kmarks);
        for km = 1:Kmarks
            d1k{km} = mat2json(dist.getParam(2+km).paramValue);
        end
        mmapSpec{'D1k'} = d1k;
        d{'mmap'} = mmapSpec;
    case 'MMPP2'
        d{'type'} = 'MMPP2';
        params = configureDictionary('string','cell');
        params{'lambda0'} = dist.getParam(1).paramValue;
        params{'lambda1'} = dist.getParam(2).paramValue;
        params{'sigma0'} = dist.getParam(3).paramValue;
        params{'sigma1'} = dist.getParam(4).paramValue;
        d{'params'} = params;
    case 'NHPP'
        d{'type'} = 'NHPP';
        params = configureDictionary('string','cell');
        % num2cell keeps a single-segment (1x1) rate vector from collapsing to
        % a JSON scalar, which the Gson reader would reject.
        nhppBp = double(dist.getBreakpoints());
        nhppRt = double(dist.getRates());
        params{'breakpoints'} = num2cell(nhppBp(:)');
        params{'rates'} = num2cell(nhppRt(:)');
        params{'cyclic'} = dist.isCyclic();
        d{'params'} = params;
    case {'MAPt','PHt'}
        % Segment matrices go out as a cell of matrices, one per segment, so a
        % single-segment schedule keeps its nesting instead of collapsing.
        d{'type'} = class(dist);
        params = configureDictionary('string','cell');
        schedBp = double(dist.getBreakpoints());
        params{'breakpoints'} = num2cell(schedBp(:)');
        if isa(dist,'MAPt')
            params{'D0'} = dist.getD0Segments();
            params{'D1'} = dist.getD1Segments();
        else
            aseg = dist.getAlphaSegments();
            for zz = 1:numel(aseg)
                aseg{zz} = aseg{zz}(:)';
            end
            params{'alpha'} = aseg;
            params{'S'} = dist.getSSegments();
        end
        params{'cyclic'} = dist.isCyclic();
        d{'params'} = params;
    case 'DiscreteSampler'
        d{'type'} = 'DiscreteSampler';
        params = configureDictionary('string','cell');
        params{'p'} = dist.getParam(1).paramValue(:)';
        params{'x'} = dist.getParam(2).paramValue(:)';
        d{'params'} = params;
    case 'Replayer'
        d{'type'} = 'Replayer';
        params = configureDictionary('string','cell');
        params{'fileName'} = dist.getParam(1).paramValue;
        try
            params{'mean'} = dist.getMean();
        catch
        end
        d{'params'} = params;
        % Save APH fit as fallback
        try
            aphDist = dist.fitAPH();
            if ~isempty(aphDist) && isa(aphDist, 'Distribution')
                ph = configureDictionary('string','cell');
                alpha = aphDist.getParam(1).paramValue;
                T = aphDist.getParam(2).paramValue;
                if isvector(alpha)
                    ph{'alpha'} = num2cell(double(full(alpha(:)')));
                else
                    ph{'alpha'} = mat2json(alpha);
                end
                ph{'T'} = mat2json(T);
                d{'ph'} = ph;
            end
        catch
        end
    case 'Prior'
        d{'type'} = 'Prior';
        if dist.isContinuousPrior()
            % getNumAlternatives is NaN here, so the discrete branch below would
            % emit an EMPTY alternative set -- see _kb/09-ldes-and-cache.md
            d{'kind'} = 'continuous';
            d{'paramDist'} = dist2json(dist.paramDist);
            d{'factory'} = prior_factory2json(dist);
        else
            d{'kind'} = 'discrete';
            alts = {};
            for ai = 1:dist.getNumAlternatives()
                altDist = dist.getAlternative(ai);
                altJson = dist2json(altDist);
                if ~isempty(altJson)
                    alts{end+1} = altJson; %#ok<AGROW>
                end
            end
            d{'distributions'} = alts;
            d{'probabilities'} = dist.probabilities(:)';
        end
    otherwise
        % Unmatched type: warn and emit real type name + mean/SCV, not a silently mislabeled Exp -- see _kb/09-ldes-and-cache.md
        line_warning(mfilename, sprintf(['Distribution "%s" has no JSON representation; ' ...
            'saving its mean and SCV only. The reloaded model will use an APH fitted ' ...
            'to those two moments.\n'], cn));
        m = dist.getMean();
        s = dist.getSCV();
        d{'type'} = cn;
        params = configureDictionary('string','cell');
        params{'mean'} = m;
        params{'scv'} = s;
        d{'params'} = params;
end
end


% =========================================================================
%  JSON encoding
% =========================================================================

function f = prior_factory2json(prior)
% Encode a continuous Prior's factory as the distribution it BUILDS.
%
% The handle theta -> Distribution cannot cross JSON, but the distribution it
% returns can, with the parameter left in the slot the handle fills: probing at
% two distinct quantiles of the parameter law identifies that slot, and a reader
% that substitutes theta there rebuilds the same alternative at ANY node count,
% so options.samples still means what it means in MATLAB. A factory that
% transforms the parameter instead of placing it (@(m) Exp(1/m)) has no such
% slot and is refused by name -- see _kb/09-ldes-and-cache.md.
tol = 1e-12;
t1 = Prior.quantile(prior.paramDist, 1/3);
t2 = Prior.quantile(prior.paramDist, 2/3);
if ~(t2 > t1 * (1 + tol))
    line_error(mfilename, ['A continuous Prior whose parameter law is concentrated at one ' ...
        'point cannot be encoded: the two probes of the factory coincide, so the parameter ' ...
        'slot cannot be identified. Use the discrete form Prior({dist}, 1) instead.']);
end
d1 = prior.distFactory(t1);
d2 = prior.distFactory(t2);
if ~isa(d1, 'Distribution') || ~isa(d2, 'Distribution')
    line_error(mfilename, 'The distribution factory of a Prior must return a Distribution object');
end
j1 = dist2json(d1);
j2 = dist2json(d2);
if ~strcmp(char(j1{'type'}), char(j2{'type'}))
    line_error(mfilename, sprintf(['The distribution factory of a Prior returned a %s at one ' ...
        'parameter value and a %s at another; a factory that switches family cannot be ' ...
        'encoded. Split the model into one Prior alternative per family.'], ...
        char(j1{'type'}), char(j2{'type'})));
end
if ~isKey(j1, 'params') || ~isKey(j2, 'params')
    line_error(mfilename, sprintf(['A %s carries its parameters as a fitted or phase-type ' ...
        'representation, not as named parameters, so a Prior factory returning one has no ' ...
        'parameter slot to encode. Use a family with explicit parameters (Exp, Erlang, Det, ' ...
        'Gamma, Uniform, Pareto, Weibull, Lognormal).'], char(j1{'type'})));
end
p1 = j1{'params'};
p2 = j2{'params'};
ks = sort(keys(p1));
slots = {};
for i = 1:numel(ks)
    key = ks(i);
    v1 = p1{key};
    v2 = p2{key};
    if isscalar(v1) && isscalar(v2) && isnumeric(v1) && isnumeric(v2) ...
            && abs(v1 - t1) <= tol * max(1, abs(t1)) && abs(v2 - t2) <= tol * max(1, abs(t2))
        slots{end+1} = char(key); %#ok<AGROW>
        continue;
    end
    same = isequal(size(v1), size(v2)) && isnumeric(v1) && isnumeric(v2) ...
        && all(abs(v1(:) - v2(:)) <= tol * max(1, abs(v1(:))));
    if ~same
        line_error(mfilename, sprintf(['The distribution factory of a Prior makes the "%s" ' ...
            'parameter of %s a function of the parameter rather than the parameter itself. ' ...
            'Only a slot filled verbatim can be encoded: reparameterize the Prior on that ' ...
            'slot, e.g. Prior(rateLaw, @(lambda) Exp(lambda)) instead of ' ...
            'Prior(meanLaw, @(m) Exp(1/m)).'], char(key), char(j1{'type'})));
    end
end
if isempty(slots)
    line_error(mfilename, sprintf(['The distribution factory of a Prior ignores its parameter: ' ...
        'every parameter of the returned %s is the same at two different parameter values. ' ...
        'Use the discrete form if the alternative does not depend on the parameter.'], ...
        char(j1{'type'})));
end
f = configureDictionary('string', 'cell');
f{'template'} = j1;
f{'slots'} = slots;
end

function s = encode_value(val, indent)
% Recursively encode a MATLAB value to JSON string.
if nargin < 2, indent = 0; end
pad  = repmat(' ', 1, indent);
pad2 = repmat(' ', 1, indent + 2);

if isa(val, 'dictionary')
    % sort: a dictionary keeps insertion order, containers.Map sorted its keys,
    % and the emitted key order is part of the on-the-wire model.json format
    ks = sort(keys(val));
    if isempty(ks)
        s = '{}';
    else
        parts = cell(1, length(ks));
        for i = 1:length(ks)
            k = ks(i);
            v = val{k};
            parts{i} = sprintf('%s"%s": %s', pad2, json_escape(char(k)), encode_value(v, indent + 2));
        end
        s = sprintf('{\n%s\n%s}', strjoin(parts, sprintf(',\n')), pad);
    end
elseif ischar(val) || isstring(val)
    s = sprintf('"%s"', json_escape(char(val)));
elseif islogical(val) && isscalar(val)
    if val, s = 'true'; else, s = 'false'; end
elseif isnumeric(val) && isscalar(val)
    if isnan(val)
        s = 'null';
    elseif isinf(val)
        if val > 0, s = '"Infinity"'; else, s = '"-Infinity"'; end
    elseif val == floor(val) && abs(val) < 1e15
        s = sprintf('%d', val);
    else
        % ROUND-TRIP OR NOTHING. %.15g drops the last two bits of a double, and
        % a reader that reloads the model then solves a DIFFERENT one: the C++
        % fluid solver answered cqn_repairmen with QLen 2.22225 from the
        % 15-digit rate 0.666666666666667 against 2.22217 from the exact
        % 0.66666666666666663. Widen only as far as the value needs, so the
        % common case stays short.
        s = sprintf('%.15g', val);
        if sscanf(s, '%f') ~= val
            s = sprintf('%.16g', val);
            if sscanf(s, '%f') ~= val
                s = sprintf('%.17g', val);
            end
        end
    end
elseif isnumeric(val) && isvector(val) && ~isscalar(val)
    parts = cell(1, length(val));
    for i = 1:length(val)
        parts{i} = encode_value(val(i), 0);
    end
    s = ['[', strjoin(parts, ', '), ']'];
elseif isnumeric(val) && ismatrix(val) && ~isvector(val)
    rows = cell(1, size(val, 1));
    for i = 1:size(val, 1)
        rows{i} = encode_value(val(i,:), 0);
    end
    s = ['[', strjoin(rows, ', '), ']'];
elseif iscell(val)
    if isempty(val)
        s = '[]';
    else
        parts = cell(1, length(val));
        for i = 1:length(val)
            parts{i} = sprintf('%s%s', pad2, encode_value(val{i}, indent + 2));
        end
        s = sprintf('[\n%s\n%s]', strjoin(parts, sprintf(',\n')), pad);
    end
elseif isstruct(val) && isscalar(val)
    fnames = fieldnames(val);
    if isempty(fnames)
        s = '{}';
    else
        parts = cell(1, length(fnames));
        for i = 1:length(fnames)
            fn = fnames{i};
            fv = val.(fn);
            parts{i} = sprintf('%s"%s": %s', pad2, json_escape(fn), encode_value(fv, indent + 2));
        end
        s = sprintf('{\n%s\n%s}', strjoin(parts, sprintf(',\n')), pad);
    end
else
    s = 'null';
end
end

function s = json_escape(str)
% Escape special characters for JSON strings.
s = strrep(str, '\', '\\');
s = strrep(s, '"', '\"');
s = strrep(s, sprintf('\n'), '\n');
s = strrep(s, sprintf('\r'), '\r');
s = strrep(s, sprintf('\t'), '\t');
end


% =========================================================================
%  Helper functions
% =========================================================================

function s = node_type_str(node)
% Get the JSON node type string for a node object.
if isa(node, 'Source'),      s = 'Source';
elseif isa(node, 'Sink'),   s = 'Sink';
elseif isa(node, 'Delay'),  s = 'Delay';
elseif isa(node, 'Cache'),  s = 'Cache';
elseif isa(node, 'Place'),  s = 'Place';
elseif isa(node, 'Transition'), s = 'Transition';
elseif isa(node, 'Queue'),  s = 'Queue';
elseif isa(node, 'Fork'),   s = 'Fork';
elseif isa(node, 'Join'),   s = 'Join';
elseif isa(node, 'Router'), s = 'Router';
elseif isa(node, 'ClassSwitch'), s = 'ClassSwitch';
elseif isa(node, 'Logger'), s = 'Logger';
else
    % No silent default. The former `else s = 'Queue'` saved a Logger as a
    % Queue -- a node with a service process and a buffer, which a Logger has
    % neither of -- so the model reloaded with an extra station and no trace.
    % Same failure mode, and same fix, as the sched_id_to_str whitelist above.
    line_error(mfilename, sprintf(['Node "%s" is of class %s, which has no model.json node ' ...
        'type; add it to node_type_str rather than letting it default.'], ...
        node.getName(), class(node)));
end
end

function s = sched_id_to_str(id)
% Map a SchedStrategy numeric ID to the wire enum name.
%
% The wire carries enum NAMES, uppercased, matching the JAR enum constants
% (jline.lang.constant.SchedStrategy) one-for-one for all 40 strategies. Do not
% reintroduce a hand-rolled whitelist here: the previous one covered 23 of 40 and
% silently degraded SRPT/FSP/EDD/EDF/FB/LCFSPI/PSJF/LRPT/SETF/LPS/... to FCFS.
% SchedStrategy.toText errors on an unknown id rather than inventing a default.
s = upper(SchedStrategy.toText(SchedStrategy.toId(id)));
end

function s = repl_to_str(id)
% Map ReplacementStrategy numeric ID to the wire enum name.
if id == ReplacementStrategy.LRU,   s = 'LRU';
elseif id == ReplacementStrategy.FIFO, s = 'FIFO';
elseif id == ReplacementStrategy.RR, s = 'RR';
elseif id == ReplacementStrategy.SFIFO, s = 'SFIFO';
elseif id == ReplacementStrategy.HLRU, s = 'HLRU';
elseif id == ReplacementStrategy.CLIMB, s = 'CLIMB';
elseif id == ReplacementStrategy.QLRU, s = 'QLRU';
else
    line_error(mfilename, sprintf('Unrecognized replacement strategy id %d.', id));
end
end

function s = depdisc_to_str(id)
% Map a DepartureDiscipline numeric ID to the wire name. The names are the JAR
% enum constants (jline.lang.constant.DepartureDiscipline), which is what the
% JAR writer emits and what its reader matches case-insensitively.
if id == DepartureDiscipline.NORMAL,    s = 'Normal';
elseif id == DepartureDiscipline.FIFO,  s = 'FIFO';
else
    line_error(mfilename, sprintf('Unrecognized departure discipline id %d.', id));
end
end

function s = droprule_to_str(id)
% Map DropStrategy numeric ID to schema-compatible string.
if id == DropStrategy.DROP,        s = 'drop';
elseif id == DropStrategy.WAITQ,   s = 'waitingQueue';
elseif id == DropStrategy.BAS,     s = 'blockingAfterService';
elseif id == DropStrategy.RETRIAL, s = 'retrial';
elseif id == DropStrategy.RETRIAL_WITH_LIMIT, s = 'retrialWithLimit';
else,                              s = '';
end
end

function rateMap = oi_rate_table(muFun, maxc, K)
% Build the OI/PAS macrostate rate table: for every per-class count vector cnt
% on the box lattice 0 <= cnt(r) <= maxc(r), evaluate mu on a canonical ordered
% microstate holding cnt(r) copies of class r (any ordering is valid since mu is
% order-independent). Keyed by the comma-joined 0-based class counts, matching
% the JAR reader (LineModelIO.oiServiceRate). The empty state is omitted.
rateMap = configureDictionary('string','cell');
shp = maxc + 1;
total = prod(shp);
for i = 1:total
    li = i - 1; cnt = zeros(1, K);
    for d = 1:K, cnt(d) = mod(li, shp(d)); li = floor(li / shp(d)); end
    if sum(cnt) == 0, continue, end
    micro = repelem(1:K, cnt);
    rate = muFun(micro);
    if ~isfinite(rate), rate = 0; end
    key = strjoin(arrayfun(@(x) sprintf('%d', x), cnt, 'UniformOutput', false), ',');
    rateMap{key} = rate;
end
end

function tbl = cd_scaling_table(beta, maxc, K)
% Materialize the class-dependence handle beta(n) over the box lattice
% 0 <= n(r) <= maxc(r). Keyed by the comma-joined 0-based per-class counts,
% matching the JAR reader (LineModelIO, "classDependence"). The handle may
% return a scalar (one scaling shared by every class) or a length-K vector of
% per-class scalings; the scalar form is broadcast to K entries here so that the
% reader is uniform and need not re-derive which form was used.
tbl = configureDictionary('string','cell');
shp = maxc + 1;
total = prod(shp);
for i = 1:total
    li = i - 1; n = zeros(1, K);
    for d = 1:K, n(d) = mod(li, shp(d)); li = floor(li / shp(d)); end
    v = beta(n);
    if isscalar(v)
        v = repmat(double(v), 1, K);
    else
        v = double(v(:)');
    end
    v(~isfinite(v)) = 0;
    key = strjoin(arrayfun(@(x) sprintf('%d', x), n, 'UniformOutput', false), ',');
    tbl{key} = num2cell(v);
end
end

function blk = gd_block(model)
% Materialize the network-level global (Whittle) dependence phi(n) onto the wire.
% n is the full (nstations x nclasses) population matrix, but only the slots a
% class can occupy carry a coordinate: a Source holds no jobs and a class with
% zero per-class capacity at a station never appears there, so those entries are
% pinned to 0. The lattice is the BOX over the remaining slots, which guarantees
% that clamping any queried population lands on a tabulated key.
sn = model.getStruct();
nodes = model.getNodes();
classes = model.getClasses();
M = sn.nstations; K = sn.nclasses;
phi = model.getGlobalDependence();
peak = model.getGlobalDependencePeak();
wcut = model.getGlobalDependenceCutoff();

stationNames = cell(1,M); classNames = cell(1,K);
for i = 1:M, stationNames{i} = nodes{sn.stationToNode(i)}.name; end
for r = 1:K, classNames{r} = classes{r}.name; end

slotSt = []; slotCl = []; cuts = [];
for i = 1:M
    if sn.nodetype(sn.stationToNode(i)) == NodeType.Source
        continue
    end
    for r = 1:K
        cap = sn.classcap(i,r);
        if ~(cap > 0)
            continue
        end
        if isfinite(sn.njobs(r))
            c = round(sn.njobs(r));
        else
            c = wcut;
        end
        if isfinite(cap)
            c = min(c, round(cap));
        end
        if c < 0, c = 0; end
        slotSt(end+1) = i; %#ok<AGROW>
        slotCl(end+1) = r; %#ok<AGROW>
        cuts(end+1) = c;   %#ok<AGROW>
    end
end

total = prod(double(cuts) + 1);
maxEntries = 200000;
if total > maxEntries
    line_error(mfilename, sprintf('The global dependence lattice has %g points (%d varying station-class slots with cutoffs %s), above the wire limit of %d. Lower the third argument of setGlobalDependence, or solve the model natively.', total, numel(cuts), mat2str(cuts), maxEntries));
end

slots = {};
for s = 1:numel(cuts)
    sm = configureDictionary('string','cell');
    sm{'station'} = stationNames{slotSt(s)};
    sm{'class'} = classNames{slotCl(s)};
    slots{end+1} = sm; %#ok<AGROW>
end

blk = configureDictionary('string','cell');
blk{'type'} = 'globalDependent';
blk{'stations'} = stationNames;
blk{'classes'} = classNames;
blk{'slots'} = slots;
blk{'cutoffs'} = num2cell(double(cuts(:)'));
blk{'cutoff'} = double(wcut);
blk{'scaling'} = gd_scaling_table(phi, slotSt, slotCl, cuts, M, K);
pk = double(peak);
blk{'peak'} = num2cell(reshape(pk', 1, []));
end

function tbl = gd_scaling_table(phi, slotSt, slotCl, cuts, M, K)
% Tabulate phi over the slot box lattice. Key: the comma-joined 0-based slot
% counts, in slot order. Value: the FULL (M x K) scaling flattened row-major, so
% the reader restores the matrix the handle returned without re-deriving which of
% the scalar / column / matrix return forms was used.
tbl = configureDictionary('string','cell');
P = numel(cuts);
shp = double(cuts) + 1;
total = prod(shp);
for i = 1:total
    li = i - 1; c = zeros(1, P);
    for d = 1:P, c(d) = mod(li, shp(d)); li = floor(li / shp(d)); end
    n = zeros(M, K);
    for d = 1:P, n(slotSt(d), slotCl(d)) = c(d); end
    v = phi(n);
    if isscalar(v)
        v = double(v) * ones(M, K);
    elseif isequal(size(v), [M, 1])
        v = repmat(double(v(:)), 1, K);
    else
        v = double(v);
    end
    v(~isfinite(v)) = 0;
    if P == 0
        key = '0';   % degenerate: no varying slot, phi is a constant
    else
        key = strjoin(arrayfun(@(x) sprintf('%d', x), c, 'UniformOutput', false), ',');
    end
    tbl{key} = num2cell(reshape(v', 1, []));
end
end

function tbl = firingdep_scaling_table(g, slotIdx, caps, nnodes, nclasses)
% Materialize the firing-rate dependence handle g(M) over the box lattice of the
% enabling (place,class) slots, 0 <= count(s) <= caps(s). M is the full
% node-indexed marking matrix (nnodes x nclasses); only the enabling slots vary,
% all other entries are held at 0. Keyed by the comma-joined 0-based slot counts,
% matching the JAR/Python readers. g returns a positive scalar multiplier.
tbl = configureDictionary('string','cell');
P = size(slotIdx, 1);
shp = caps + 1;
total = prod(shp);
for i = 1:total
    li = i - 1; c = zeros(1, P);
    for d = 1:P, c(d) = mod(li, shp(d)); li = floor(li / shp(d)); end
    M = zeros(nnodes, nclasses);
    for d = 1:P, M(slotIdx(d,1), slotIdx(d,2)) = c(d); end
    v = double(g(M));
    if ~isscalar(v), v = v(1); end
    if ~isfinite(v), v = 0; end
    key = strjoin(arrayfun(@(x) sprintf('%d', x), c, 'UniformOutput', false), ',');
    tbl{key} = v;
end
end

function s = prectype_to_str(id)
% Map ActivityPrecedenceType numeric ID to JAR-compatible string.
if id == ActivityPrecedenceType.PRE_SEQ,        s = 'pre';
elseif id == ActivityPrecedenceType.PRE_AND,    s = 'pre-AND';
elseif id == ActivityPrecedenceType.PRE_OR,     s = 'pre-OR';
elseif id == ActivityPrecedenceType.POST_SEQ,   s = 'post';
elseif id == ActivityPrecedenceType.POST_AND,   s = 'post-AND';
elseif id == ActivityPrecedenceType.POST_OR,    s = 'post-OR';
elseif id == ActivityPrecedenceType.POST_LOOP,  s = 'post-LOOP';
elseif id == ActivityPrecedenceType.POST_CACHE, s = 'post-CACHE';
else,                                           s = 'pre';
end
end

function c = mat2json(M)
% C = MAT2JSON(M)
% Cell of row cells, so jsonencode writes a 2-D ARRAY at every size. A bare
% numeric matrix collapses to a scalar at 1x1 and to a flat array at 1xn, and
% the readers on the other side of the interchange (line-cli, the Python
% loader) require nested arrays: a one-phase APH written as "T": -1 aborts the
% C++ reader with 'type must be array, but is number'.
M = double(full(M));
c = cell(1, size(M, 1));
for rr = 1:size(M, 1)
    c{rr} = num2cell(M(rr, :));
end
end

function tf = isIdentityClassSwitch(csm, K)
% True if the classSwitchMatrix CSM is the K x K identity (or empty/unset,
% which defaults to identity). A non-identity matrix carries the class switch
% and must be collapsed to same-class routing on export; an identity matrix
% means the switch is expressed through the routing and must be kept verbatim.
tf = true;
if isempty(csm)
    return;
end
for rr = 1:min(K, size(csm,1))
    for ss = 1:min(K, size(csm,2))
        if abs(csm(rr,ss) - double(rr==ss)) > 1e-12
            tf = false;
            return;
        end
    end
end
end
