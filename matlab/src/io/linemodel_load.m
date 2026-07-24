function model = linemodel_load(filename)
% LINEMODEL_LOAD Load a LINE model from JSON.
%
%   MODEL = LINEMODEL_LOAD(FILENAME) loads a model from the specified JSON
%   file (conforming to line-model.schema.json) and returns a Network,
%   LayeredNetwork, Workflow, or Environment object.
%
% Parameters:
%   filename - path to a .json file
%
% Returns:
%   model - Network, LayeredNetwork, Workflow, or Environment object
%
% Example:
%   model = linemodel_load('mm1.json');
%   solver = SolverMVA(model);
%   AvgTable = solver.getAvgTable();
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

jsonText = fileread(filename);
doc = jsondecode(jsonText);

if ~isfield(doc, 'model')
    error('linemodel_load:noModel', 'JSON file does not contain a "model" field.');
end

data = doc.model;
mtype = data.type;

switch mtype
    case 'Network'
        model = json2network(data, jsonText);
    case 'LayeredNetwork'
        model = json2layered(data);
    case 'Workflow'
        model = json2workflow(data);
    case 'Environment'
        model = json2environment(data, jsonText);
    otherwise
        error('linemodel_load:unknownType', 'Unsupported model type: %s', mtype);
end
end


% =========================================================================
%  Network deserialization
% =========================================================================

function model = json2network(data, rawJson)
% Reconstruct a Network from decoded JSON struct.
% rawJson is the original text, used for parsing routing keys with commas.

modelName = 'model';
if isfield(data, 'name')
    modelName = data.name;
end
model = Network(modelName);

% --- Create nodes (before classes, since ClosedClass needs refstat) ---
nodeList = {};
if isfield(data, 'nodes')
    nds = data.nodes;
    if isstruct(nds)
        nds = num2cell(nds);
    end
    for i = 1:length(nds)
        nd = nds{i};
        if isstruct(nd)
            nd_name = nd.name;
            nd_type = nd.type;
        else
            nd_name = nd('name');
            nd_type = nd('type');
        end
        node = create_node(model, nd, nd_name, nd_type);
        nodeList{end+1} = node; %#ok<AGROW>
    end
end
node_map = containers.Map();
for i = 1:length(nodeList)
    node_map(nodeList{i}.name) = nodeList{i};
end

% --- Deferred node linking (Fork/Join, Fork tasksPerLink) ---
if isfield(data, 'nodes')
    nds2 = data.nodes;
    if isstruct(nds2)
        nds2 = num2cell(nds2);
    end
    for i = 1:length(nds2)
        nd2 = nds2{i};
        nd_name2 = nd2.name;
        node2 = node_map(nd_name2);
        % Join: link to paired Fork
        if isfield(nd2, 'forkNode') && isa(node2, 'Join')
            if node_map.isKey(nd2.forkNode)
                node2.joinOf = node_map(nd2.forkNode);
            end
        end
        % Fork: set tasksPerLink
        if isfield(nd2, 'tasksPerLink') && isa(node2, 'Fork')
            node2.setTasksPerLink(nd2.tasksPerLink);
        end
    end
end

% --- Create classes ---
classesList = {};
if isfield(data, 'classes')
    cls = data.classes;
    if isstruct(cls)
        cls = num2cell(cls);
    end
    for i = 1:length(cls)
        cd = cls{i};
        cname = cd.name;
        ctype = cd.type;
        switch ctype
            case 'Open'
                prio = 0;
                if isfield(cd, 'priority'), prio = cd.priority; end
                jc = OpenClass(model, cname, prio);
            case {'Closed', 'SelfLooping'}
                pop = cd.population;
                refNode = [];
                if isfield(cd, 'refNode') && node_map.isKey(cd.refNode)
                    refNode = node_map(cd.refNode);
                end
                prio = 0;
                if isfield(cd, 'priority'), prio = cd.priority; end
                if isempty(refNode)
                    error('linemodel_load:noRefNode', ...
                        '%s class "%s" has no valid refNode.', ctype, cname);
                end
                if strcmp(ctype, 'SelfLooping')
                    jc = SelfLoopingClass(model, cname, pop, refNode, prio);
                else
                    jc = ClosedClass(model, cname, pop, refNode, prio);
                end
            case 'Signal'
                prio = 0;
                if isfield(cd, 'priority'), prio = cd.priority; end
                sigType = SignalType.NEGATIVE;
                if isfield(cd, 'signalType')
                    sigType = SignalType.fromText(cd.signalType);
                end
                if ~isfield(cd, 'openOrClosed')
                    % Bare Signal: neither open nor closed, marked by omitting 'openOrClosed' -- see _kb/09-ldes-and-cache.md
                    jc = Signal(model, cname, sigType, prio);
                elseif strcmp(cd.openOrClosed, 'Closed')
                    refNode = [];
                    if isfield(cd, 'refNode') && node_map.isKey(cd.refNode)
                        refNode = node_map(cd.refNode);
                    end
                    if isempty(refNode)
                        error('linemodel_load:noRefNode', ...
                            'ClosedSignal "%s" has no valid refNode.', cname);
                    end
                    jc = ClosedSignal(model, cname, sigType, refNode, prio);
                else
                    jc = OpenSignal(model, cname, sigType, prio);
                end
                % Removal distribution
                if isfield(cd, 'removalDistribution')
                    remDist = json2dist(cd.removalDistribution);
                    if ~isempty(remDist)
                        jc.setRemovalDistribution(remDist);
                    end
                end
                % Removal policy
                if isfield(cd, 'removalPolicy')
                    jc.setRemovalPolicy(RemovalPolicy.fromText(cd.removalPolicy));
                end
            otherwise
                jc = OpenClass(model, cname);
        end
        if isfield(cd, 'deadline') && isfinite(cd.deadline)
            jc.deadline = cd.deadline;
        end
        if isfield(cd, 'isReferenceClass') && cd.isReferenceClass
            jc.setReferenceClass(true);
        end
        % Class-level (global) patience. A node-scoped 'patience' entry, restored
        % later, overrides this for the station it names.
        if isfield(cd, 'patience') && ~isempty(cd.patience)
            patDist = json2dist(cd.patience);
            if ~isempty(patDist) && ~isa(patDist, 'Disabled')
                if isfield(cd, 'impatienceType')
                    jc.setPatience(str_to_impatience(cd.impatienceType), patDist);
                else
                    jc.setPatience(patDist);
                end
            end
        end
        classesList{end+1} = jc; %#ok<AGROW>
    end
end
class_map = containers.Map();
for i = 1:length(classesList)
    class_map(classesList{i}.name) = classesList{i};
end

% --- Resolve signal targetClass associations ---
if isfield(data, 'classes')
    cls2 = data.classes;
    if isstruct(cls2), cls2 = num2cell(cls2); end
    for i = 1:length(cls2)
        cd2 = cls2{i};
        if isfield(cd2, 'type') && strcmp(cd2.type, 'Signal') && isfield(cd2, 'targetClass')
            if class_map.isKey(cd2.name) && class_map.isKey(cd2.targetClass)
                sigCls = class_map(cd2.name);
                sigCls.forJobClass(class_map(cd2.targetClass));
            end
        end
        % Reply signal binding (sn.syncreply). Resolved in a second pass because
        % the reply class may be declared after the class that references it.
        if isfield(cd2, 'replySignalClass') && class_map.isKey(cd2.name) ...
                && class_map.isKey(cd2.replySignalClass)
            class_map(cd2.name).setReplySignalClass(class_map(cd2.replySignalClass));
        end
    end
end

% --- Set service/arrival distributions ---
if isfield(data, 'nodes')
    nds = data.nodes;
    if isstruct(nds)
        nds = num2cell(nds);
    end
    for i = 1:length(nds)
        nd = nds{i};
        nd_name = nd.name;
        node = node_map(nd_name);

        % OI/PAS queues are parameterized by oiServiceRate below, not per-class distributions -- see _kb/09-ldes-and-cache.md
        isOIQueue = isa(node, 'Queue') && ~isa(node, 'Delay') && ...
            (node.schedStrategy == SchedStrategy.PAS || node.schedStrategy == SchedStrategy.OI);
        if isfield(nd, 'service') && ~isempty(nd.service) && ~isOIQueue
            svc = nd.service;
            svcFields = fieldnames(svc);
            for f = 1:length(svcFields)
                cname = svcFields{f};
                distJson = svc.(cname);
                if ~class_map.isKey(cname)
                    continue;
                end
                jc = class_map(cname);
                dist = json2dist(distJson);
                if ~isempty(dist)
                    if isa(node, 'Source')
                        node.setArrival(jc, dist);
                    elseif isa(node, 'Place')
                        % setService turns the Place into a queueing Place, installing the server section from create_node's strategy
                        node.setService(jc, dist);
                    elseif isa(node, 'Queue') || isa(node, 'Delay')
                        % Pass DPS/GPS weight if available
                        weight = 1;
                        if isfield(nd, 'schedParams')
                            sp_tmp = nd.schedParams;
                            if isfield(sp_tmp, cname)
                                weight = sp_tmp.(cname);
                            end
                        end
                        node.setService(jc, dist, weight);
                    end
                end
            end
        end

        % Batch arrivals, applied after setArrival because setArrivalBatch
        % validates the batch law independently of the interarrival process.
        if isfield(nd, 'arrivalBatch') && ~isempty(nd.arrivalBatch) && isa(node, 'Source')
            batchFields = fieldnames(nd.arrivalBatch);
            for f = 1:numel(batchFields)
                cname = batchFields{f};
                if ~class_map.isKey(cname)
                    continue;
                end
                bdist = json2dist(nd.arrivalBatch.(cname));
                if ~isempty(bdist)
                    node.setArrivalBatch(class_map(cname), bdist);
                end
            end
        end

        % Marked (MMAP) arrival: rebind the shared MarkedMAP so mark k drives the k-th listed class
        if isfield(nd, 'markedClasses') && ~isempty(nd.markedClasses) && isa(node, 'Source')
            markedNames = nd.markedClasses;
            if ischar(markedNames), markedNames = {markedNames}; end
            markedList = {};
            for f = 1:numel(markedNames)
                mnm = markedNames{f};
                if class_map.isKey(mnm)
                    markedList{end+1} = class_map(mnm); %#ok<AGROW>
                end
            end
            if ~isempty(markedList)
                firstDist = node.getArrivalProcess(markedList{1}.index);
                if isa(firstDist, 'MarkedMAP')
                    node.setMarkedArrival(firstDist, markedList);
                end
            end
        end

        % ClassSwitch matrix (dict format: classSwitchMatrix)
        if isfield(nd, 'classSwitchMatrix') && isa(node, 'ClassSwitch')
            csm_data = nd.classSwitchMatrix;
            classes = model.getClasses();
            K = length(classes);
            mat = zeros(K);
            class_idx = containers.Map();
            for ci = 1:K
                class_idx(classes{ci}.name) = ci;
            end
            fromFields = fieldnames(csm_data);
            for fi = 1:length(fromFields)
                fromName = fromFields{fi};
                if ~class_idx.isKey(fromName), continue; end
                ri = class_idx(fromName);
                toStruct = csm_data.(fromName);
                toFields = fieldnames(toStruct);
                for ti = 1:length(toFields)
                    toName = toFields{ti};
                    if ~class_idx.isKey(toName), continue; end
                    ci = class_idx(toName);
                    mat(ri, ci) = toStruct.(toName);
                end
            end
            node.server = node.server.updateClassSwitch(mat);
        % Legacy 2D array format: csMatrix (from older JAR saves)
        elseif isfield(nd, 'csMatrix') && isa(node, 'ClassSwitch')
            mat = nd.csMatrix;
            if iscell(mat)
                mat = cell2mat(mat);
            end
            node.server = node.server.updateClassSwitch(mat);
        end

        % DPS scheduling parameters (weights already passed via setService above)

        % Departure discipline read after the service loop, which creates the departureDiscipline slots
        if isfield(nd, 'departureDiscipline') && isa(node, 'Place')
            ddData = nd.departureDiscipline;
            ddFields = fieldnames(ddData);
            for ddi = 1:length(ddFields)
                cname = ddFields{ddi};
                if class_map.isKey(cname)
                    node.setDepartureDiscipline(class_map(cname), ...
                        str_to_depdisc(ddData.(cname)));
                end
            end
        end

        % Per-class buffer capacity. Station, not Queue: a queueing Place has a
        % per-class capacity too and Place does not extend Queue.
        if isfield(nd, 'classCap') && isa(node, 'Station')
            ccData = nd.classCap;
            ccFields = fieldnames(ccData);
            for cci = 1:length(ccFields)
                cname = ccFields{cci};
                if class_map.isKey(cname)
                    jc = class_map(cname);
                    % Find class index
                    cls = model.getClasses();
                    for ci = 1:length(cls)
                        if strcmp(cls{ci}.name, cname)
                            node.classCap(ci) = ccData.(cname);
                            break;
                        end
                    end
                end
            end
        end

        % Drop rules
        if isfield(nd, 'dropRule') && isa(node, 'Station')
            drData = nd.dropRule;
            drFields = fieldnames(drData);
            for dri = 1:length(drFields)
                cname = drFields{dri};
                if class_map.isKey(cname)
                    cls = model.getClasses();
                    for ci = 1:length(cls)
                        if strcmp(cls{ci}.name, cname)
                            node.dropRule(ci) = str_to_droprule(drData.(cname));
                            break;
                        end
                    end
                end
            end
        end

        % OI/PAS mu(c) macrostate table + cutoffs (+ swap graph for PAS) -- see _kb/09-ldes-and-cache.md
        if isfield(nd, 'oiServiceRate') && isa(node, 'Queue') && ...
                (node.schedStrategy == SchedStrategy.PAS || node.schedStrategy == SchedStrategy.OI)
            % OI queues keep a fixed zero swap graph (setSwapGraph rejects them).
            if isfield(nd, 'swapGraph') && node.schedStrategy == SchedStrategy.PAS
                node.setSwapGraph(json2mat(nd.swapGraph));
            end
            if isfield(nd, 'oiCutoffs')
                oicut = round(double(oi_cutoffs_vec(nd.oiCutoffs)));
            else
                oicut = [];
            end
            node.setServiceRateFunction(oi_table_to_handle(nd.oiServiceRate, oicut));
        end

        % Load-dependent scaling
        if isfield(nd, 'loadDependence') && isa(node, 'Queue')
            ld = nd.loadDependence;
            if isfield(ld, 'type') && strcmp(ld.type, 'loadDependent') && isfield(ld, 'scaling')
                scaling = ld.scaling(:)';
                node.setLoadDependence(scaling);
            end
        end

        % Class-dependent scaling beta_{i,r}(n): rebuild the handle from the
        % materialized lattice table (see cd_scaling_table in linemodel_save).
        if isfield(nd, 'classDependence') && isa(node, 'Station')
            cdep = nd.classDependence;
            if isfield(cdep, 'type') && strcmp(cdep.type, 'classDependent') ...
                    && isfield(cdep, 'scaling')
                if isfield(cdep, 'cutoffs')
                    cdcut = round(double(cdep.cutoffs(:)'));
                else
                    cdcut = [];
                end
                cdHandle = cd_table_to_handle(cdep.scaling, cdcut);
                if isfield(cdep, 'peak') && ~isempty(cdep.peak)
                    if iscell(cdep.peak)
                        cdPeak = double(cell2mat(cdep.peak));
                    else
                        cdPeak = double(cdep.peak(:)');
                    end
                else
                    % Legacy JSON without an explicit peak: derive from the handle over the lattice, matching Util = T*S/peak
                    cdPeak = cd_peak_scaling(cdHandle, cdcut, numel(cdcut));
                end
                node.setLimitedClassDependence(cdHandle, cdPeak);
            end
        end

        % Joint-dependent scaling eta_i(n) (non-product-form): rebuild the
        % handle from the materialized lattice table, twin of classDependence.
        if isfield(nd, 'jointDependence') && isa(node, 'Station')
            jdep = nd.jointDependence;
            if isfield(jdep, 'type') && strcmp(jdep.type, 'jointDependent') ...
                    && isfield(jdep, 'scaling')
                if isfield(jdep, 'cutoffs')
                    jdcut = round(double(jdep.cutoffs(:)'));
                else
                    jdcut = [];
                end
                jdHandle = cd_table_to_handle(jdep.scaling, jdcut);
                if isfield(jdep, 'peak') && ~isempty(jdep.peak)
                    if iscell(jdep.peak)
                        jdPeak = double(cell2mat(jdep.peak));
                    else
                        jdPeak = double(jdep.peak(:)');
                    end
                else
                    jdPeak = cd_peak_scaling(jdHandle, jdcut, numel(jdcut));
                end
                node.setLimitedJointDependence(jdHandle, jdPeak);
            end
        end

        % Join strategy and quorum
        if isa(node, 'Join')
            if isfield(nd, 'joinStrategy')
                jsStr = nd.joinStrategy;
                classes = model.getClasses();
                for ci = 1:length(classes)
                    switch jsStr
                        case 'STD'
                            node.input.setStrategy(classes{ci}, JoinStrategy.STD);
                        case {'PARTIAL', 'QUORUM', 'Quorum'}
                            node.input.setStrategy(classes{ci}, JoinStrategy.PARTIAL);
                    end
                end
            end
            if isfield(nd, 'joinQuorum')
                jq = nd.joinQuorum;
                classes = model.getClasses();
                for ci = 1:length(classes)
                    node.input.setRequired(classes{ci}, jq);
                end
            end
        end

        % Cache hit/miss/popularity: accept both nested MATLAB (nd.cache.*) and flat JAR/Python schema
        if isa(node, 'Cache')
            if isfield(nd, 'cache')
                cc = nd.cache;
            else
                cc = struct();
                if isfield(nd, 'hitClass'), cc.hitClass = nd.hitClass; end
                if isfield(nd, 'missClass'), cc.missClass = nd.missClass; end
                if isfield(nd, 'popularity'), cc.popularity = nd.popularity; end
                if isfield(nd, 'accessGraph'), cc.accessGraph = nd.accessGraph; end
                if isfield(nd, 'accessProb'), cc.accessProb = nd.accessProb; end
                if isfield(nd, 'initialState'), cc.initialState = nd.initialState; end
            end
            % Prefer the nested copy, fall back to flat, so any of the three bridges' output loads the same way
            if ~isfield(cc, 'admissionProb') && isfield(nd, 'admissionProb')
                cc.admissionProb = nd.admissionProb;
            end
            % q-LRU admission probability on a miss
            if isfield(cc, 'admissionProb')
                node.setAdmissionProb(cc.admissionProb);
            end
            % Hit class mapping
            if isfield(cc, 'hitClass')
                hcData = cc.hitClass;
                hcFields = fieldnames(hcData);
                for hci = 1:length(hcFields)
                    inName = hcFields{hci};
                    outName = hcData.(inName);
                    if class_map.isKey(inName) && class_map.isKey(outName)
                        node.setHitClass(class_map(inName), class_map(outName));
                    end
                end
            end
            % Miss class mapping
            if isfield(cc, 'missClass')
                mcData = cc.missClass;
                mcFields = fieldnames(mcData);
                for mci = 1:length(mcFields)
                    inName = mcFields{mci};
                    outName = mcData.(inName);
                    if class_map.isKey(inName) && class_map.isKey(outName)
                        node.setMissClass(class_map(inName), class_map(outName));
                    end
                end
            end
            % Popularity distributions (setRead)
            if isfield(cc, 'popularity')
                popData = cc.popularity;
                popFields = fieldnames(popData);
                for pfi = 1:length(popFields)
                    cname = popFields{pfi};
                    if class_map.isKey(cname)
                        popDist = json2dist(popData.(cname));
                        if ~isempty(popDist) && ~isa(popDist, 'Disabled')
                            node.setRead(class_map(cname), popDist);
                        end
                    end
                end
            end
            % Access-cost (list-move) structure: per-item graph shared by all
            % classes, or full per-class accessProb matrices
            if isfield(cc, 'accessGraph')
                node.graph = json2matcell(cc.accessGraph);
            elseif isfield(cc, 'accessProb')
                apData = cc.accessProb;
                if ~iscell(apData)
                    % jsondecode collapsed the uniform [class][item][r][c]
                    % nesting into a 4D numeric array; re-split per class
                    apData = arrayfun(@(k1) squeeze(apData(k1, :, :, :)), ...
                        1:size(apData, 1), 'UniformOutput', false);
                end
                Kap = numel(apData);
                rows = cellfun(@json2matcell, apData, 'UniformOutput', false);
                Nap = max(cellfun(@numel, rows));
                R = cell(Kap, Nap);
                for k1 = 1:Kap
                    R(k1, 1:numel(rows{k1})) = rows{k1};
                end
                node.setAccessProb(R);
            end
            % Initial cache state [class counts | contents | retrieval bitmap]
            if isfield(cc, 'initialState')
                node.setState(cc.initialState(:)');
            end
        end
        % Heterogeneous server types
        if isa(node, 'Queue') && isfield(nd, 'serverTypes')
            stArr = nd.serverTypes;
            if isstruct(stArr), stArr = num2cell(stArr); end
            for si = 1:length(stArr)
                stData = stArr{si};
                stName = stData.name;
                stCount = stData.count;
                st = ServerType(stName, stCount);
                % Compatible classes
                if isfield(stData, 'compatibleClasses')
                    ccList = stData.compatibleClasses;
                    if ~iscell(ccList), ccList = {ccList}; end
                    for cci = 1:length(ccList)
                        if class_map.isKey(ccList{cci})
                            st.addCompatible(class_map(ccList{cci}));
                        end
                    end
                end
                node.addServerType(st);
                % Per-class service distributions
                if isfield(stData, 'service')
                    svcData = stData.service;
                    svcFields = fieldnames(svcData);
                    for fi = 1:length(svcFields)
                        cname = svcFields{fi};
                        if class_map.isKey(cname)
                            jc = class_map(cname);
                            dist = json2dist(svcData.(cname));
                            if ~isempty(dist)
                                node.setHeteroService(jc, st, dist);
                            end
                        end
                    end
                end
            end
            % Scheduling policy
            if isfield(nd, 'heteroSchedPolicy')
                policy = HeteroSchedPolicy.fromText(nd.heteroSchedPolicy);
                node.setHeteroSchedPolicy(policy);
            end
        end
    end
end

% --- Restore Balking, Retrial, Patience ---
if isfield(data, 'nodes')
    ndsImp = data.nodes;
    if isstruct(ndsImp), ndsImp = num2cell(ndsImp); end
    for i = 1:length(ndsImp)
        ndImp = ndsImp{i};
        if ~node_map.isKey(ndImp.name), continue; end
        node = node_map(ndImp.name);
        % Assigned directly, not via setStatePrior: its node.space length check cannot pass before the solver builds the state space
        if isfield(ndImp, 'statePrior') && isa(node, 'StatefulNode')
            node.statePrior = double(ndImp.statePrior(:));
        end
        if ~isa(node, 'Queue'), continue; end
        % Immediate feedback on self-loops, per class
        if isfield(ndImp, 'immediateFeedback') && ~isempty(ndImp.immediateFeedback)
            ifData = ndImp.immediateFeedback;
            ifNames = fieldnames(ifData);
            for fi = 1:length(ifNames)
                cname = ifNames{fi};
                if class_map.isKey(cname) && ifData.(cname)
                    node.setImmediateFeedback(class_map(cname));
                end
            end
        end
        % Orbit impatience (abandonment from the retrial orbit), per class
        if isfield(ndImp, 'orbitImpatience') && ~isempty(ndImp.orbitImpatience)
            orbData = ndImp.orbitImpatience;
            orbNames = fieldnames(orbData);
            for fi = 1:length(orbNames)
                cname = orbNames{fi};
                if ~class_map.isKey(cname), continue; end
                orbDist = json2dist(orbData.(cname));
                if ~isempty(orbDist) && ~isa(orbDist, 'Disabled')
                    node.setOrbitImpatience(class_map(cname), orbDist);
                end
            end
        end
        % Batch rejection probability (retrial queues), per class
        if isfield(ndImp, 'batchRejectProb') && ~isempty(ndImp.batchRejectProb)
            brpData = ndImp.batchRejectProb;
            brpNames = fieldnames(brpData);
            for fi = 1:length(brpNames)
                cname = brpNames{fi};
                if ~class_map.isKey(cname), continue; end
                node.setBatchRejectProbability(class_map(cname), brpData.(cname));
            end
        end
        % Balking
        if isfield(ndImp, 'balking') && ~isempty(ndImp.balking)
            balkData = ndImp.balking;
            fnames = fieldnames(balkData);
            for fi = 1:length(fnames)
                className = fnames{fi};
                if ~class_map.isKey(className), continue; end
                jc = class_map(className);
                bjc = balkData.(className);
                % Parse strategy
                switch bjc.strategy
                    case 'QUEUE_LENGTH', strategy = BalkingStrategy.QUEUE_LENGTH;
                    case 'EXPECTED_WAIT', strategy = BalkingStrategy.EXPECTED_WAIT;
                    case 'COMBINED', strategy = BalkingStrategy.COMBINED;
                    otherwise, continue;
                end
                % Parse thresholds
                thData = bjc.thresholds;
                if isstruct(thData), thData = num2cell(thData); end
                thresholds = {};
                for ti = 1:length(thData)
                    td = thData{ti};
                    maxJobs = td.maxJobs;
                    if maxJobs < 0, maxJobs = Inf; end
                    thresholds{end+1} = {td.minJobs, maxJobs, td.probability};
                end
                node.setBalking(jc, strategy, thresholds);
            end
        end
        % Retrial
        if isfield(ndImp, 'retrial') && ~isempty(ndImp.retrial)
            retData = ndImp.retrial;
            fnames = fieldnames(retData);
            for fi = 1:length(fnames)
                className = fnames{fi};
                if ~class_map.isKey(className), continue; end
                jc = class_map(className);
                rjc = retData.(className);
                delayDist = json2dist(rjc.delay);
                maxAttempts = -1;
                if isfield(rjc, 'maxAttempts')
                    maxAttempts = rjc.maxAttempts;
                end
                node.setRetrial(jc, delayDist, maxAttempts);
            end
        end
        % Patience
        if isfield(ndImp, 'patience') && ~isempty(ndImp.patience)
            patData = ndImp.patience;
            fnames = fieldnames(patData);
            for fi = 1:length(fnames)
                className = fnames{fi};
                if ~class_map.isKey(className), continue; end
                jc = class_map(className);
                pjc = patData.(className);
                patDist = json2dist(pjc.distribution);
                if isfield(pjc, 'impatienceType')
                    impType = str_to_impatience(pjc.impatienceType);
                else
                    impType = ImpatienceType.RENEGING;
                end
                node.setPatience(jc, impType, patDist);
            end
        end
    end
end

% --- Configure Transition modes ---
if isfield(data, 'nodes')
    nds3 = data.nodes;
    if isstruct(nds3)
        nds3 = num2cell(nds3);
    end
    for i = 1:length(nds3)
        nd3 = nds3{i};
        if ~isfield(nd3, 'modes'), continue; end
        if ~strcmp(nd3.type, 'Transition'), continue; end
        tnode = node_map(nd3.name);
        modesData = nd3.modes;
        if isstruct(modesData)
            modesData = num2cell(modesData);
        end
        for mi = 1:length(modesData)
            md = modesData{mi};
            modeName = 'Mode';
            if isfield(md, 'name'), modeName = md.name; end
            mode = tnode.addMode(modeName);
            % Distribution
            if isfield(md, 'distribution') && ~isempty(md.distribution)
                dist = json2dist(md.distribution);
                if ~isempty(dist)
                    tnode.setDistribution(mode, dist);
                end
            end
            % Timing strategy
            if isfield(md, 'timingStrategy')
                if strcmp(md.timingStrategy, 'IMMEDIATE')
                    tnode.setTimingStrategy(mode, TimingStrategy.IMMEDIATE);
                else
                    tnode.setTimingStrategy(mode, TimingStrategy.TIMED);
                end
            end
            % Number of servers
            if isfield(md, 'numServers')
                nsVal = md.numServers;
                if ischar(nsVal) || isstring(nsVal)
                    if strcmpi(nsVal, 'Infinity'), nsVal = Inf; else, nsVal = str2double(nsVal); end
                end
                if nsVal > 1
                    tnode.setNumberOfServers(mode, nsVal);
                end
            end
            % Firing priority
            if isfield(md, 'firingPriority')
                tnode.setFiringPriorities(mode, md.firingPriority);
            end
            % Firing weight
            if isfield(md, 'firingWeight')
                tnode.setFiringWeights(mode, md.firingWeight);
            end
            % Enabling conditions
            if isfield(md, 'enablingConditions')
                ecList = md.enablingConditions;
                if isstruct(ecList), ecList = num2cell(ecList); end
                for ei = 1:length(ecList)
                    ec = ecList{ei};
                    if node_map.isKey(ec.node) && class_map.isKey(ec.class)
                        tnode.setEnablingConditions(mode, class_map(ec.class), node_map(ec.node), ec.count);
                    end
                end
            end
            % Inhibiting conditions
            if isfield(md, 'inhibitingConditions')
                icList = md.inhibitingConditions;
                if isstruct(icList), icList = num2cell(icList); end
                for ii = 1:length(icList)
                    ic = icList{ii};
                    if node_map.isKey(ic.node) && class_map.isKey(ic.class)
                        tnode.setInhibitingConditions(mode, class_map(ic.class), node_map(ic.node), ic.count);
                    end
                end
            end
            % Firing outcomes
            if isfield(md, 'firingOutcomes')
                foList = md.firingOutcomes;
                if isstruct(foList), foList = num2cell(foList); end
                for fi = 1:length(foList)
                    fo = foList{fi};
                    if node_map.isKey(fo.node) && class_map.isKey(fo.class)
                        tnode.setFiringOutcome(mode, class_map(fo.class), node_map(fo.node), fo.count);
                    end
                end
            end
            % Marking-dependent firing-rate multiplier (after enabling/timing/
            % distribution are set so the setter guard sees the final mode state).
            if isfield(md, 'firingRateDependence') && ~isempty(md.firingRateDependence)
                g = firingdep_table_to_handle(md.firingRateDependence, node_map, class_map);
                if ~isempty(g)
                    tnode.setFiringRateDependence(mode, g);
                end
            end
        end
    end
end

% --- Restore initial state for Place nodes ---
if isfield(data, 'nodes')
    nds4 = data.nodes;
    if isstruct(nds4), nds4 = num2cell(nds4); end
    for i = 1:length(nds4)
        nd4 = nds4{i};
        if isfield(nd4, 'initialState') && node_map.isKey(nd4.name)
            nodeObj = node_map(nd4.name);
            if isa(nodeObj, 'Place')
                stVal = nd4.initialState;
                stVal = stVal(:)';  % Ensure row vector (jsondecode returns column vectors)
                nodeObj.setState(stVal);
            end
        end
    end
end

% --- Build routing ---
if isfield(data, 'routing') && isfield(data.routing, 'type') && strcmp(data.routing.type, 'matrix')
    P = model.initRoutingMatrix();
    K = length(classesList);
    M = length(nodeList);

    % Build node/class index maps
    nodeIdx = containers.Map();
    for i = 1:M
        nodeIdx(nodeList{i}.name) = i;
    end
    classIdx = containers.Map();
    for i = 1:K
        classIdx(classesList{i}.name) = i;
    end

    % Parse routing keys from raw JSON to preserve commas
    routingEntries = parse_routing_keys(rawJson, class_map, node_map);

    for e = 1:length(routingEntries)
        re = routingEntries{e};
        r = classIdx(re.className1);
        s = classIdx(re.className2);
        ii = nodeIdx(re.fromNode);
        jj = nodeIdx(re.toNode);
        P{r,s}(ii, jj) = re.prob;
    end

    model.link(P);
end

% --- Restore routing strategies ---
if isfield(data, 'routingStrategies')
    stratMap = containers.Map();
    stratMap('RAND') = RoutingStrategy.RAND;
    stratMap('RROBIN') = RoutingStrategy.RROBIN;
    stratMap('WRROBIN') = RoutingStrategy.WRROBIN;
    stratMap('JSQ') = RoutingStrategy.JSQ;
    stratMap('SQ') = RoutingStrategy.SQ;
    stratMap('FIRING') = RoutingStrategy.FIRING;
    stratMap('RL') = RoutingStrategy.RL;
    stratMap('DISABLED') = RoutingStrategy.DISABLED;

    rsFields = fieldnames(data.routingStrategies);
    for fi = 1:length(rsFields)
        nodeName = rsFields{fi};
        if node_map.isKey(nodeName)
            nodeObj = node_map(nodeName);
            classStrats = data.routingStrategies.(nodeName);
            csFields = fieldnames(classStrats);
            for ci = 1:length(csFields)
                className = csFields{ci};
                stratName = classStrats.(className);
                if class_map.isKey(className) && stratMap.isKey(stratName)
                    % Skip RAND, PROB: already handled by routing matrix
                    % Skip WRROBIN: handled separately in routingWeights section
                    rs = stratMap(stratName);
                    if rs ~= RoutingStrategy.RAND && rs ~= RoutingStrategy.PROB && rs ~= RoutingStrategy.WRROBIN
                        nodeObj.setRouting(class_map(className), rs);
                    end
                end
            end
        end
    end
end

% --- Restore routing weights (WRROBIN) ---
if isfield(data, 'routingWeights')
    rwFields = fieldnames(data.routingWeights);
    for fi = 1:length(rwFields)
        nodeName = rwFields{fi};
        if node_map.isKey(nodeName)
            nodeObj = node_map(nodeName);
            classWeights = data.routingWeights.(nodeName);
            cwFields = fieldnames(classWeights);
            for ci = 1:length(cwFields)
                className = cwFields{ci};
                destWeights = classWeights.(className);
                if class_map.isKey(className)
                    % Clear existing routing entries for this class
                    % (link() may have set PROB entries that would accumulate)
                    classIdx = class_map(className).index;
                    if length(nodeObj.output.outputStrategy) >= classIdx && ...
                            length(nodeObj.output.outputStrategy{1, classIdx}) >= 3
                        nodeObj.output.outputStrategy{1, classIdx}{3} = {};
                    end
                    dwFields = fieldnames(destWeights);
                    for di = 1:length(dwFields)
                        destName = dwFields{di};
                        weight = destWeights.(destName);
                        if node_map.isKey(destName)
                            nodeObj.setRouting(class_map(className), RoutingStrategy.WRROBIN, node_map(destName), weight);
                        end
                    end
                end
            end
        end
    end
end

% --- Restore setup / delay-off, polling type and switchover times ---
% jsondecode yields a cell array when node objects carry different field sets; normalized to cells here -- see _kb/12-interfaces-and-docs.md
ndsSo = data.nodes;
if isstruct(ndsSo), ndsSo = num2cell(ndsSo); end
for ni = 1:length(ndsSo)
    nd = ndsSo{ni};
    if ~node_map.isKey(nd.name)
        continue;
    end

    % Setup / delay-off. The writer emits the two maps together, keyed by
    % class name, because setDelayOff requires both distributions.
    if isfield(nd, 'setupTime') && ~isempty(nd.setupTime) && ...
            isfield(nd, 'delayOffTime') && ~isempty(nd.delayOffTime)
        nodeObj = node_map(nd.name);
        suNames = fieldnames(nd.setupTime);
        for fi = 1:length(suNames)
            cname = suNames{fi};
            if ~class_map.isKey(cname) || ~isfield(nd.delayOffTime, cname)
                continue;
            end
            suDist = json2dist(nd.setupTime.(cname));
            doffDist = json2dist(nd.delayOffTime.(cname));
            if ~isempty(suDist) && ~isempty(doffDist)
                nodeObj.setDelayOff(class_map(cname), suDist, doffDist);
            end
        end
    end

    % Polling type restored by name, must precede switchover restore (setPollingType resets it to Immediate)
    if isfield(nd, 'pollingType') && ~isempty(nd.pollingType)
        nodeObj = node_map(nd.name);
        ptId = PollingType.fromName(nd.pollingType);
        if ptId == PollingType.KLIMITED
            if isfield(nd, 'pollingPar') && ~isempty(nd.pollingPar)
                nodeObj.setPollingType(ptId, nd.pollingPar);
            else
                nodeObj.setPollingType(ptId, 1);
            end
        else
            nodeObj.setPollingType(ptId);
        end
    end

    % Switchover times: entries without a "to" field carry the per-class
    % polling form, entries with one the (from,to) pair form.
    if isfield(nd, 'switchoverTimes') && ~isempty(nd.switchoverTimes)
        nodeObj = node_map(nd.name);
        soArr = nd.switchoverTimes;
        if ~iscell(soArr)
            soArr = num2cell(soArr);
        end
        for si = 1:length(soArr)
            so = soArr{si};
            if ~class_map.isKey(so.from)
                continue;
            end
            fromCls = class_map(so.from);
            dist = json2dist(so.distribution);
            if isempty(dist)
                continue;
            end
            if isfield(so, 'to') && ~isempty(so.to)
                if ~class_map.isKey(so.to)
                    continue;
                end
                nodeObj.setSwitchover(fromCls, class_map(so.to), dist);
            else
                nodeObj.setSwitchover(fromCls, dist);
            end
        end
    end
end

% --- Restore finite capacity regions ---
if isfield(data, 'finiteCapacityRegions')
    fcrArr = data.finiteCapacityRegions;
    if ~iscell(fcrArr)
        fcrArr = {fcrArr};
    end
    classes = model.getClasses();
    for ri = 1:length(fcrArr)
        rj = fcrArr{ri};
        regNodes = {};
        % Support both old ("nodes" list) and new ("stations" array) format;
        % jsondecode yields a STRUCT ARRAY here, so multi-station regions need stArr(si) indexing -- see _kb/12-interfaces-and-docs.md
        if isfield(rj, 'stations')
            stArr = rj.stations;
            if iscell(stArr)
                for si = 1:length(stArr)
                    nodeName = stArr{si}.node;
                    if node_map.isKey(nodeName)
                        regNodes{end+1} = node_map(nodeName); %#ok<AGROW>
                    end
                end
            else
                for si = 1:length(stArr)
                    nodeName = stArr(si).node;
                    if node_map.isKey(nodeName)
                        regNodes{end+1} = node_map(nodeName); %#ok<AGROW>
                    end
                end
            end
        elseif isfield(rj, 'nodes')
            nodeNames = rj.nodes;
            if ~iscell(nodeNames), nodeNames = {nodeNames}; end
            for ni = 1:length(nodeNames)
                if node_map.isKey(nodeNames{ni})
                    regNodes{end+1} = node_map(nodeNames{ni}); %#ok<AGROW>
                end
            end
        end
        maxJobs = FiniteCapacityRegion.UNBOUNDED;
        if isfield(rj, 'globalMaxJobs')
            maxJobs = rj.globalMaxJobs;
        end
        if ~isempty(regNodes)
            try
                region = model.addRegion(regNodes);
                if isfield(rj, 'name') && ~isempty(rj.name)
                    region.setName(rj.name);
                end
                if maxJobs ~= FiniteCapacityRegion.UNBOUNDED
                    region.setGlobalMaxJobs(maxJobs);
                end
                % globalMaxMemory
                if isfield(rj, 'globalMaxMemory')
                    region.globalMaxMemory = rj.globalMaxMemory;
                end
                % classMaxJobs
                if isfield(rj, 'classMaxJobs')
                    cmj = rj.classMaxJobs;
                    cmjFields = fieldnames(cmj);
                    for ci = 1:length(cmjFields)
                        cname = cmjFields{ci};
                        if class_map.isKey(cname)
                            jc = class_map(cname);
                            region.classMaxJobs(jc.index) = cmj.(cname);
                        end
                    end
                end
                % dropRule
                if isfield(rj, 'dropRule')
                    drData = rj.dropRule;
                    drFields = fieldnames(drData);
                    for di = 1:length(drFields)
                        cname = drFields{di};
                        if class_map.isKey(cname)
                            jc = class_map(cname);
                            region.dropRule(jc.index) = str_to_droprule(drData.(cname));
                        end
                    end
                end
                % Per-station classWeight and classSize from stations array
                if isfield(rj, 'stations')
                    stArr2 = rj.stations;
                    if ~iscell(stArr2), stArr2 = {stArr2}; end
                    for si = 1:length(stArr2)
                        sj2 = stArr2{si};
                        if isfield(sj2, 'classWeight')
                            cwData = sj2.classWeight;
                            cwFields = fieldnames(cwData);
                            for ci = 1:length(cwFields)
                                cname = cwFields{ci};
                                if class_map.isKey(cname)
                                    jc = class_map(cname);
                                    region.classWeight(jc.index) = cwData.(cname);
                                end
                            end
                        end
                        if isfield(sj2, 'classSize')
                            csData = sj2.classSize;
                            csFields = fieldnames(csData);
                            for ci = 1:length(csFields)
                                cname = csFields{ci};
                                if class_map.isKey(cname)
                                    jc = class_map(cname);
                                    region.classSize(jc.index) = csData.(cname);
                                end
                            end
                        end
                    end
                end
                % Linear constraints A * x <= b
                if isfield(rj, 'constraintA') && isfield(rj, 'constraintB')
                    Adata = rj.constraintA;
                    bdata = rj.constraintB;
                    if iscell(Adata)
                        nrows = length(Adata);
                        K = length(region.classes);
                        A = zeros(nrows, K);
                        for ar = 1:nrows
                            row = Adata{ar};
                            A(ar, 1:length(row)) = row(:)';
                        end
                    else
                        A = Adata;
                    end
                    region.setConstraint(A, bdata(:));
                end
            catch
            end
        end
    end
end

% --- Rewards ---
if isfield(data, 'rewards')
    rewardsArr = data.rewards;
    if isstruct(rewardsArr), rewardsArr = num2cell(rewardsArr); end
    for i = 1:length(rewardsArr)
        rw = rewardsArr{i};
        if ~isfield(rw, 'name') || ~isfield(rw, 'type')
            line_warning(mfilename, 'Ignoring a reward entry without a "name" or "type" field.');
            continue;
        end
        rname = rw.name;
        rtype = rw.type;
        if ~isfield(rw, 'node') || isempty(rw.node) || ~node_map.isKey(rw.node)
            line_warning(mfilename, sprintf(['Reward "%s" refers to node "%s", which is not defined in this ' ...
                'model; the reward is ignored.'], rname, char(getfield_default(rw, 'node', ''))));
            continue;
        end
        rnode = node_map(rw.node);
        rclass = [];
        if isfield(rw, 'class') && ~isempty(rw.class)
            if ~class_map.isKey(rw.class)
                line_warning(mfilename, sprintf(['Reward "%s" refers to class "%s", which is not defined in ' ...
                    'this model; the reward is ignored.'], rname, rw.class));
                continue;
            end
            rclass = class_map(rw.class);
        end
        switch rtype
            case 'QLen'
                if isempty(rclass)
                    model.setReward(rname, Reward.queueLength(rnode));
                else
                    model.setReward(rname, Reward.queueLength(rnode, rclass));
                end
            case 'Util'
                if isempty(rclass)
                    model.setReward(rname, Reward.utilization(rnode));
                else
                    model.setReward(rname, Reward.utilization(rnode, rclass));
                end
            case 'Blocking'
                model.setReward(rname, Reward.blocking(rnode));
            otherwise
                line_warning(mfilename, sprintf(['Reward "%s" has type "%s", for which no reward template is ' ...
                    'implemented; the reward is ignored.'], rname, rtype));
        end
    end
end
end


function v = getfield_default(s, fieldName, defaultValue)
% Return s.(fieldName) when present, otherwise defaultValue.
if isfield(s, fieldName)
    v = s.(fieldName);
else
    v = defaultValue;
end
end


function node = create_node(model, nd, name, ntype)
% Create a node from JSON data.
switch ntype
    case 'Source'
        node = Source(model, name);
    case 'Sink'
        node = Sink(model, name);
    case 'Delay'
        node = Delay(model, name);
    case 'Queue'
        schedStr = 'FCFS';
        if isfield(nd, 'scheduling')
            schedStr = nd.scheduling;
        end
        schedId = str_to_sched_id(schedStr);
        node = Queue(model, name, schedId);
        if isfield(nd, 'servers')
            ns = servers_from_json(nd.servers);
            if isinf(ns) || ns > 1
                node.setNumberOfServers(ns);
            end
        end
        if isfield(nd, 'buffer') && isfinite(nd.buffer)
            node.cap = nd.buffer;
        end
    case 'Fork'
        node = Fork(model, name);
    case 'Join'
        node = Join(model, name);
    case 'Router'
        node = Router(model, name);
    case 'ClassSwitch'
        node = ClassSwitch(model, name);
    case 'Cache'
        % Accept both nested MATLAB (nd.cache.items/capacity/replacement) and flat JAR/Python (nd.numItems/...) schema
        cc = struct();
        if isfield(nd, 'cache')
            cc = nd.cache;
        end
        nitems = 10;
        if isfield(cc, 'items'), nitems = cc.items;
        elseif isfield(nd, 'numItems'), nitems = nd.numItems; end
        cap = 1;
        if isfield(cc, 'capacity'), cap = cc.capacity;
        elseif isfield(nd, 'itemLevelCap'), cap = nd.itemLevelCap(:)'; end
        replStr = 'LRU';
        if isfield(cc, 'replacement'), replStr = cc.replacement;
        elseif isfield(nd, 'replacementStrategy'), replStr = nd.replacementStrategy; end
        replId = str_to_repl_id(replStr);
        node = Cache(model, name, nitems, cap, replId);
    case 'Place'
        % Scheduling strategy supplied to the ctor: installQueueServer (first setService) picks the server section from it
        if isfield(nd, 'scheduling')
            node = Place(model, name, str_to_sched_id(nd.scheduling));
        else
            node = Place(model, name);
        end
        if isfield(nd, 'servers')
            node.numberOfServers = servers_from_json(nd.servers);
        end
        if isfield(nd, 'buffer') && isfinite(nd.buffer)
            node.cap = nd.buffer;
        end
    case 'Transition'
        node = Transition(model, name);
    otherwise
        node = Queue(model, name, SchedStrategy.FCFS);
end
end


function ns = servers_from_json(v)
% Decode a server count. Inf crosses the wire as the string "Infinity" (the form
% Place numServers already used); a numeric value is taken verbatim.
if ischar(v) || isstring(v)
    if strcmpi(v, 'Infinity')
        ns = Inf;
    else
        ns = str2double(v);
    end
else
    ns = double(v);
end
end


% =========================================================================
%  LayeredNetwork deserialization
% =========================================================================

function mult = mult_from_json(v)
% Decode a JSON multiplicity, mapping the infinite-multiplicity sentinel
% (Java's Integer.MAX_VALUE, as written by the JAR and by linemodel_save) back
% to Inf. Negative values are also treated as infinite: earlier versions of the
% Python writer emitted -1 for infinite-server hosts, so files in that format
% must keep loading rather than silently becoming single-server stations.
if isempty(v)
    mult = 1;
elseif v >= 2147483647 || v < 0
    mult = Inf;
else
    mult = v;
end
end

function model = json2layered(data)
% Reconstruct a LayeredNetwork from decoded JSON struct.

modelName = 'model';
if isfield(data, 'name')
    modelName = data.name;
end
model = LayeredNetwork(modelName);

% --- Processors (Python schema: "processors", JAR schema: "hosts") ---
proc_map = containers.Map();
if isfield(data, 'processors')
    procs = data.processors;
elseif isfield(data, 'hosts')
    procs = data.hosts;
else
    procs = [];
end
if ~isempty(procs)
    if isstruct(procs), procs = num2cell(procs); end
    for i = 1:length(procs)
        pd = procs{i};
        pname = pd.name;
        mult = 1;
        if isfield(pd, 'multiplicity'), mult = mult_from_json(pd.multiplicity); end
        schedStr = 'INF';
        if isfield(pd, 'scheduling'), schedStr = pd.scheduling; end
        schedId = str_to_sched_id(schedStr);
        quantum = 0.001;
        if isfield(pd, 'quantum'), quantum = pd.quantum; end
        sf = 1.0;
        if isfield(pd, 'speedFactor'), sf = pd.speedFactor; end
        proc = Host(model, pname, mult, schedId, quantum, sf);
        if isfield(pd, 'replication') && pd.replication > 1
            proc.setReplication(pd.replication);
        end
        proc_map(pname) = proc;
    end
end

% --- Tasks ---
task_map = containers.Map();
if isfield(data, 'tasks')
    tsks = data.tasks;
    if isstruct(tsks), tsks = num2cell(tsks); end
    for i = 1:length(tsks)
        td = tsks{i};
        tname = td.name;
        mult = 1;
        if isfield(td, 'multiplicity'), mult = mult_from_json(td.multiplicity); end
        schedStr = 'INF';
        if isfield(td, 'scheduling'), schedStr = td.scheduling; end
        schedId = str_to_sched_id(schedStr);
        taskType = 'Task';
        if isfield(td, 'taskType'), taskType = td.taskType; end
        if strcmp(taskType, 'FunctionTask')
            task = FunctionTask(model, tname, mult, schedId);
        elseif strcmp(taskType, 'CacheTask')
            totalItems = 1;
            if isfield(td, 'totalItems'), totalItems = td.totalItems; end
            cacheCap = 1;
            if isfield(td, 'cacheCapacity'), cacheCap = td.cacheCapacity; end
            rsStr = 'FIFO';
            if isfield(td, 'replacementStrategy'), rsStr = td.replacementStrategy; end
            rsMap = containers.Map({'RR','FIFO','SFIFO','LRU'}, ...
                {ReplacementStrategy.RR, ReplacementStrategy.FIFO, ...
                 ReplacementStrategy.SFIFO, ReplacementStrategy.LRU});
            if rsMap.isKey(upper(rsStr))
                rs = rsMap(upper(rsStr));
            else
                rs = ReplacementStrategy.FIFO;
            end
            task = CacheTask(model, tname, totalItems, cacheCap, rs, mult, schedId);
        else
            task = Task(model, tname, mult, schedId);
        end
        % Assign to processor (Python schema: "processor", JAR schema: "host")
        procRef = '';
        if isfield(td, 'processor'), procRef = td.processor;
        elseif isfield(td, 'host'), procRef = td.host;
        end
        if ~isempty(procRef) && proc_map.isKey(procRef)
            task.on(proc_map(procRef));
        end
        % Think time (Python schema: "thinkTime" as dist, JAR schema: "thinkTimeMean"/"thinkTimeSCV")
        if isfield(td, 'thinkTime')
            dist = json2dist(td.thinkTime);
            if ~isempty(dist)
                task.setThinkTime(dist);
            end
        elseif isfield(td, 'thinkTimeMean') && td.thinkTimeMean > 0
            task.setThinkTime(Exp(1.0 / td.thinkTimeMean));
        end
        % Setup time
        if isfield(td, 'setupTime')
            dist = json2dist(td.setupTime);
            if ~isempty(dist)
                task.setSetupTime(dist);
            end
        elseif isfield(td, 'setupTimeMean') && td.setupTimeMean > 1e-8
            task.setSetupTime(Exp(1.0 / td.setupTimeMean));
        end
        % Delay-off time
        if isfield(td, 'delayOffTime')
            dist = json2dist(td.delayOffTime);
            if ~isempty(dist)
                task.setDelayOffTime(dist);
            end
        elseif isfield(td, 'delayOffTimeMean') && td.delayOffTimeMean > 1e-8
            task.setDelayOffTime(Exp(1.0 / td.delayOffTimeMean));
        end
        % Fan in
        if isfield(td, 'fanIn') && isstruct(td.fanIn)
            fnames = fieldnames(td.fanIn);
            for fi = 1:length(fnames)
                task.setFanIn(fnames{fi}, td.fanIn.(fnames{fi}));
            end
        end
        % Fan out
        if isfield(td, 'fanOut') && isstruct(td.fanOut)
            fnames = fieldnames(td.fanOut);
            for fi = 1:length(fnames)
                task.setFanOut(fnames{fi}, td.fanOut.(fnames{fi}));
            end
        end
        % Replication
        if isfield(td, 'replication') && td.replication > 1
            task.setReplication(td.replication);
        end
        task_map(tname) = task;
    end
end

% --- Entries ---
entry_map = containers.Map();
if isfield(data, 'entries')
    ents = data.entries;
    if isstruct(ents), ents = num2cell(ents); end
    for i = 1:length(ents)
        ed = ents{i};
        ename = ed.name;
        entryType = 'Entry';
        if isfield(ed, 'entryType'), entryType = ed.entryType; end
        if strcmp(entryType, 'ItemEntry')
            totalItems = 1;
            if isfield(ed, 'totalItems'), totalItems = ed.totalItems; end
            accessProb = [];
            if isfield(ed, 'accessProb')
                ap = ed.accessProb;
                if isstruct(ap)
                    accessProb = json2dist(ap);
                elseif isnumeric(ap)
                    accessProb = DiscreteSampler(ap);
                end
            end
            if isempty(accessProb)
                % Default uniform distribution
                accessProb = DiscreteSampler(ones(1, totalItems) / totalItems);
            end
            entry = ItemEntry(model, ename, totalItems, accessProb);
        else
            entry = Entry(model, ename);
        end
        if isfield(ed, 'task') && task_map.isKey(ed.task)
            entry.on(task_map(ed.task));
        end
        % Entry arrival distribution
        if isfield(ed, 'arrival')
            dist = json2dist(ed.arrival);
            if ~isempty(dist)
                entry.setArrival(dist);
            end
        end
        entry_map(ename) = entry;
    end
end

% --- Activities ---
act_map = containers.Map();
if isfield(data, 'activities')
    acts = data.activities;
    if isstruct(acts), acts = num2cell(acts); end
    for i = 1:length(acts)
        ad = acts{i};
        aname = ad.name;

        % Host demand
        hd = GlobalConstants.FineTol;
        if isfield(ad, 'hostDemand')
            hdDist = json2dist(ad.hostDemand);
            if ~isempty(hdDist)
                hd = hdDist;
            end
        end

        % Bound to entry (Python schema: "boundTo", JAR schema: "boundToEntry")
        bte = '';
        if isfield(ad, 'boundTo')
            bte = ad.boundTo;
        elseif isfield(ad, 'boundToEntry')
            bte = ad.boundToEntry;
        end

        act = Activity(model, aname, hd, bte);

        % Assign to task
        if isfield(ad, 'task') && task_map.isKey(ad.task)
            act.on(task_map(ad.task));
        end

        % Replies to entry
        if isfield(ad, 'repliesTo') && entry_map.isKey(ad.repliesTo)
            act.repliesTo(entry_map(ad.repliesTo));
        end

        % Synch calls (Python schema: "entry", JAR schema: "dest")
        if isfield(ad, 'synchCalls')
            scs = ad.synchCalls;
            if isstruct(scs), scs = num2cell(scs); end
            for j = 1:length(scs)
                sc = scs{j};
                if isfield(sc, 'entry'), ename = sc.entry;
                elseif isfield(sc, 'dest'), ename = sc.dest;
                else, continue;
                end
                meanCalls = 1.0;
                if isfield(sc, 'mean'), meanCalls = sc.mean; end
                if entry_map.isKey(ename)
                    act.synchCall(entry_map(ename), meanCalls);
                end
            end
        end

        % Asynch calls (Python schema: "entry", JAR schema: "dest")
        if isfield(ad, 'asynchCalls')
            acs = ad.asynchCalls;
            if isstruct(acs), acs = num2cell(acs); end
            for j = 1:length(acs)
                ac = acs{j};
                if isfield(ac, 'entry'), ename = ac.entry;
                elseif isfield(ac, 'dest'), ename = ac.dest;
                else, continue;
                end
                meanCalls = 1.0;
                if isfield(ac, 'mean'), meanCalls = ac.mean; end
                if entry_map.isKey(ename)
                    act.asynchCall(entry_map(ename), meanCalls);
                end
            end
        end

        act_map(aname) = act;
    end
end

% --- Precedences (Python schema: "type"/"activities", JAR schema: "preActs"/"postActs"/"preType"/"postType") ---
if isfield(data, 'precedences')
    precs = data.precedences;
    if isstruct(precs), precs = num2cell(precs); end
    for i = 1:length(precs)
        pd = precs{i};
        if ~isfield(pd, 'task') || ~task_map.isKey(pd.task)
            continue;
        end
        task = task_map(pd.task);

        if isfield(pd, 'preActs') || isfield(pd, 'postActs')
            % JAR schema
            preNames = {};
            postNames = {};
            if isfield(pd, 'preActs')
                preNames = pd.preActs;
                if ischar(preNames), preNames = {preNames}; end
            end
            if isfield(pd, 'postActs')
                postNames = pd.postActs;
                if ischar(postNames), postNames = {postNames}; end
            end
            preType = 'pre';
            postType = 'post';
            if isfield(pd, 'preType'), preType = pd.preType; end
            if isfield(pd, 'postType'), postType = pd.postType; end

            % Normalize JAR naming convention to Python convention
            switch postType
                case 'post-AND', postType = 'and-fork';
                case 'post-OR', postType = 'or-fork';
                case 'post-LOOP', postType = 'loop';
            end
            switch preType
                case 'pre-AND', preType = 'and-join';
                case 'pre-OR', preType = 'or-join';
            end

            % Extract postParams (JAR schema: probabilities/loopCount)
            postParams = [];
            if isfield(pd, 'postParams')
                postParams = pd.postParams;
                if iscell(postParams), postParams = cell2mat(postParams); end
            end

            preActs = {};
            for ai = 1:length(preNames)
                if act_map.isKey(preNames{ai})
                    preActs{end+1} = act_map(preNames{ai}); %#ok<AGROW>
                end
            end
            postActs = {};
            for ai = 1:length(postNames)
                if act_map.isKey(postNames{ai})
                    postActs{end+1} = act_map(postNames{ai}); %#ok<AGROW>
                end
            end

            if strcmp(preType, 'pre') && strcmp(postType, 'post')
                if length(preActs) == 1 && length(postActs) == 1
                    ap = ActivityPrecedence.Serial(preActs{1}, postActs{1});
                    task.addPrecedence(ap);
                end
            elseif strcmp(preType, 'pre') && strcmp(postType, 'and-fork')
                if ~isempty(preActs) && ~isempty(postActs)
                    ap = ActivityPrecedence.AndFork(preActs{1}, postActs);
                    task.addPrecedence(ap);
                end
            elseif strcmp(preType, 'and-join') && strcmp(postType, 'post')
                if ~isempty(preActs) && ~isempty(postActs)
                    ap = ActivityPrecedence.AndJoin(preActs, postActs{1});
                    task.addPrecedence(ap);
                end
            elseif strcmp(preType, 'pre') && strcmp(postType, 'or-fork')
                if ~isempty(preActs) && ~isempty(postActs)
                    probs = [];
                    if isfield(pd, 'probabilities')
                        probs = pd.probabilities;
                        if isstruct(probs), probs = cell2mat(struct2cell(probs)); end
                    end
                    if isempty(probs) && ~isempty(postParams)
                        probs = postParams(:)';
                    end
                    if isempty(probs)
                        n = length(postActs);
                        probs = ones(1, n) / n;
                    end
                    ap = ActivityPrecedence.OrFork(preActs{1}, postActs, probs);
                    task.addPrecedence(ap);
                end
            elseif strcmp(preType, 'or-join') && strcmp(postType, 'post')
                if ~isempty(preActs) && ~isempty(postActs)
                    ap = ActivityPrecedence.OrJoin(preActs, postActs{1});
                    task.addPrecedence(ap);
                end
            elseif strcmp(preType, 'pre') && strcmp(postType, 'loop')
                count = 1.0;
                if isfield(pd, 'loopCount'), count = pd.loopCount; end
                if count == 1.0 && ~isempty(postParams)
                    count = postParams(1);
                end
                if ~isempty(preActs) && ~isempty(postActs)
                    if length(postActs) > 1
                        ap = ActivityPrecedence.Loop(preActs{1}, postActs(1:end-1), postActs{end}, count);
                    else
                        ap = ActivityPrecedence.Loop(preActs{1}, postActs, count);
                    end
                    task.addPrecedence(ap);
                end
            elseif strcmp(preType, 'pre') && strcmp(postType, 'post-CACHE')
                if ~isempty(preActs) && ~isempty(postActs)
                    ap = ActivityPrecedence.CacheAccess(preActs{1}, postActs);
                    task.addPrecedence(ap);
                end
            end
        else
            % Python schema
            ptype = pd.type;
            actNames = pd.activities;
            if ischar(actNames), actNames = {actNames}; end

            % Resolve activity names to objects
            actObjs = {};
            for ai = 1:length(actNames)
                an = actNames{ai};
                if act_map.isKey(an)
                    actObjs{end+1} = act_map(an); %#ok<AGROW>
                end
            end
            % A Loop's trigger is separate ('preActivity'), so its body may be a single activity -- see _kb/04-networkstruct.md
            isLoopWithPre = strcmp(ptype, 'Loop') && isfield(pd, 'preActivity') ...
                && act_map.isKey(pd.preActivity) && ~isempty(actObjs);
            if length(actObjs) < 2 && ~isLoopWithPre
                continue;
            end

            switch ptype
                case 'Serial'
                    ap = ActivityPrecedence.Serial(actObjs{:});
                    task.addPrecedence(ap);
                case 'AndFork'
                    ap = ActivityPrecedence.AndFork(actObjs{1}, actObjs(2:end));
                    task.addPrecedence(ap);
                case 'AndJoin'
                    ap = ActivityPrecedence.AndJoin(actObjs(1:end-1), actObjs{end});
                    task.addPrecedence(ap);
                case 'OrFork'
                    probs = [];
                    if isfield(pd, 'probabilities')
                        probs = pd.probabilities;
                        if isstruct(probs), probs = cell2mat(struct2cell(probs)); end
                    end
                    if isempty(probs)
                        n = length(actObjs) - 1;
                        probs = ones(1, n) / n;
                    end
                    ap = ActivityPrecedence.OrFork(actObjs{1}, actObjs(2:end), probs);
                    task.addPrecedence(ap);
                case 'OrJoin'
                    ap = ActivityPrecedence.OrJoin(actObjs(1:end-1), actObjs{end});
                    task.addPrecedence(ap);
                case 'Loop'
                    count = 1.0;
                    if isfield(pd, 'loopCount'), count = pd.loopCount; end
                    % Check for explicit preActivity field (new format)
                    if isfield(pd, 'preActivity') && act_map.isKey(pd.preActivity)
                        preAct = act_map(pd.preActivity);
                        ap = ActivityPrecedence.Loop(preAct, actObjs, count);
                    elseif length(actObjs) >= 3
                        % Legacy format: first is pre, rest is body+end
                        ap = ActivityPrecedence.Loop(actObjs{1}, actObjs(2:end-1), actObjs{end}, count);
                    else
                        ap = ActivityPrecedence.Loop(actObjs{1}, actObjs(2:end), count);
                    end
                    task.addPrecedence(ap);
                case 'CacheAccess'
                    if length(actObjs) >= 2
                        ap = ActivityPrecedence.CacheAccess(actObjs{1}, actObjs(2:end));
                        task.addPrecedence(ap);
                    end
            end
        end
    end
end
end


% =========================================================================
%  Workflow deserialization
% =========================================================================

function model = json2workflow(data)
% Reconstruct a Workflow from decoded JSON struct.

modelName = 'workflow';
if isfield(data, 'name')
    modelName = data.name;
end
model = Workflow(modelName);

% --- Activities ---
if isfield(data, 'activities')
    acts = data.activities;
    if isstruct(acts), acts = num2cell(acts); end
    for i = 1:length(acts)
        ad = acts{i};
        actName = ad.name;
        if isfield(ad, 'hostDemand') && ~isempty(ad.hostDemand)
            dist = json2dist(ad.hostDemand);
            model.addActivity(actName, dist);
        else
            model.addActivity(actName, 1.0);
        end
    end
end

% --- Precedences ---
if isfield(data, 'precedences')
    precs = data.precedences;
    if isstruct(precs), precs = num2cell(precs); end
    for i = 1:length(precs)
        pd = precs{i};

        % preActs
        if isfield(pd, 'preActs')
            preActs = cellify_string_array(pd.preActs);
        else
            preActs = {};
        end

        % postActs
        if isfield(pd, 'postActs')
            postActs = cellify_string_array(pd.postActs);
        else
            postActs = {};
        end

        % preType / postType - convert JAR strings to numeric IDs
        preType = ActivityPrecedenceType.PRE_SEQ;
        if isfield(pd, 'preType')
            preType = str_to_prectype(pd.preType);
        end
        postType = ActivityPrecedenceType.POST_SEQ;
        if isfield(pd, 'postType')
            postType = str_to_prectype(pd.postType);
        end

        % preParams / postParams
        preParams = [];
        if isfield(pd, 'preParams') && ~isempty(pd.preParams)
            preParams = pd.preParams(:)';
        end
        postParams = [];
        if isfield(pd, 'postParams') && ~isempty(pd.postParams)
            postParams = pd.postParams(:)';
        end

        ap = ActivityPrecedence(preActs, postActs, preType, postType, preParams, postParams);
        model.addPrecedence(ap);
    end
end
end


% =========================================================================
%  Environment deserialization
% =========================================================================

function model = json2environment(data, rawJson)
% Reconstruct an Environment from decoded JSON struct.

modelName = 'env';
if isfield(data, 'name')
    modelName = data.name;
end
numStages = 0;
if isfield(data, 'numStages')
    numStages = data.numStages;
end
model = Environment(modelName, numStages);

% --- Node failures ---
% "nodeFailures" has two roles depending on whether DOWN_<node> stages are already declared -- see _kb/09-ldes-and-cache.md
nfArr = {};
if isfield(data, 'nodeFailures')
    nfArr = data.nodeFailures;
    if isstruct(nfArr), nfArr = num2cell(nfArr); end
end

stages = {};
if isfield(data, 'stages')
    stages = data.stages;
    if isstruct(stages), stages = num2cell(stages); end
end

declaredNames = cell(1, length(stages));
for i = 1:length(stages)
    if isfield(stages{i}, 'name')
        declaredNames{i} = stages{i}.name;
    else
        declaredNames{i} = sprintf('Stage%d', i);
    end
end

% Expand the macro form only when no DOWN stage is declared for any entry.
macroMode = ~isempty(nfArr);
for k = 1:length(nfArr)
    if ~isfield(nfArr{k}, 'node')
        line_error(mfilename, 'A "nodeFailures" entry is missing the required "node" field.');
    end
    if any(strcmp(declaredNames, sprintf('DOWN_%s', nfArr{k}.node)))
        macroMode = false;
        break;
    end
end

if macroMode
    if length(stages) ~= 1
        line_error(mfilename, ['"nodeFailures" expands the base model into the UP and DOWN_<node> stages, ' ...
            'so "stages" must declare exactly one stage, holding the base (UP) model.']);
    end
    if isfield(data, 'transitions') && ~isempty(data.transitions)
        line_error(mfilename, ['"nodeFailures" implies the breakdown and repair transitions; "transitions" ' ...
            'must not be declared alongside it.']);
    end
    if ~isfield(stages{1}, 'model') || isempty(stages{1}.model)
        line_error(mfilename, '"nodeFailures" requires the base stage to carry a "model".');
    end
    baseModel = json2network(stages{1}.model, rawJson);
    for k = 1:length(nfArr)
        nf = nfArr{k};
        [breakdownDist, repairDist, downServiceDist, resetB, resetR] = nodefailure_fields(nf);
        if isempty(repairDist)
            model.addNodeBreakdown(baseModel, nf.node, breakdownDist, downServiceDist, resetB);
        else
            model.addNodeFailureRepair(baseModel, nf.node, breakdownDist, repairDist, downServiceDist, ...
                resetB, resetR);
        end
    end
else
    % --- Stages ---
    stageNames = {};
    for i = 1:length(stages)
        sd = stages{i};
        stageName = declaredNames{i};
        stageNames{end+1} = stageName; %#ok<AGROW>

        stageType = '';
        if isfield(sd, 'type')
            stageType = sd.type;
        end

        stageModel = [];
        if isfield(sd, 'model') && ~isempty(sd.model)
            stageModel = json2network(sd.model, rawJson);
        end

        if ~isempty(stageModel)
            model.addStage(stageName, stageType, stageModel);
        end
    end

    % --- Transitions ---
    if isfield(data, 'transitions')
        trans = data.transitions;
        if isstruct(trans), trans = num2cell(trans); end
        for i = 1:length(trans)
            td = trans{i};
            fromIdx = td.from + 1;  % Convert from 0-indexed (JAR) to 1-indexed (MATLAB)
            toIdx = td.to + 1;      % Convert from 0-indexed (JAR) to 1-indexed (MATLAB)
            if isfield(td, 'distribution') && ~isempty(td.distribution)
                dist = json2dist(td.distribution);
                if ~isempty(dist) && ~isa(dist, 'Disabled')
                    % Use stage names for MATLAB Environment API
                    if fromIdx <= length(stageNames) && toIdx <= length(stageNames)
                        model.addTransition(stageNames{fromIdx}, stageNames{toIdx}, dist);
                    end
                end
            end
        end
    end

    % Re-attach the node-failure descriptors and their reset policies to the
    % stages just built, so that the environment serializes back identically.
    for k = 1:length(nfArr)
        nf = nfArr{k};
        [breakdownDist, repairDist, downServiceDist, resetB, resetR] = nodefailure_fields(nf);
        model.registerNodeFailure(nf.node, breakdownDist, repairDist, downServiceDist, resetB, resetR);
    end
end

model.init();
end


function [breakdownDist, repairDist, downServiceDist, resetB, resetR] = nodefailure_fields(nf)
% Decode the distributions and reset policies of one "nodeFailures" entry.
% Note that breakdownRate/repairRate carry full distributions, not scalar rates.
if ~isfield(nf, 'breakdownRate') || isempty(nf.breakdownRate)
    line_error(mfilename, sprintf('Node failure on "%s" is missing the required "breakdownRate" field.', nf.node));
end
if ~isfield(nf, 'downService') || isempty(nf.downService)
    line_error(mfilename, sprintf('Node failure on "%s" is missing the required "downService" field.', nf.node));
end
breakdownDist = json2dist(nf.breakdownRate);
downServiceDist = json2dist(nf.downService);
repairDist = [];
if isfield(nf, 'repairRate') && ~isempty(nf.repairRate)
    repairDist = json2dist(nf.repairRate);
end
resetB = 'keep';
if isfield(nf, 'breakdownResetPolicy') && ~isempty(nf.breakdownResetPolicy)
    resetB = nf.breakdownResetPolicy;
end
resetR = 'keep';
if isfield(nf, 'repairResetPolicy') && ~isempty(nf.repairResetPolicy)
    resetR = nf.repairResetPolicy;
end
end


% =========================================================================
%  Distribution deserialization
% =========================================================================

function dist = json2dist(d)
% Convert a JSON distribution struct to a MATLAB Distribution object.
if isempty(d)
    dist = [];
    return;
end

dtype = d.type;

switch dtype
    case 'Disabled'
        dist = Disabled.getInstance();
        return;
    case 'Immediate'
        dist = Immediate.getInstance();
        return;
    case 'Expolynomial'
        % Nested object (the density is a Sirio expression string); "Inf" carries
        % an unbounded latest firing time.
        ep = d.expolynomial;
        if ischar(ep.lft) || isstring(ep.lft)
            lft = Inf;
        else
            lft = double(ep.lft);
        end
        dist = Expolynomial(ep.density, double(ep.eft), lft);
        return;
end

% Direct params
if isfield(d, 'params') && ~isempty(d.params)
    p = d.params;
    switch dtype
        case 'Exp'
            % 'lambda' is canonical; 'rate' is accepted as a read-side alias
            % (the manual's published LQN example uses it, and Python honours it).
            if isfield(p, 'lambda')
                lam = p.lambda;
            elseif isfield(p, 'rate')
                lam = p.rate;
            else
                line_error(mfilename, 'Exp distribution has neither "lambda" nor "rate".');
            end
            dist = Exp(lam);
            return;
        case 'Det'
            dist = Det(p.value);
            return;
        case 'Erlang'
            dist = Erlang(p.lambda, p.k);
            return;
        case 'HyperExp'
            % Two wire forms: 2-phase (scalar p) and n-phase (vector p,lambda) -- see _kb/09-ldes-and-cache.md
            pv = p.p(:)';
            lv = p.lambda(:)';
            if isscalar(pv)
                dist = HyperExp(pv, lv(1), lv(2));
            elseif numel(pv) == 2 && numel(lv) == 2
                dist = HyperExp(pv(1), lv(1), lv(2));
            else
                if numel(pv) ~= numel(lv)
                    line_error(mfilename, sprintf(...
                        'HyperExp has %d probabilities but %d rates.', numel(pv), numel(lv)));
                end
                dist = HyperExp(pv, lv);
            end
            return;
        case 'Gamma'
            dist = Gamma(p.alpha, p.beta);
            return;
        case 'Lognormal'
            dist = Lognormal(p.mu, p.sigma);
            return;
        case 'Uniform'
            dist = Uniform(p.a, p.b);
            return;
        case 'Zipf'
            dist = Zipf(p.s, p.n);
            return;
        case 'Pareto'
            dist = Pareto(p.alpha, p.scale);
            return;
        case 'Weibull'
            % JSON: alpha = scale (getParam(1)), beta = shape (getParam(2))
            % Constructor: Weibull(shape, scale)
            dist = Weibull(p.beta, p.alpha);
            return;
        case 'Normal'
            dist = Normal(p.mu, p.sigma);
            return;
        case 'Geometric'
            dist = Geometric(p.p);
            return;
        case 'Binomial'
            dist = Binomial(p.n, p.p);
            return;
        case 'Poisson'
            dist = Poisson(p.lambda);
            return;
        case 'Bernoulli'
            dist = Bernoulli(p.p);
            return;
        case 'DiscreteUniform'
            dist = DiscreteUniform(p.min, p.max);
            return;
        case 'DiscreteSampler'
            pv = p.p(:)';
            xv = p.x(:)';
            dist = DiscreteSampler(pv, xv);
            return;
        case 'Coxian'
            mu = p.mu(:)';
            phi = p.phi(:)';
            dist = Coxian(mu, phi);
            return;
        case 'MMPP2'
            dist = MMPP2(p.lambda0, p.lambda1, p.sigma0, p.sigma1);
            return;
        case 'MMDP2'
            dist = MMDP2(p.r0, p.r1, p.sigma0, p.sigma1);
            return;
        case 'NHPP'
            % Absent 'cyclic' means cyclic, matching the constructor default.
            if isfield(p, 'cyclic')
                cyc = logical(p.cyclic);
            else
                cyc = true;
            end
            dist = NHPP(p.breakpoints(:)', p.rates(:)', cyc);
            return;
        case 'ME'
            dist = ME(p.alpha(:)', json2mat(p.A));
            return;
        case 'RAP'
            dist = RAP(json2mat(p.H0), json2mat(p.H1));
            return;
        case 'DMAP'
            dist = DMAP(json2mat(p.D0), json2mat(p.D1));
            return;
        case 'BMAP'
            % D = {D0, D1, ..., Dk}, Dk driving batches of size k
            dist = BMAP(json2matcell(p.D));
            return;
        case 'MarkedMMPP'
            % D = {D0, D11, ..., D1K}; the ctor rebuilds the aggregate D1 from
            % the K == length(D)-1 form
            Dcell = json2matcell(p.D);
            if isfield(p, 'K')
                Kmarks = p.K;
            else
                Kmarks = numel(Dcell) - 1;
            end
            dist = MarkedMMPP(Dcell, Kmarks);
            return;
        case 'EmpiricalCDF'
            dist = EmpiricalCDF(p.x(:), p.F(:));
            return;
        case 'Replayer'
            % Try to load from file first
            if isfield(p, 'fileName') && exist(p.fileName, 'file') == 2
                dist = Replayer(p.fileName);
                return;
            end
            % Fallback to APH if available
            if isfield(d, 'ph') && ~isempty(d.ph)
                ph = d.ph;
                alpha = ph.alpha;
                T = ph.T;
                if ~isvector(alpha), alpha = alpha(:)'; end
                dist = PH(alpha, T);
                return;
            end
            % Fallback to Exp with stored mean
            m = 1.0;
            if isfield(p, 'mean'), m = p.mean; end
            dist = Exp(1.0 / m);
            return;
    end
end

% Prior distribution (mixture of alternatives with prior probabilities)
if strcmp(dtype, 'Prior')
    if isfield(d, 'distributions') && isfield(d, 'probabilities')
        altJsons = d.distributions;
        probs = d.probabilities;
        if ~iscell(altJsons)
            % jsondecode may return struct array instead of cell
            altJsons = num2cell(altJsons);
        end
        alts = cell(1, length(altJsons));
        for ai = 1:length(altJsons)
            alts{ai} = json2dist(altJsons{ai});
        end
        probs = probs(:)';
        dist = Prior(alts, probs);
        return;
    end
end

% PH/APH representation
if isfield(d, 'ph') && ~isempty(d.ph)
    ph = d.ph;
    alpha = ph.alpha;
    T = ph.T;
    alpha = alpha(:)';  % Ensure row vector (jsondecode returns column vectors)
    if strcmp(dtype, 'APH')
        dist = APH(alpha, T);
    else
        dist = PH(alpha, T);
    end
    return;
end

% MAP representation
if isfield(d, 'map') && ~isempty(d.map)
    mapSpec = d.map;
    D0 = mapSpec.D0;
    D1 = mapSpec.D1;
    dist = MAP(D0, D1);
    return;
end

% Marked MAP representation: {D0, per-mark D1k}; the aggregate D1 is
% rebuilt by the MarkedMAP constructor (K == length(D)-1 form)
if isfield(d, 'mmap') && ~isempty(d.mmap)
    mmapSpec = d.mmap;
    D0 = mmapSpec.D0;
    d1k = mmapSpec.D1k;
    if isnumeric(d1k)
        % jsondecode collapses equally-sized matrices into an ND array:
        % (K x n x n) with marks along the first dimension
        Kmarks = size(d1k, 1);
        Dcell = cell(1, 1 + Kmarks);
        Dcell{1} = D0;
        for km = 1:Kmarks
            Dcell{1+km} = squeeze(d1k(km, :, :));
        end
    else
        if ~iscell(d1k), d1k = num2cell(d1k); end
        Dcell = [{D0}, d1k(:)'];
    end
    dist = MarkedMAP(Dcell, numel(Dcell)-1);
    return;
end

% Fit specification
if isfield(d, 'fit') && ~isempty(d.fit)
    fit = d.fit;
    method = fit.method;
    switch method
        case 'fitMean'
            m = fit.mean;
            switch dtype
                case 'Exp'
                    dist = Exp(1.0 / m);
                case 'Det'
                    dist = Det(m);
                otherwise
                    dist = Exp(1.0 / m);
            end
            return;
        case 'fitMeanAndSCV'
            m = fit.mean;
            scv = fit.scv;
            switch dtype
                case 'Erlang'
                    dist = Erlang.fitMeanAndSCV(m, scv);
                case 'HyperExp'
                    dist = HyperExp.fitMeanAndSCV(m, scv);
                otherwise
                    dist = Exp(1.0 / m);
            end
            return;
        case 'fitMeanAndOrder'
            m = fit.mean;
            order = fit.order;
            switch dtype
                case 'Erlang'
                    dist = Erlang.fitMeanAndOrder(m, order);
                otherwise
                    dist = Exp(1.0 / m);
            end
            return;
    end
end

% Unrecognized type: rebuild an APH matching the writer's fallback mean/SCV, warn -- see _kb/09-ldes-and-cache.md
if isfield(d, 'params') && isfield(d.params, 'mean') && isfield(d.params, 'scv')
    line_warning(mfilename, sprintf(['Distribution type "%s" is not supported on load; ' ...
        'reconstructing an APH fitted to its mean and SCV.\n'], dtype));
    dist = APH.fitMeanAndSCV(d.params.mean, d.params.scv);
    return;
end
if isfield(d, 'params') && isfield(d.params, 'mean')
    line_warning(mfilename, sprintf(['Distribution type "%s" is not supported on load and ' ...
        'carries no SCV; reconstructing an Exp with its mean.\n'], dtype));
    dist = Exp(1.0 / d.params.mean);
    return;
end
line_error(mfilename, sprintf(['Distribution type "%s" is not supported on load and carries ' ...
    'no moments to fit.'], dtype));
end


% =========================================================================
%  Routing parser (handles comma keys in JSON)
% =========================================================================

function entries = parse_routing_keys(rawJson, class_map, node_map)
% Parse routing matrix from raw JSON text to handle keys with commas.
% Returns a cell array of structs with fields:
%   className1, className2, fromNode, toNode, prob
entries = {};

% Build reverse mapping: jsondecode-sanitized name -> original node name
% jsondecode uses matlab.lang.makeValidName which replaces spaces etc.
nodeNames = node_map.keys();
sanitized_map = containers.Map();
for ni = 1:length(nodeNames)
    origName = nodeNames{ni};
    sanitized = matlab.lang.makeValidName(origName);
    sanitized_map(sanitized) = origName;
end

classNames = class_map.keys();

% For each pair of class names, try to find the corresponding key in the JSON
for ri = 1:length(classNames)
    for si = 1:length(classNames)
        cn1 = classNames{ri};
        cn2 = classNames{si};
        keyStr = ['"', cn1, ',', cn2, '"'];

        % Find this key in the raw JSON
        pos = strfind(rawJson, keyStr);
        if isempty(pos)
            continue;
        end

        % For each occurrence, extract the nested from -> to -> prob structure
        for pidx = 1:length(pos)
            startPos = pos(pidx) + length(keyStr);
            % Skip whitespace and colon
            idx = startPos;
            while idx <= length(rawJson) && (rawJson(idx) == ' ' || rawJson(idx) == ':' || rawJson(idx) == newline || rawJson(idx) == char(13) || rawJson(idx) == char(9))
                idx = idx + 1;
            end
            if idx > length(rawJson) || rawJson(idx) ~= '{'
                continue;
            end
            % Extract the JSON object using brace counting
            objStr = extract_json_object(rawJson, idx);
            if isempty(objStr)
                continue;
            end
            % Parse the from -> to -> prob structure
            try
                fromTo = jsondecode(objStr);
                fromNames = fieldnames(fromTo);
                for fi = 1:length(fromNames)
                    fromField = fromNames{fi};
                    toStruct = fromTo.(fromField);
                    toNames = fieldnames(toStruct);
                    % Resolve sanitized field names back to original node names
                    if sanitized_map.isKey(fromField)
                        fromName = sanitized_map(fromField);
                    else
                        fromName = fromField;
                    end
                    for ti = 1:length(toNames)
                        toField = toNames{ti};
                        prob = toStruct.(toField);
                        if sanitized_map.isKey(toField)
                            toName = sanitized_map(toField);
                        else
                            toName = toField;
                        end
                        % Verify names exist in the model
                        if node_map.isKey(fromName) && node_map.isKey(toName)
                            re = struct();
                            re.className1 = cn1;
                            re.className2 = cn2;
                            re.fromNode = fromName;
                            re.toNode = toName;
                            re.prob = prob;
                            entries{end+1} = re; %#ok<AGROW>
                        end
                    end
                end
            catch
                % Skip if parsing fails
            end
        end
    end
end
end


function objStr = extract_json_object(str, startIdx)
% Extract a JSON object string starting at startIdx (must be '{').
if str(startIdx) ~= '{'
    objStr = '';
    return;
end
depth = 0;
inString = false;
escaped = false;
for i = startIdx:length(str)
    c = str(i);
    if escaped
        escaped = false;
        continue;
    end
    if c == '\'
        escaped = true;
        continue;
    end
    if c == '"'
        inString = ~inString;
        continue;
    end
    if ~inString
        if c == '{'
            depth = depth + 1;
        elseif c == '}'
            depth = depth - 1;
            if depth == 0
                objStr = str(startIdx:i);
                return;
            end
        end
    end
end
objStr = '';
end


% =========================================================================
%  Helper functions
% =========================================================================

function id = str_to_sched_id(str)
% Map a wire scheduling enum name to a SchedStrategy numeric ID.
%
% SchedStrategy.fromText case-folds and resolves the FCFSPRIO/HOL, LAS/FB and
% SET/SETF aliases, and errors on an unknown name. Do not reintroduce a
% hand-rolled whitelist here: the previous one covered 23 of 40 strategies,
% silently degraded the rest to FCFS, and lacked the PAS/OI cases that
% linemodel_save itself emits.
if iscell(str)
    str = str{1};
end
switch upper(str)
    case 'RAND'
        % Legacy alias kept for files written before SIRO was named: JMT calls
        % the same discipline "RAND". Not a SchedStrategy.fromText case.
        id = SchedStrategy.SIRO;
    otherwise
        id = SchedStrategy.fromText(lower(char(str)));
end
end


function id = str_to_depdisc(str)
% Map a departure discipline name to a DepartureDiscipline numeric ID. Matched
% case-insensitively, as the JAR reader does.
switch lower(char(str))
    case 'normal', id = DepartureDiscipline.NORMAL;
    case 'fifo',   id = DepartureDiscipline.FIFO;
    otherwise
        line_error(mfilename, sprintf('Unrecognized departure discipline "%s".', str));
end
end


function id = str_to_impatience(str)
% Map an impatience type name to an ImpatienceType numeric ID.
switch lower(char(str))
    case 'reneging', id = ImpatienceType.RENEGING;
    case 'balking',  id = ImpatienceType.BALKING;
    case 'retrial',  id = ImpatienceType.RETRIAL;
    otherwise
        line_error(mfilename, sprintf('Unrecognized impatience type "%s".', str));
end
end


function id = str_to_repl_id(str)
% Map replacement strategy string to ReplacementStrategy numeric ID.
switch upper(str)
    case 'LRU',   id = ReplacementStrategy.LRU;
    case 'FIFO',  id = ReplacementStrategy.FIFO;
    case 'RR',    id = ReplacementStrategy.RR;
    case 'SFIFO', id = ReplacementStrategy.SFIFO;
    case 'HLRU',  id = ReplacementStrategy.HLRU;
    case 'CLIMB', id = ReplacementStrategy.CLIMB;
    case 'QLRU',  id = ReplacementStrategy.QLRU;
    otherwise
        line_error(mfilename, sprintf('Unrecognized replacement strategy "%s".', str));
end
end


function id = str_to_prectype(str)
% Map JAR precedence type string to MATLAB ActivityPrecedenceType numeric ID.
switch str
    case 'pre',        id = ActivityPrecedenceType.PRE_SEQ;
    case 'pre-AND',    id = ActivityPrecedenceType.PRE_AND;
    case 'pre-OR',     id = ActivityPrecedenceType.PRE_OR;
    case 'post',       id = ActivityPrecedenceType.POST_SEQ;
    case 'post-AND',   id = ActivityPrecedenceType.POST_AND;
    case 'post-OR',    id = ActivityPrecedenceType.POST_OR;
    case 'post-LOOP',  id = ActivityPrecedenceType.POST_LOOP;
    case 'post-CACHE', id = ActivityPrecedenceType.POST_CACHE;
    otherwise,         id = ActivityPrecedenceType.PRE_SEQ;
end
end


function id = str_to_droprule(str)
% Map drop rule string to DropStrategy numeric ID.
switch str
    case 'drop',                  id = DropStrategy.DROP;
    case 'waitingQueue',          id = DropStrategy.WAITQ;
    case 'blockingAfterService',  id = DropStrategy.BAS;
    case 'retrial',               id = DropStrategy.RETRIAL;
    case 'retrialWithLimit',      id = DropStrategy.RETRIAL_WITH_LIMIT;
    otherwise,                    id = DropStrategy.WAITQ;
end
end


function m = json2mat(arr)
% Convert a JSON 2D array to a numeric matrix. jsondecode returns a numeric
% matrix for a rectangular array of numbers, but a cell array of row vectors
% when the rows differ in length or when the array was written through a cell
% wrapper.
if iscell(arr)
    m = cell2mat(cellfun(@(r) double(r(:)'), arr(:), 'UniformOutput', false));
else
    m = double(arr);
end
end

function c = json2matcell(arr)
% Convert a JSON array of 2D matrices (which jsondecode returns as a cell
% array of matrices, or collapses into a 3D numeric array when all matrices
% have equal size) into a cell array of 2D matrices.
if iscell(arr)
    c = cell(1, numel(arr));
    for ci = 1:numel(arr)
        m = arr{ci};
        if iscell(m) % array of row arrays (ragged rows)
            c{ci} = cell2mat(cellfun(@(r) r(:)', m(:), 'UniformOutput', false));
        else
            c{ci} = m;
        end
    end
elseif ndims(arr) == 3
    c = cell(1, size(arr, 1));
    for ci = 1:size(arr, 1)
        c{ci} = squeeze(arr(ci, :, :));
    end
elseif ismatrix(arr)
    c = {arr}; % single matrix
else
    c = {};
end
end

function c = cellify_string_array(arr)
% Convert a JSON string array (which may be decoded as a char, cell, or
% struct array) into a cell array of character vectors.
if ischar(arr)
    c = {arr};
elseif isstring(arr)
    c = cellstr(arr);
elseif iscell(arr)
    c = arr;
else
    % jsondecode can return a struct array or char matrix for string arrays
    c = cellstr(arr);
end
end

function v = oi_cutoffs_vec(c)
% Decode the oiCutoffs array. The writer wraps it in a cell so that a
% single-class model still encodes as a JSON array rather than a bare scalar.
if iscell(c)
    v = cell2mat(c(:)');
else
    v = double(c(:)');
end
end

function muFun = oi_table_to_handle(tblStruct, cutoffs)
% Rebuild an OI/PAS total-service-rate handle mu(c) from the materialized
% macrostate table written by OI_RATE_TABLE (linemodel_save). The table is keyed
% by the per-class counts, which is lossless because mu is order-independent;
% the handle therefore reduces the ordered microstate vector c (a list of class
% indices) to its class counts before looking up. Counts are clamped to CUTOFFS,
% so mu saturates beyond the tabulated range exactly as the table intends. As in
% CD_TABLE_TO_HANDLE, jsondecode mangles the JSON keys ("1,1") into valid MATLAB
% identifiers ("x1_1"), so the counts are parsed back out of the field names.
map = containers.Map('KeyType', 'char', 'ValueType', 'double');
fn = fieldnames(tblStruct);
K = numel(cutoffs);
for i = 1:numel(fn)
    nm = fn{i};
    parts = strsplit(nm(2:end), '_');   % drop the 'x' prefix jsondecode prepends
    n = cellfun(@str2double, parts);
    if K == 0
        K = numel(n);
    end
    map(cd_state_key(n)) = double(tblStruct.(nm));
end
muFun = @(c) oi_table_eval(c, map, cutoffs, K);
end

function rate = oi_table_eval(c, map, cutoffs, K)
c = round(double(c(:)'));
cnt = zeros(1, K);
for i = 1:numel(c)
    ci = c(i);
    if ci >= 1 && ci <= K
        cnt(ci) = cnt(ci) + 1;
    end
end
if ~isempty(cutoffs)
    m = min(numel(cnt), numel(cutoffs));
    cnt(1:m) = min(cnt(1:m), cutoffs(1:m));
end
if sum(cnt) == 0
    rate = 0;   % the empty state is omitted from the table: an idle queue
    return;
end
k = cd_state_key(cnt);
if isKey(map, k)
    rate = map(k);
else
    rate = 0;
end
end

function beta = cd_table_to_handle(tblStruct, cutoffs)
% Rebuild a class-dependence handle beta(n) from the materialized lattice table
% written by CD_SCALING_TABLE (linemodel_save). jsondecode mangles the JSON keys
% ("1,1") into valid MATLAB identifiers ("x1_1"), so the per-class counts are
% parsed back out of the field names rather than reconstructed from them. The
% population is clamped to CUTOFFS, so beta saturates beyond the tabulated range
% exactly as the table intends.
map = containers.Map('KeyType', 'char', 'ValueType', 'any');
fn = fieldnames(tblStruct);
for i = 1:numel(fn)
    nm = fn{i};
    parts = strsplit(nm(2:end), '_');   % drop the 'x' prefix jsondecode prepends
    n = cellfun(@str2double, parts);
    v = double(tblStruct.(nm));
    map(cd_state_key(n)) = v(:)';
end
beta = @(ni) cd_table_eval(ni, map, cutoffs);
end

function k = cd_state_key(n)
k = strjoin(arrayfun(@(x) sprintf('%d', x), n(:)', 'UniformOutput', false), ',');
end

function g = firingdep_table_to_handle(frm, node_map, class_map)
% Rebuild the firing-rate dependence handle g(M) from the materialized lattice
% written by FIRINGMOD_SCALING_TABLE (linemodel_save). M is the node-indexed
% marking matrix; the enabling (place,class) slots are read off, clamped to the
% cutoffs, and the tabulated multiplier is looked up (default 1 outside range).
g = [];
if ~isfield(frm, 'slots') || ~isfield(frm, 'scaling'), return; end
slotsData = frm.slots;
if isstruct(slotsData), slotsData = num2cell(slotsData); end
slotIdx = zeros(numel(slotsData), 2);
for s = 1:numel(slotsData)
    sm = slotsData{s};
    if ~node_map.isKey(sm.node) || ~class_map.isKey(sm.class), return; end
    slotIdx(s,:) = [node_map(sm.node).index, class_map(sm.class).index];
end
cutoffs = [];
if isfield(frm, 'cutoffs')
    cutoffs = double(cell2mat_or_vec(frm.cutoffs));
end
map = containers.Map('KeyType', 'char', 'ValueType', 'any');
fn = fieldnames(frm.scaling);
for i = 1:numel(fn)
    nm = fn{i};
    parts = strsplit(nm(2:end), '_');   % drop the 'x' prefix jsondecode prepends
    c = cellfun(@str2double, parts);
    map(cd_state_key(c)) = double(frm.scaling.(nm));
end
g = @(M) firingdep_table_eval(M, slotIdx, map, cutoffs);
end

function v = firingdep_table_eval(M, slotIdx, map, cutoffs)
P = size(slotIdx, 1);
c = zeros(1, P);
for s = 1:P
    c(s) = round(M(slotIdx(s,1), slotIdx(s,2)));
end
c(c < 0) = 0;
if ~isempty(cutoffs)
    m = min(numel(c), numel(cutoffs));
    c(1:m) = min(c(1:m), cutoffs(1:m));
end
k = cd_state_key(c);
if isKey(map, k)
    v = map(k);
else
    v = 1;   % a marking absent from the table is neutral (no dependence)
end
end

function v = cell2mat_or_vec(x)
% jsondecode yields a numeric column for a homogeneous array but a cell when the
% writer emitted num2cell; normalize both to a row vector.
if iscell(x)
    v = cell2mat(x(:)');
else
    v = double(x(:)');
end
end

function v = cd_table_eval(ni, map, cutoffs)
n = round(double(ni(:)'));
n(n < 0) = 0;
if ~isempty(cutoffs)
    m = min(numel(n), numel(cutoffs));
    n(1:m) = min(n(1:m), cutoffs(1:m));
end
k = cd_state_key(n);
if isKey(map, k)
    v = map(k);
else
    v = 1;   % a state absent from the table is neutral (no scaling)
end
end
