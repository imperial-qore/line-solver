function ok = buildLayersPH(self, mode, flat)
% OK = BUILDLAYERSPH(SELF, MODE, FLAT) Build the ensemble of a PH encoding
%
% MODE 'probe' only asks whether the model CAN be served: it runs the feature
% gate and composes the per-entry workflows, which is where a precedence graph
% with no series-parallel reduction is found, and answers true or false without
% touching the ensemble. The alias method 'srvn' uses it to choose between
% 'srvn.ph' and 'srvn.cs'. MODE 'build' (the default) builds the layers, reusing
% the workflows a preceding probe already composed.
%
% FLAT false (the default) is method 'srvn.ph': one layer per served element,
% each a two-station cycle Delay('Clients') + Queue(server), with one closed
% class per caller task. FLAT true is method 'flat.ph': ONE layer holding a
% station for every processor and every called task, with the same one closed
% class per caller task, which now visits each of the servers it uses once per
% invocation instead of meeting them through surrogate delays.
%
% Either way the sequencing that the routing encoding writes as one class per
% entry, per activity and per call, plus Fork, Join, Router and ClassSwitch
% nodes, is composed instead into a single phase-type service law per
% (server, caller), by the exact series-parallel reduction of Workflow.toPH. The
% activity graph survives as a distribution rather than as a topology.
%
% see _kb/06-solver-catalog.md (LN section) for the layering taxonomy
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2 || isempty(mode)
    mode = 'build';
end
if nargin < 3 || isempty(flat)
    flat = false;
end
ok = true;
lqn = self.lqn;
nelem = lqn.nhosts + lqn.ntasks;
mname = 'srvn.ph';
if flat
    mname = 'flat.ph';
end

if strcmp(mode, 'probe')
    % Answer the feasibility question without disturbing the solver: both the
    % feature gate and the series-parallel reduction can refuse, and the second
    % only finds out by composing the workflows, which is work the build then
    % reuses.
    try
        assertSrvnPHSupported(self, lqn, flat);
        self.ph = phInitLaws(self, lqn);
    catch ME
        self.ph = [];
        line_debug('LN: method=%s cannot serve this model (%s)', mname, ME.message);
        ok = false;
    end
    return
end

if isempty(self.ph)
    assertSrvnPHSupported(self, lqn, flat);
end

% The interlock correction rewrites the populations of the call classes, which
% this method does not create: its callers reach the server in one class each
if self.options.config.interlocking
    line_debug('LN: method=%s builds no call classes, interlocking correction disabled', mname);
    self.options.config.interlocking = false;
end

self.ensemble = cell(nelem,1);
self.servt_classes_updmap = cell(nelem,1);
self.call_classes_updmap = cell(nelem,1);
self.arvproc_classes_updmap = cell(nelem,1);
self.thinkt_classes_updmap = cell(nelem,1);
self.actthinkt_classes_updmap = cell(nelem,1);
self.route_prob_updmap = cell(nelem,1);
self.singleReplicaTasks = [];
for i = 1:nelem
    self.servt_classes_updmap{i} = zeros(0,4);
    self.call_classes_updmap{i} = zeros(0,4);
    self.arvproc_classes_updmap{i} = zeros(0,4);
    self.thinkt_classes_updmap{i} = zeros(0,4);
    self.actthinkt_classes_updmap{i} = zeros(0,4);
    self.route_prob_updmap{i} = zeros(0,7);
end

%% Per-entry workflows, and the processor-demand law, which is iteration-invariant.
%% A preceding probe has already composed them; they do not depend on the iterate.
if isempty(self.ph)
    self.ph = phInitLaws(self, lqn);
end

%% Seed the fixed point with the static demands, then compose the entry laws
self.residt = zeros(lqn.nidx,1);
self.servt = zeros(lqn.nidx,1);
self.callservt = zeros(lqn.ncalls,1);
self.callresidt = zeros(lqn.ncalls,1);
for aidx = (lqn.ashift+1):(lqn.ashift+lqn.nacts)
    self.residt(aidx) = lqn.hostdem_mean(aidx);
end
for cidx = 1:lqn.ncalls
    if lqn.calltype(cidx) == CallType.SYNC || lqn.calltype(cidx) == CallType.ASYNC
        eidx = lqn.callpair(cidx,2);
        self.callservt(cidx) = lqn.callproc_mean(cidx) * self.ph.hostmean(eidx);
        self.callresidt(cidx) = self.callservt(cidx);
    end
end
self.phComposeEntryLaws();

if flat
    %% ONE subnetwork holding every processor and every called task
    servers = buildPHFlatLayer(self, lqn);
else
    %% One subnetwork per processor
    for hidx = 1:lqn.nhosts
        if self.ignore(hidx)
            self.ensemble{hidx} = [];
            continue
        end
        callers = hostLayerCallers(self, lqn, hidx);
        if isempty(callers)
            self.ensemble{hidx} = [];
            continue
        end
        buildPHLayer(self, lqn, hidx, callers, true);
    end

    %% One subnetwork per called task
    for t = 1:lqn.ntasks
        tidx = lqn.tshift + t;
        if self.ignore(tidx) || lqn.isref(tidx)
            self.ensemble{tidx} = [];
            continue
        end
        callers = taskLayerCallers(self, lqn, tidx);
        asyncCalls = asyncCallsInto(lqn, tidx);
        if isempty(callers) && isempty(asyncCalls)
            self.ensemble{tidx} = [];
            continue
        end
        buildPHLayer(self, lqn, tidx, callers, false);
    end
end

self.thinkt_classes_updmap = cell2mat(self.thinkt_classes_updmap);
self.actthinkt_classes_updmap = cell2mat(self.actthinkt_classes_updmap);
self.call_classes_updmap = cell2mat(self.call_classes_updmap);
self.servt_classes_updmap = cell2mat(self.servt_classes_updmap);
self.arvproc_classes_updmap = cell2mat(self.arvproc_classes_updmap);
self.route_prob_updmap = cell2mat(self.route_prob_updmap);

if flat
    % every server resolves to the single flat layer, which is at once the host
    % layer and the task layer -- the same convention buildLayers uses
    self.ensemble = self.ensemble(1);
    self.idxhash = NaN(nelem,1);
    self.idxhash(servers) = 1;
    self.hostLayerIndices = 1;
    self.taskLayerIndices = 1;
else
    emptymodels = cellfun(@isempty, self.ensemble);
    self.ensemble(emptymodels) = [];
    self.idxhash = (1:length(emptymodels))' - cumsum(emptymodels);
    self.idxhash(emptymodels) = NaN;

    self.hostLayerIndices = [];
    self.taskLayerIndices = [];
    for hidx = 1:lqn.nhosts
        if ~isnan(self.idxhash(hidx))
            self.hostLayerIndices(end+1) = self.idxhash(hidx);
        end
    end
    for t = 1:lqn.ntasks
        tidx = lqn.tshift + t;
        if ~isnan(self.idxhash(tidx))
            self.taskLayerIndices(end+1) = self.idxhash(tidx);
        end
    end
end

self.layerHasRegion = false(length(self.ensemble),1);
self.layerChains = cell(length(self.ensemble),1);

% install the initial laws, so that iteration 1 sees the seeded demands rather
% than the placeholders the stations were created with
self.thinkt = zeros(lqn.nidx,1);
self.tput = zeros(lqn.nidx,1);
self.updateLayersPH(0);

self.model.ensemble = self.ensemble;
end

% ------------------------------------------------------------------------
function buildPHLayer(self, lqn, idx, callers, ishost)
% Build the two-station layer of server element IDX

model = Network(lqn.hashnames{idx});
model.setChecks(false);
model.attribute = struct('hosts',[],'tasks',[],'entries',[],'activities',[],'calls',[],'serverIdx',0);
model.attribute.serverIdxOf = NaN(lqn.nidx,1);
model.attribute.clientIdx = 1;
model.attribute.serverIdx = 2;
model.attribute.sourceIdx = NaN;

clientDelay = Delay(model, 'Clients');

nreplicas = phReplicaCount(self, lqn, idx, callers, ishost);
mult = lqn.maxmult;
srvStation = cell(1, nreplicas);
for m = 1:nreplicas
    if m == 1
        srvStation{m} = Queue(model, lqn.hashnames{idx}, lqn.sched(idx));
    else
        srvStation{m} = Queue(model, [lqn.hashnames{idx},'.',num2str(m)], lqn.sched(idx));
    end
    srvStation{m}.setNumberOfServers(mult(idx));
    srvStation{m}.attribute.ishost = ishost;
    srvStation{m}.attribute.idx = idx;
end
model.attribute.serverIdxOf(idx) = model.getNodeIndex(srvStation{1});
if ishost
    model.attribute.hostStations = model.attribute.serverIdxOf(idx);
    model.attribute.taskStations = [];
    model.attribute.hosts(end+1,:) = [NaN, model.attribute.serverIdxOf(idx)];
else
    model.attribute.hostStations = [];
    model.attribute.taskStations = model.attribute.serverIdxOf(idx);
    model.attribute.tasks(end+1,:) = [NaN, model.attribute.serverIdxOf(idx)];
end

% --- closed class per caller task
classOfCaller = zeros(lqn.nidx,1);
for c = callers(:)'
    njobs = phLayerPopulation(self, lqn, idx, c, nreplicas);
    self.njobs(c, idx) = njobs;
    cls = ClosedClass(model, lqn.hashnames{c}, njobs, clientDelay);
    cls.setReferenceClass(true);
    cls.attribute = [LayeredNetworkElement.TASK, c];
    classOfCaller(c) = cls.index;
    model.attribute.tasks(end+1,:) = [cls.index, c];
    clientDelay.setService(cls, Exp.fitMean(max(GlobalConstants.FineTol, lqn_ref_thinktime(lqn, c))));
    for m = 1:nreplicas
        srvStation{m}.setService(cls, Exp.fitMean(GlobalConstants.FineTol));
    end
    % every layer must be refreshed after a law change: post() resets the
    % layers named by the think-time map
    self.thinkt_classes_updmap{idx}(end+1,:) = [idx, c, 1, cls.index];
    self.servt_classes_updmap{idx}(end+1,:) = [idx, c, model.getNodeIndex(srvStation{1}), cls.index];
end

% --- open classes: entry arrivals on a host layer, async calls on a task layer
openArrivals = zeros(0,2);  % [class index, entry index] or [class index, -call index]
hasSource = false;
sourceStation = [];
sinkStation = [];
if ishost
    for c = callers(:)'
        % A task no other task calls has no task layer, so updateThinkTimesPH
        % never gives its caller class a surrogate delay: the class cycles against
        % an Immediate one and an open stream on top of it doubles the load. The
        % chain is the representation that honours the thread pool, so it is kept
        % and closed on the arrival rate instead -- see updateThinkTimesPH.
        if phOpenArrivalOnly(lqn, c)
            continue
        end
        for eidx = lqn.entriesof{c}
            if ~phHasOpenArrival(lqn, eidx)
                continue
            end
            [hasSource, sourceStation, sinkStation, model] = ensureSource(hasSource, sourceStation, sinkStation, model);
            ocls = OpenClass(model, [lqn.hashnames{eidx},'.Open'], 0);
            ocls.attribute = [LayeredNetworkElement.ENTRY, eidx];
            sourceStation.setArrival(ocls, lqn.arrival{eidx});
            clientDelay.setService(ocls, Disabled.getInstance());
            for m = 1:nreplicas
                srvStation{m}.setService(ocls, Exp.fitMean(max(GlobalConstants.FineTol, self.ph.hostmean(eidx))));
            end
            openArrivals(end+1,:) = [ocls.index, eidx]; %#ok<AGROW>
            model.attribute.entries(end+1,:) = [ocls.index, eidx];
            self.arvproc_classes_updmap{idx}(end+1,:) = [idx, -eidx, model.getNodeIndex(sourceStation), ocls.index];
        end
    end
else
    for cidx = asyncCallsInto(lqn, idx)
        [hasSource, sourceStation, sinkStation, model] = ensureSource(hasSource, sourceStation, sinkStation, model);
        ocls = OpenClass(model, lqn.callhashnames{cidx}, 0);
        ocls.attribute = [LayeredNetworkElement.CALL, cidx];
        sourceStation.setArrival(ocls, Immediate.getInstance());
        clientDelay.setService(ocls, Disabled.getInstance());
        eidx = lqn.callpair(cidx,2);
        for m = 1:nreplicas
            srvStation{m}.setService(ocls, Exp.fitMean(max(GlobalConstants.FineTol, self.ph.entrymean(eidx))));
        end
        openArrivals(end+1,:) = [ocls.index, -cidx]; %#ok<AGROW>
        model.attribute.calls(end+1,:) = [ocls.index, cidx, lqn.callpair(cidx,1), lqn.callpair(cidx,2)];
        self.arvproc_classes_updmap{idx}(end+1,:) = [idx, cidx, model.getNodeIndex(sourceStation), ocls.index];
        self.call_classes_updmap{idx}(end+1,:) = [idx, cidx, model.getNodeIndex(srvStation{1}), ocls.index];
    end
end

if hasSource
    nClasses = model.getNumberOfClasses();
    for k = 1:nClasses
        if k > length(sourceStation.input.sourceClasses) || isempty(sourceStation.input.sourceClasses{k})
            sourceStation.input.sourceClasses{k} = {[], ServiceStrategy.LI, Disabled.getInstance()};
        end
        if k > length(sourceStation.arrivalProcess) || isempty(sourceStation.arrivalProcess{k})
            sourceStation.arrivalProcess{k} = Disabled.getInstance();
        end
    end
end

%% Routing: one visit to the server per client cycle. The number of calls is
%% carried by the service law, not by a visit ratio, so no arc ever changes
P = model.initRoutingMatrix;
for c = callers(:)'
    k = classOfCaller(c);
    cls = model.classes{k};
    for m = 1:nreplicas
        P{cls, cls}(clientDelay, srvStation{m}) = 1 / nreplicas;
        P{cls, cls}(srvStation{m}, clientDelay) = 1;
    end
end
for r = 1:size(openArrivals,1)
    cls = model.classes{openArrivals(r,1)};
    for m = 1:nreplicas
        P{cls, cls}(sourceStation, srvStation{m}) = 1 / nreplicas;
        P{cls, cls}(srvStation{m}, sinkStation) = 1;
    end
end
model.link(P);

self.ph.layer{idx} = struct('idx', idx, 'ishost', ishost, ...
    'callers', callers(:)', 'classOfCaller', classOfCaller, ...
    'nreplicas', nreplicas, 'clientstation', 1, 'qstations', 1 + (1:nreplicas), ...
    'svcmeanByClass', zeros(model.getNumberOfClasses(),1), 'openArrivals', openArrivals, ...
    'npop', phModelPopulation(self, idx, callers));

self.ensemble{idx} = model;
end

% ------------------------------------------------------------------------
function n = phModelPopulation(self, idx, callers)
% Closed population of the model a server belongs to, i.e. how many jobs a job
% can queue behind. Under 'srvn.ph' that is the callers of the layer's single
% server; under 'flat.ph' it is every caller of the one model.
n = 0;
for c = callers(:)'
    v = self.njobs(c, idx);
    if isfinite(v) && v > 0
        n = n + v;
    end
end
if n < 1
    n = 1;
end
end

% ------------------------------------------------------------------------
function servers = buildPHFlatLayer(self, lqn)
% Build the ONE layer of method 'flat.ph': a client delay plus a station for
% every processor and every called task.
%
% A caller task is one closed class, and it visits each server it uses ONCE per
% invocation, carrying there the composed law of the demand it places on that
% server -- the same law method 'srvn.ph' installs in the server's own layer.
% What changes is that the servers now contend inside one network instead of
% seeing each other through surrogate delays, so the client delay keeps only the
% think times and whatever of the cycle this model does not hold. That is the
% whole difference between the two encodings of the PH composition, and it is
% why the reconstruction passes are shared verbatim.

servers = phFlatServerSet(self, lqn);

model = Network('FlatPH');
model.setChecks(false);
model.attribute = struct('hosts',[],'tasks',[],'entries',[],'activities',[],'calls',[],'serverIdx',0);
model.attribute.serverIdxOf = NaN(lqn.nidx,1);
model.attribute.clientIdx = 1;
model.attribute.sourceIdx = NaN;
model.attribute.hostStations = [];
model.attribute.taskStations = [];

clientDelay = Delay(model, 'Clients');

nsrv = numel(servers);
mult = lqn.maxmult;
srvStation = cell(1, nsrv);
for s = 1:nsrv
    idx = servers(s);
    ishost = idx <= lqn.nhosts;
    q = Queue(model, lqn.hashnames{idx}, lqn.sched(idx));
    q.setNumberOfServers(mult(idx));
    q.attribute.ishost = ishost;
    q.attribute.idx = idx;
    srvStation{s} = q;
    stn = model.getNodeIndex(q);
    model.attribute.serverIdxOf(idx) = stn;
    if ishost
        model.attribute.hostStations(end+1) = stn;
        model.attribute.hosts(end+1,:) = [NaN, stn];
    else
        model.attribute.taskStations(end+1) = stn;
        model.attribute.tasks(end+1,:) = [NaN, stn];
    end
end
% the scalar fallback of stationIdxOf, which no served element ever reaches here
model.attribute.serverIdx = model.attribute.serverIdxOf(servers(1));

%% Callers of each server, and the union of them, which becomes the class set
callersOf = cell(lqn.nidx,1);
for s = 1:nsrv
    idx = servers(s);
    if idx <= lqn.nhosts
        callersOf{idx} = hostLayerCallers(self, lqn, idx);
    else
        callersOf{idx} = taskLayerCallers(self, lqn, idx);
    end
end
allCallers = [];
for s = 1:nsrv
    allCallers = union(allCallers, callersOf{servers(s)});
end
allCallers = reshape(allCallers, 1, []);

%% One closed class per caller task
classOfCaller = zeros(lqn.nidx,1);
npop = 0;
for c = allCallers
    % phFlatServerSet has refused every replicated element, so the per-replica
    % reduction the srvn builder makes is the identity here
    njobs = phLayerPopulation(self, lqn, servers(1), c, 1);
    cls = ClosedClass(model, lqn.hashnames{c}, njobs, clientDelay);
    cls.setReferenceClass(true);
    cls.attribute = [LayeredNetworkElement.TASK, c];
    classOfCaller(c) = cls.index;
    model.attribute.tasks(end+1,:) = [cls.index, c];
    npop = npop + njobs;
    clientDelay.setService(cls, Exp.fitMean(max(GlobalConstants.FineTol, lqn_ref_thinktime(lqn, c))));
    % A station this caller never reaches must say so with Disabled, NOT with a
    % tiny placeholder law. An FCFS station carries ONE service law across its
    % classes, so a placeholder is not inert there: it is mixed into the
    % multiserver correction and invents waiting where there is none. On the
    % two-task probe with a 5-server callee and 5 callers, which cannot queue at
    % all, the placeholder reported R/S = 1.2513 instead of 1.0000, and the
    % error compounded through the fixed point into 7 to 13 per cent on the
    % reference throughput. Under 'srvn.ph' the question never arises, since
    % every class of a layer visits that layer's single server.
    for s = 1:nsrv
        srvStation{s}.setService(cls, Disabled.getInstance());
    end
    for s = 1:nsrv
        idx = servers(s);
        if ~any(callersOf{idx} == c)
            continue
        end
        srvStation{s}.setService(cls, Exp.fitMean(GlobalConstants.FineTol));
        self.njobs(c, idx) = njobs;
        % every layer must be refreshed after a law change: post() resets the
        % layers named by the think-time map
        self.thinkt_classes_updmap{idx}(end+1,:) = [idx, c, 1, cls.index];
        self.servt_classes_updmap{idx}(end+1,:) = [idx, c, model.getNodeIndex(srvStation{s}), cls.index];
    end
end

%% Open classes: entry arrivals on a processor station, async calls on a task one
openArrivalsOf = cell(lqn.nidx,1);
for s = 1:nsrv
    openArrivalsOf{servers(s)} = zeros(0,2);
end
hasSource = false;
sourceStation = [];
sinkStation = [];
for s = 1:nsrv
    hidx = servers(s);
    if hidx > lqn.nhosts
        continue
    end
    for c = callersOf{hidx}
        % A task no other task calls has no task station, so updateThinkTimesPH
        % never gives its caller class a surrogate delay: the class cycles
        % against an Immediate one and an open stream on top of it doubles the
        % load. The chain is the representation that honours the thread pool, so
        % it is kept and closed on the arrival rate instead.
        if phOpenArrivalOnly(lqn, c)
            continue
        end
        for eidx = lqn.entriesof{c}
            if ~phHasOpenArrival(lqn, eidx)
                continue
            end
            [hasSource, sourceStation, sinkStation, model] = ensureSource(hasSource, sourceStation, sinkStation, model);
            ocls = OpenClass(model, [lqn.hashnames{eidx},'.Open'], 0);
            ocls.attribute = [LayeredNetworkElement.ENTRY, eidx];
            sourceStation.setArrival(ocls, lqn.arrival{eidx});
            clientDelay.setService(ocls, Disabled.getInstance());
            % Disabled, not a placeholder, at every station this stream misses
            for s2 = 1:nsrv
                srvStation{s2}.setService(ocls, Disabled.getInstance());
            end
            srvStation{s}.setService(ocls, Exp.fitMean(max(GlobalConstants.FineTol, self.ph.hostmean(eidx))));
            openArrivalsOf{hidx}(end+1,:) = [ocls.index, eidx];
            model.attribute.entries(end+1,:) = [ocls.index, eidx];
            self.arvproc_classes_updmap{hidx}(end+1,:) = [hidx, -eidx, model.getNodeIndex(sourceStation), ocls.index];
        end
    end
end
for s = 1:nsrv
    tidx = servers(s);
    if tidx <= lqn.nhosts
        continue
    end
    for cidx = asyncCallsInto(lqn, tidx)
        [hasSource, sourceStation, sinkStation, model] = ensureSource(hasSource, sourceStation, sinkStation, model);
        ocls = OpenClass(model, lqn.callhashnames{cidx}, 0);
        ocls.attribute = [LayeredNetworkElement.CALL, cidx];
        sourceStation.setArrival(ocls, Immediate.getInstance());
        clientDelay.setService(ocls, Disabled.getInstance());
        % Disabled, not a placeholder, at every station this stream misses
        for s2 = 1:nsrv
            srvStation{s2}.setService(ocls, Disabled.getInstance());
        end
        eidx = lqn.callpair(cidx,2);
        srvStation{s}.setService(ocls, Exp.fitMean(max(GlobalConstants.FineTol, self.ph.entrymean(eidx))));
        openArrivalsOf{tidx}(end+1,:) = [ocls.index, -cidx];
        model.attribute.calls(end+1,:) = [ocls.index, cidx, lqn.callpair(cidx,1), lqn.callpair(cidx,2)];
        self.arvproc_classes_updmap{tidx}(end+1,:) = [tidx, cidx, model.getNodeIndex(sourceStation), ocls.index];
        self.call_classes_updmap{tidx}(end+1,:) = [tidx, cidx, model.getNodeIndex(srvStation{s}), ocls.index];
    end
end

if hasSource
    nClasses = model.getNumberOfClasses();
    for k = 1:nClasses
        if k > length(sourceStation.input.sourceClasses) || isempty(sourceStation.input.sourceClasses{k})
            sourceStation.input.sourceClasses{k} = {[], ServiceStrategy.LI, Disabled.getInstance()};
        end
        if k > length(sourceStation.arrivalProcess) || isempty(sourceStation.arrivalProcess{k})
            sourceStation.arrivalProcess{k} = Disabled.getInstance();
        end
    end
end

%% Routing: one visit per server the caller uses, in server order. The number of
%% calls is carried by the service law, not by a visit ratio, so no arc ever moves
P = model.initRoutingMatrix;
for c = allCallers
    cls = model.classes{classOfCaller(c)};
    prev = clientDelay;
    visited = false;
    for s = 1:nsrv
        if ~any(callersOf{servers(s)} == c)
            continue
        end
        P{cls, cls}(prev, srvStation{s}) = 1;
        prev = srvStation{s};
        visited = true;
    end
    if visited
        P{cls, cls}(prev, clientDelay) = 1;
    end
end
for s = 1:nsrv
    oa = openArrivalsOf{servers(s)};
    for r = 1:size(oa,1)
        cls = model.classes{oa(r,1)};
        P{cls, cls}(sourceStation, srvStation{s}) = 1;
        P{cls, cls}(srvStation{s}, sinkStation) = 1;
    end
end
model.link(P);

nclasses = model.getNumberOfClasses();
for s = 1:nsrv
    idx = servers(s);
    self.ph.layer{idx} = struct('idx', idx, 'ishost', idx <= lqn.nhosts, ...
        'callers', reshape(callersOf{idx},1,[]), 'classOfCaller', classOfCaller, ...
        'nreplicas', 1, 'clientstation', 1, ...
        'qstations', model.getNodeIndex(srvStation{s}), ...
        'svcmeanByClass', zeros(nclasses,1), 'openArrivals', openArrivalsOf{idx}, ...
        'npop', max(1, npop));
end

self.ensemble{1} = model;
end

% ------------------------------------------------------------------------
function servers = phFlatServerSet(self, lqn)
% Processors and called tasks that become stations of the flat layer.
%
% The set is the elements the srvn builder would have given a layer of their
% own, so 'flat.ph' and 'srvn.ph' place the SAME stations and differ only in how
% many networks hold them. The squashing refusals (a replica, a cache task, a
% setup task: per-layer state that one submodel cannot hold) are stated in
% ln_method_refusal, which assertSrvnPHSupported asks before this runs.
servers = [];
for hidx = 1:lqn.nhosts
    if self.ignore(hidx) || isempty(lqn.tasksof{hidx})
        continue
    end
    if isempty(hostLayerCallers(self, lqn, hidx))
        continue
    end
    servers(end+1) = hidx; %#ok<AGROW>
end
for t = 1:lqn.ntasks
    tidx = lqn.tshift + t;
    if self.ignore(tidx) || lqn.isref(tidx)
        continue
    end
    if isempty(taskLayerCallers(self, lqn, tidx)) && isempty(asyncCallsInto(lqn, tidx))
        continue
    end
    servers(end+1) = tidx; %#ok<AGROW>
end
if isempty(servers)
    line_error(mfilename, 'method=''flat.ph'' found no server: the model has no processor with tasks.');
end
end

% ------------------------------------------------------------------------
function [hasSource, sourceStation, sinkStation, model] = ensureSource(hasSource, sourceStation, sinkStation, model)
if ~hasSource
    hasSource = true;
    model.attribute.sourceIdx = length(model.nodes)+1;
    sourceStation = Source(model,'Source');
    sinkStation = Sink(model,'Sink');
end
end

% ------------------------------------------------------------------------
function ph = phInitLaws(self, lqn)
% Build the per-entry workflows and the iteration-invariant processor law
ph = struct();
ph.wf = cell(lqn.nidx,1);
ph.wfhost = cell(lqn.nidx,1);
ph.execs = cell(lqn.nidx,1);
ph.callexecs = cell(lqn.nidx,1);
ph.hostalpha = cell(lqn.nidx,1);
ph.hostT = cell(lqn.nidx,1);
ph.hostmean = zeros(lqn.nidx,1);
ph.entryalpha = cell(lqn.nidx,1);
ph.entryT = cell(lqn.nidx,1);
ph.entrymean = zeros(lqn.nidx,1);
ph.entryscv = ones(lqn.nidx,1);
ph.share = zeros(lqn.nidx,1);
ph.xdemand = zeros(lqn.nidx,1); % requests the callers ask of each task
ph.ncalls = zeros(lqn.nidx, lqn.nidx); % [caller task, called entry] calls per invocation
ph.layer = cell(lqn.nhosts+lqn.ntasks,1);
ph.callresp = zeros(max(lqn.ncalls,1),1);
ph.callscv = ones(max(lqn.ncalls,1),1);

for e = 1:lqn.nentries
    eidx = lqn.eshift + e;
    tidx = lqn.parent(eidx);
    if self.ignore(tidx)
        continue
    end
    [ph.wf{eidx}, ~, ~, ph.execs{eidx}, ph.callexecs{eidx}] = ...
        lqn_entry_workflow(self.model, lqn, eidx, true);
    ph.wfhost{eidx} = lqn_entry_workflow(self.model, lqn, eidx, false);
    % the processor sees the WORK of concurrent branches, not their elapsed
    % time, so the host law serialises an AND fork -- see lqn_ph_serial_law
    [ph.hostalpha{eidx}, ph.hostT{eidx}] = lqn_ph_serial_law(ph.wfhost{eidx});
    ph.hostmean(eidx) = lqn_ph_moments(ph.hostalpha{eidx}, ph.hostT{eidx});
end

% until the first iteration reports throughputs, a task splits its requests
% evenly over its entries
for t = 1:lqn.ntasks
    tidx = lqn.tshift + t;
    entries = lqn.entriesof{tidx};
    if isempty(entries)
        continue
    end
    for eidx = entries
        ph.share(eidx) = 1 / numel(entries);
    end
end
end

% ------------------------------------------------------------------------
function callers = hostLayerCallers(self, lqn, hidx)
% Tasks that run on processor HIDX and reach it with requests
callers = [];
for tidx = lqn.tasksof{hidx}
    if self.ignore(tidx)
        continue
    end
    if lqn.isref(tidx)
        callers(end+1) = tidx; %#ok<AGROW>
        continue
    end
    served = false;
    for eidx = lqn.entriesof{tidx}
        if any(full(lqn.issynccaller(:, eidx))) || any(full(lqn.isasynccaller(:, eidx))) ...
                || phHasOpenArrival(lqn, eidx)
            served = true;
            break
        end
    end
    if served
        callers(end+1) = tidx; %#ok<AGROW>
    end
end
end

% ------------------------------------------------------------------------
function callers = taskLayerCallers(self, lqn, tidx)
% Tasks issuing a synchronous call to an entry of TIDX
callers = [];
for c = lqn.tshift + (1:lqn.ntasks)
    if c == tidx || self.ignore(c)
        continue
    end
    if any(full(lqn.issynccaller(c, lqn.entriesof{tidx})))
        callers(end+1) = c; %#ok<AGROW>
    end
end
end

% ------------------------------------------------------------------------
function cidxs = asyncCallsInto(lqn, tidx)
% Asynchronous calls whose target entry belongs to TIDX
cidxs = [];
if lqn.ncalls == 0
    return
end
targets = lqn.entriesof{tidx};
for cidx = 1:lqn.ncalls
    if lqn.calltype(cidx) == CallType.ASYNC && any(targets == lqn.callpair(cidx,2))
        cidxs(end+1) = cidx; %#ok<AGROW>
    end
end
end

% ------------------------------------------------------------------------
function tf = phOpenArrivalOnly(lqn, tidx)
% True when an entry arrival is the ONLY way requests reach task TIDX. 'srvn.ph'
% refuses forwarding calls outright, so sync/async callers are the whole test.
tf = false;
if lqn.isref(tidx)
    return
end
for eidx = lqn.entriesof{tidx}
    if any(full(lqn.issynccaller(:, eidx))) || any(full(lqn.isasynccaller(:, eidx)))
        return
    end
end
for eidx = lqn.entriesof{tidx}
    if phHasOpenArrival(lqn, eidx)
        tf = true;
        return
    end
end
end

% ------------------------------------------------------------------------
function tf = phHasOpenArrival(lqn, eidx)
tf = isfield(lqn,'arrival') && ~isempty(lqn.arrival) && iscell(lqn.arrival) && ...
    eidx <= numel(lqn.arrival) && ~isempty(lqn.arrival{eidx});
end

% ------------------------------------------------------------------------
function nreplicas = phReplicaCount(self, lqn, idx, callers, ishost)
% Replicas of the server station, with the same fan-out reduction as the
% default builder: a caller that reaches every replica sees one representative
rawReplicas = lqn.repl(idx);
nreplicas = rawReplicas;
if rawReplicas <= 1 || isempty(callers)
    nreplicas = max(1, rawReplicas);
    return
end
reduceFanout = false;
if ~ishost && isfield(lqn,'fanout') && ~isempty(lqn.fanout)
    reduceFanout = true;
    for c = callers(:)'
        if lqn.fanout(c, idx) < rawReplicas
            reduceFanout = false;
            break
        end
    end
elseif ishost
    reduceFanout = true;
    for c = callers(:)'
        if lqn.repl(c) ~= rawReplicas
            reduceFanout = false;
            break
        end
    end
end
if reduceFanout
    nreplicas = 1;
    if ~ishost
        self.singleReplicaTasks(end+1) = idx;
    end
end
end

% ------------------------------------------------------------------------
function njobs = phLayerPopulation(self, lqn, idx, c, nreplicas)
% Threads of caller C present in the layer of IDX
mult = lqn.maxmult;
callerIsSingleReplica = nreplicas == 1 && lqn.repl(idx) > 1;
callerIsSingleReplica = callerIsSingleReplica || any(self.singleReplicaTasks == c);
if callerIsSingleReplica
    njobs = mult(c);
else
    njobs = mult(c) * lqn.repl(c);
end
if isinf(njobs)
    callers_of_c = find(lqn.taskgraph(:,c));
    njobs = sum(mult(callers_of_c)); %#ok<FNDSB>
    if isinf(njobs)
        njobs = min(sum(mult(isfinite(mult)) .* lqn.repl(isfinite(mult))), 1000);
    end
end
end

% ------------------------------------------------------------------------
function assertSrvnPHSupported(self, lqn, flat) %#ok<INUSL>
% Features the composed law cannot represent are refused by name rather than
% silently degraded -- see _kb/06-solver-catalog.md (LN section). The rules
% (a second phase, a forwarding call, a cache task, a setup on an INF task, a
% routed call group, an admission constraint, a queue-dependent rate, and under
% 'flat.ph' the squashing refusals) live in LN_METHOD_REFUSAL, the predicate
% SolverLN.supportsModelMethod asks, so the gate and this run speak one sentence.
if nargin < 3 || isempty(flat)
    flat = false;
end
mname = 'srvn.ph';
if flat
    mname = 'flat.ph';
end
[ok, reason] = ln_method_refusal(lqn, mname);
if ~ok
    line_error(mfilename, reason);
end
end
