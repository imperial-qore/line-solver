function lsn = getStruct(self, regenerate)
% LSN = GETSTRUCT(SELF, regenerate)
%
%
% Copyright 2012-2026, Imperial College London
if nargin<2
    regenerate = false;
end
if ~isempty(self.lsn) && ~regenerate
    lsn = self.lsn;
    return
end
lsn = LayeredNetworkStruct();
lsn.nidx = 0;  % total number of hosts, tasks, entries, and activities, except the reference tasks
lsn.hshift = 0;
lsn.nhosts = length(self.hosts);
lsn.ntasks = length(self.tasks);
lsn.nentries = length(self.entries);
lsn.nacts = length(self.activities);
lsn.tshift = lsn.nhosts;
lsn.eshift = lsn.nhosts + lsn.ntasks;
lsn.ashift = lsn.nhosts + lsn.ntasks + lsn.nentries;
lsn.cshift = lsn.nhosts + lsn.ntasks + lsn.nentries + lsn.nacts;

%% analyze static properties
lsn.nidx = lsn.nhosts + lsn.ntasks + lsn.nentries + lsn.nacts;
idx = 1;
lsn.tasksof = cell(lsn.nhosts,1);
lsn.entriesof = cell(lsn.nhosts+lsn.ntasks,1);
lsn.actsof = cell(lsn.nhosts+lsn.ntasks+lsn.nentries,1);
lsn.callsof = cell(lsn.nacts,1);
lsn.hostdem = {};
lsn.hostdem_type = zeros(lsn.nidx, 1);
lsn.hostdem_params = cell(lsn.nidx, 1);
lsn.hostdem_mean = nan(lsn.nidx, 1);
lsn.hostdem_scv = nan(lsn.nidx, 1);
lsn.hostdem_proc = cell(lsn.nidx, 1);

lsn.actthink = {};
lsn.actthink_type = zeros(lsn.nidx, 1);
lsn.actthink_params = cell(lsn.nidx, 1);
lsn.actthink_mean = nan(lsn.nidx, 1);
lsn.actthink_scv = nan(lsn.nidx, 1);
lsn.actthink_proc = cell(lsn.nidx, 1);

lsn.think = {};
lsn.think_type = zeros(lsn.nidx, 1);
lsn.think_params = cell(lsn.nidx, 1);
lsn.think_mean = nan(lsn.nidx, 1);
lsn.think_scv = nan(lsn.nidx, 1);
lsn.think_proc = cell(lsn.nidx, 1);
lsn.sched = [];
lsn.names = {};
lsn.hashnames = {};
%lsn.shortnames = {};
lsn.mult = zeros(lsn.nhosts+lsn.ntasks,1);
lsn.maxmult = zeros(lsn.nhosts+lsn.ntasks,1);
lsn.repl = zeros(lsn.nhosts+lsn.ntasks,1);
lsn.type = zeros(lsn.nidx,1);
lsn.graph = zeros(lsn.nidx,lsn.nidx);
loop_back_edges = false(lsn.nidx,lsn.nidx);
%lsn.replies = [];
lsn.replygraph = false(lsn.nacts,lsn.nentries);
lsn.actphase = ones(lsn.nacts,1);  % Phase for each activity (default=1)

lsn.nitems = zeros(lsn.nhosts+lsn.ntasks+lsn.nentries,1);
lsn.itemcap  = {};
lsn.itemproc = {};
lsn.itemproc_type = zeros(lsn.nidx, 1);
lsn.itemproc_params = cell(lsn.nidx, 1);
lsn.itemproc_mean = nan(lsn.nidx, 1);
lsn.itemproc_scv = nan(lsn.nidx, 1);
lsn.itemproc_proc = cell(lsn.nidx, 1);

lsn.iscache = zeros(lsn.nhosts+lsn.ntasks,1);
lsn.setuptime = cell(lsn.nhosts+lsn.ntasks,1);
lsn.setuptime_type = zeros(lsn.nidx, 1);
lsn.setuptime_params = cell(lsn.nidx, 1);
lsn.setuptime_mean = nan(lsn.nidx, 1);
lsn.setuptime_scv = nan(lsn.nidx, 1);
lsn.setuptime_proc = cell(lsn.nidx, 1);

lsn.delayofftime = {};
lsn.delayofftime_type = zeros(lsn.nidx, 1);
lsn.delayofftime_params = cell(lsn.nidx, 1);
lsn.delayofftime_mean = nan(lsn.nidx, 1);
lsn.delayofftime_scv = nan(lsn.nidx, 1);
lsn.delayofftime_proc = cell(lsn.nidx, 1);

lsn.isfunction = zeros(lsn.nhosts+lsn.ntasks,1);

% Open arrival distributions
lsn.arrival = {};
lsn.arrival_type = zeros(lsn.nidx, 1);
lsn.arrival_params = cell(lsn.nidx, 1);
lsn.arrival_mean = nan(lsn.nidx, 1);
lsn.arrival_scv = nan(lsn.nidx, 1);
lsn.arrival_proc = cell(lsn.nidx, 1);

lsn.parent = [];
for p=1:lsn.nhosts  % for every processor, scheduling, multiplicity, replication, names, type
    lsn.sched(idx,1) = SchedStrategy.fromText(self.hosts{p}.scheduling);
    lsn.mult(idx,1) = self.hosts{p}.multiplicity;
    lsn.repl(idx,1) = self.hosts{p}.replication;
    lsn.names{idx,1} = self.hosts{p}.name;
    lsn.hashnames{idx,1} = ['P:',lsn.names{idx,1}];
    %lsn.shortnames{idx,1} = ['P',num2str(p)];
    lsn.type(idx,1) = LayeredNetworkElement.HOST; % processor
    idx = idx + 1;
end

for t=1:lsn.ntasks
    lsn.sched(idx,1) = SchedStrategy.fromText(self.tasks{t}.scheduling);
    % hostdem for tasks is Immediate
    immDist = Immediate.getInstance();
    lsn.hostdem{idx,1} = immDist;
    [lsn.hostdem_type(idx), lsn.hostdem_params{idx}, lsn.hostdem_mean(idx), lsn.hostdem_scv(idx), lsn.hostdem_proc{idx}] = extractDistParams(immDist);
    % think time
    thinkDist = self.tasks{t}.thinkTime;
    lsn.think{idx,1} = thinkDist;
    [lsn.think_type(idx), lsn.think_params{idx}, lsn.think_mean(idx), lsn.think_scv(idx), lsn.think_proc{idx}] = extractDistParams(thinkDist);
    lsn.mult(idx,1) = self.tasks{t}.multiplicity;
    lsn.repl(idx,1) = self.tasks{t}.replication;
    lsn.names{idx,1} = self.tasks{t}.name;
    switch lsn.sched(idx,1)
        case SchedStrategy.REF
            lsn.hashnames{idx,1} = ['R:',lsn.names{idx,1}];
            %lsn.shortnames{idx,1} = ['R',num2str(idx-tshift)];
        otherwise
            lsn.hashnames{idx,1} = ['T:',lsn.names{idx,1}];
            %lsn.shortnames{idx,1} = ['T',num2str(idx-tshift)];
    end
    switch class(self.tasks{t})
        case 'CacheTask'
            lsn.nitems(idx,1) = self.tasks{t}.items;
            lsn.itemcap{idx,1} = self.tasks{t}.itemLevelCap;
            lsn.replacestrat(idx,1) = self.tasks{t}.replacestrategy;
            lsn.hashnames{idx,1} = ['C:',lsn.names{idx,1}];
        case 'FunctionTask'
            setupDist = self.tasks{t}.SetupTime;
            lsn.setuptime{idx,1} = setupDist;
            [lsn.setuptime_type(idx), lsn.setuptime_params{idx}, lsn.setuptime_mean(idx), lsn.setuptime_scv(idx), lsn.setuptime_proc{idx}] = extractDistParams(setupDist);
            delayoffDist = self.tasks{t}.DelayOffTime;
            lsn.delayofftime{idx,1} = delayoffDist;
            [lsn.delayofftime_type(idx), lsn.delayofftime_params{idx}, lsn.delayofftime_mean(idx), lsn.delayofftime_scv(idx), lsn.delayofftime_proc{idx}] = extractDistParams(delayoffDist);
            lsn.hashnames{idx,1} = ['F:',lsn.names{idx,1}];
            %lsn.shortnames{idx,1} = ['C',num2str(idx-tshift)];
    end
    pidx = find(cellfun(@(x) strcmp(x.name, self.tasks{t}.parent.name), self.hosts));
    lsn.parent(idx) = pidx;
    lsn.graph(idx, pidx) = 1;
    lsn.type(idx) = LayeredNetworkElement.TASK; % task
    idx = idx + 1;
end

for p=1:lsn.nhosts  % for every processor
    pidx = p;
    lsn.tasksof{pidx} = find(lsn.parent == pidx);
end

for e=1:lsn.nentries
    lsn.names{idx,1} = self.entries{e}.name;

    % Extract open arrival distribution if present
    if ~isempty(self.entries{e}.arrival) && isa(self.entries{e}.arrival, 'Distribution')
        arrDist = self.entries{e}.arrival;
        lsn.arrival{idx,1} = arrDist;
        [lsn.arrival_type(idx), lsn.arrival_params{idx}, lsn.arrival_mean(idx), lsn.arrival_scv(idx), lsn.arrival_proc{idx}] = extractDistParams(arrDist);
    end

    switch class(self.entries{e})
        case 'Entry'
            lsn.hashnames{idx,1} = ['E:',lsn.names{idx,1}];
            for a=1:length(self.entries{e}.replyActivity)
                ractname = self.entries{e}.replyActivity{a};
                ractidx = find(cellfun(@(x) strcmp(x.name, ractname), self.activities));
                lsn.replygraph(ractidx,e)=true;
            end
            %lsn.shortnames{idx,1} = ['E',num2str(idx-eshift)];
        case 'ItemEntry'
            lsn.hashnames{idx,1} = ['I:',lsn.names{idx,1}];
            %lsn.shortnames{idx,1} = ['I',num2str(idx-eshift)];
            lsn.nitems(idx,1) = self.entries{e}.cardinality;
            itemDist = self.entries{e}.popularity;
            lsn.itemproc{idx,1} = itemDist;
            [lsn.itemproc_type(idx), lsn.itemproc_params{idx}, lsn.itemproc_mean(idx), lsn.itemproc_scv(idx), lsn.itemproc_proc{idx}] = extractDistParams(itemDist);
    end
    % hostdem for entries is Immediate
    immDist = Immediate.getInstance();
    lsn.hostdem{idx,1} = immDist;
    [lsn.hostdem_type(idx), lsn.hostdem_params{idx}, lsn.hostdem_mean(idx), lsn.hostdem_scv(idx), lsn.hostdem_proc{idx}] = extractDistParams(immDist);
    tidx = lsn.nhosts + find(cellfun(@(x) strcmp(x.name, self.entries{e}.parent.name), self.tasks));
    lsn.parent(idx) = tidx;
    lsn.graph(tidx,idx) = 1;
    lsn.entriesof{tidx}(end+1) = idx;
    lsn.type(idx) = LayeredNetworkElement.ENTRY; % entries
    idx = idx + 1;
end

for a=1:lsn.nacts
    lsn.names{idx,1} = self.activities{a}.name;
    lsn.hashnames{idx,1} = ['A:',lsn.names{idx,1}];
    %lsn.shortnames{idx,1} = ['A',num2str(idx - lsn.ashift)];
    % hostDemand for activities
    hostdemDist = self.activities{a}.hostDemand;
    lsn.hostdem{idx,1} = hostdemDist;
    [lsn.hostdem_type(idx), lsn.hostdem_params{idx}, lsn.hostdem_mean(idx), lsn.hostdem_scv(idx), lsn.hostdem_proc{idx}] = extractDistParams(hostdemDist);
    % thinkTime for activities
    actthinkDist = self.activities{a}.thinkTime;
    lsn.actthink{idx,1} = actthinkDist;
    [lsn.actthink_type(idx), lsn.actthink_params{idx}, lsn.actthink_mean(idx), lsn.actthink_scv(idx), lsn.actthink_proc{idx}] = extractDistParams(actthinkDist);
    tidx = lsn.nhosts + find(cellfun(@(x) strcmp(x.name, self.activities{a}.parent.name), self.tasks));
    lsn.parent(idx) = tidx;
    lsn.actsof{tidx}(end+1) = idx;
    lsn.type(idx) = LayeredNetworkElement.ACTIVITY; % activities
    lsn.actphase(a) = self.activities{a}.phase;  % Store activity phase
    idx = idx + 1;
end

nidx = idx - 1; % number of indices
lsn.graph(nidx,nidx) = 0;

tasks = self.tasks;
%% now analyze calls
cidx = 0;
lsn.calltype = sparse([],lsn.nidx,1);
lsn.iscaller = sparse(lsn.nidx,lsn.nidx);
lsn.issynccaller = sparse(lsn.nidx,lsn.nidx);
lsn.isasynccaller = sparse(lsn.nidx,lsn.nidx);
lsn.callpair = [];
lsn.callproc = {};
lsn.callproc_type = [];
lsn.callproc_params = {};
lsn.callproc_mean = [];
lsn.callproc_scv = [];
lsn.callproc_proc = {};
lsn.callnames = {};
lsn.callhashnames = {};
%lsn.callshortnames = {};
lsn.taskgraph = sparse(lsn.tshift+lsn.ntasks, lsn.tshift+lsn.ntasks);
lsn.actpretype = sparse(lsn.nidx,1);
lsn.actposttype = sparse(lsn.nidx,1);

% Track boundToEntry mappings to validate uniqueness
boundEntries = {};
boundActivities = {};

for t = 1:lsn.ntasks
    tidx = lsn.tshift+t;
    for a=1:length(self.tasks{t}.activities)
        aidx = findstring(lsn.hashnames, ['A:',tasks{t}.activities(a).name]);
        lsn.callsof{aidx} = [];
        boundToEntry = tasks{t}.activities(a).boundToEntry;
        %for b=1:length(boundToEntry)
        eidx = findstring(lsn.hashnames, ['E:',boundToEntry]);
        if eidx<0
            eidx = findstring(lsn.hashnames, ['I:',boundToEntry]);
        end
        if eidx>0
            lsn.graph(eidx, aidx) = 1;
            
            % Check if this entry is already bound to another activity
            activityName = tasks{t}.activities(a).name;
            existingIdx = find(strcmp(boundEntries, boundToEntry), 1);
            if ~isempty(existingIdx)
                line_error(mfilename, sprintf('Multiple activities (%s, %s) are bound to the same entry: %s', ...
                    boundActivities{existingIdx}, activityName, boundToEntry));
            else
                boundEntries{end+1} = boundToEntry;
                boundActivities{end+1} = activityName;
            end
        end
        %end

        for s=1:length(tasks{t}.activities(a).syncCallDests)
            target_eidx = findstring(lsn.hashnames, ['E:',tasks{t}.activities(a).syncCallDests{s}]);
            if target_eidx < 0
                target_eidx = findstring(lsn.hashnames, ['I:',tasks{t}.activities(a).syncCallDests{s}]);
            end
            target_tidx = lsn.parent(target_eidx);
            cidx = cidx + 1;
            lsn.calltype(cidx,1) = CallType.SYNC;
            lsn.callpair(cidx,1:2) = [aidx,target_eidx];
            if tidx == target_tidx
                line_error(mfilename, 'An entry on a task cannot call another entry on the same task.');
            end
            lsn.callnames{cidx,1} = [lsn.names{aidx},'=>',lsn.names{target_eidx}];
            lsn.callhashnames{cidx,1} = [lsn.hashnames{aidx},'=>',lsn.hashnames{target_eidx}];
            %lsn.callshortnames{cidx,1} = [lsn.shortnames{aidx},'=>',lsn.shortnames{target_eidx}];
            callDist = Geometric(1/tasks{t}.activities(a).syncCallMeans(s)); % synch
            lsn.callproc{cidx,1} = callDist;
            [lsn.callproc_type(cidx), lsn.callproc_params{cidx}, lsn.callproc_mean(cidx), lsn.callproc_scv(cidx), lsn.callproc_proc{cidx}] = extractDistParams(callDist);
            lsn.callsof{aidx}(end+1) = cidx;
            lsn.iscaller(tidx, target_tidx) = true;
            lsn.iscaller(aidx, target_tidx) = true;
            lsn.iscaller(tidx, target_eidx) = true;
            lsn.iscaller(aidx, target_eidx) = true;
            lsn.issynccaller(tidx, target_tidx) = true;
            lsn.issynccaller(aidx, target_tidx) = true;
            lsn.issynccaller(tidx, target_eidx) = true;
            lsn.issynccaller(aidx, target_eidx) = true;
            lsn.taskgraph(tidx, target_tidx) = 1;
            lsn.graph(aidx, target_eidx) = 1;
        end

        for s=1:length(tasks{t}.activities(a).asyncCallDests)
            target_entry_name = tasks{t}.activities(a).asyncCallDests{s};
            target_eidx = findstring(lsn.hashnames,['E:',target_entry_name]);
            if target_eidx < 0
                target_eidx = findstring(lsn.hashnames, ['I:',target_entry_name]);
            end
            % Validate that the target entry exists
            if target_eidx <= 0
                line_error(mfilename, sprintf('Activity "%s" has an async call to non-existent entry "%s".', tasks{t}.activities(a).name, target_entry_name));
            end
            target_tidx = lsn.parent(target_eidx);
            % Check for self-referential async calls (task calling itself)
            if tidx == target_tidx
                line_error(mfilename, sprintf('Activity "%s" in task "%s" has an async call to an entry on the same task. Async self-calls are not supported.', tasks{t}.activities(a).name, tasks{t}.name));
            end
            cidx = cidx + 1;
            lsn.callpair(cidx,1:2) = [aidx,target_eidx];
            lsn.calltype(cidx,1) = CallType.ASYNC; % async
            lsn.callnames{cidx,1} = [lsn.names{aidx},'->',lsn.names{target_eidx}];
            lsn.callhashnames{cidx,1} = [lsn.hashnames{aidx},'->',lsn.hashnames{target_eidx}];
            %lsn.callshortnames{cidx,1} = [lsn.shortnames{aidx},'->',lsn.shortnames{target_eidx}];
            callDist = Geometric(1/tasks{t}.activities(a).asyncCallMeans(s)); % asynch
            lsn.callproc{cidx,1} = callDist;
            [lsn.callproc_type(cidx), lsn.callproc_params{cidx}, lsn.callproc_mean(cidx), lsn.callproc_scv(cidx), lsn.callproc_proc{cidx}] = extractDistParams(callDist);
            lsn.callsof{aidx}(end+1) = cidx;
            lsn.iscaller(aidx, target_tidx) = true;
            lsn.iscaller(aidx, target_eidx) = true;
            lsn.iscaller(tidx, target_tidx) = true;
            lsn.iscaller(tidx, target_eidx) = true;
            lsn.isasynccaller(tidx, target_tidx) = true;
            lsn.isasynccaller(tidx, target_eidx) = true;
            lsn.isasynccaller(aidx, target_tidx) = true;
            lsn.isasynccaller(aidx, target_eidx) = true;
            lsn.taskgraph(tidx, target_tidx) = 1;
            lsn.graph(aidx, target_eidx) = 1;
        end
    end

    for ap=1:length(tasks{t}.precedences)
        pretype = tasks{t}.precedences(ap).preType;
        posttype = tasks{t}.precedences(ap).postType;
        preacts = tasks{t}.precedences(ap).preActs;
        postacts = tasks{t}.precedences(ap).postActs;

        % Validate PRE_AND activities exist before processing
        if pretype == ActivityPrecedenceType.PRE_AND
            if isempty(preacts)
                line_error(mfilename, sprintf('PRE_AND precedence in task "%s" has no pre activities.', tasks{t}.name));
            end
            for prea = 1:length(preacts)
                preaidx = findstring(lsn.hashnames, ['A:',preacts{prea}]);
                if preaidx <= 0
                    line_error(mfilename, sprintf('PRE_AND precedence references non-existent activity "%s" in task "%s".', preacts{prea}, tasks{t}.name));
                end
                if preaidx > 0 && lsn.parent(preaidx) ~= tidx
                    line_error(mfilename, sprintf('PRE_AND precedence in task "%s" references activity "%s" from a different task.', tasks{t}.name, preacts{prea}));
                end
            end
        end

        % Validate POST_AND activities exist before processing
        if posttype == ActivityPrecedenceType.POST_AND
            if isempty(postacts)
                line_error(mfilename, sprintf('POST_AND precedence in task "%s" has no post activities.', tasks{t}.name));
            end
            for posta = 1:length(postacts)
                postaidx = findstring(lsn.hashnames, ['A:',postacts{posta}]);
                if postaidx <= 0
                    line_error(mfilename, sprintf('POST_AND precedence references non-existent activity "%s" in task "%s".', postacts{posta}, tasks{t}.name));
                end
                if postaidx > 0 && lsn.parent(postaidx) ~= tidx
                    line_error(mfilename, sprintf('POST_AND precedence in task "%s" references activity "%s" from a different task.', tasks{t}.name, postacts{posta}));
                end
            end
        end

        for prea = 1:length(preacts)
            preaidx = findstring(lsn.hashnames, ['A:',tasks{t}.precedences(ap).preActs{prea}]);
            switch pretype
                case ActivityPrecedenceType.PRE_AND
                    quorum = tasks{t}.precedences(ap).preParams;
                    if isempty(quorum)
                        preParam = 1.0;
                    else
                        preParam = quorum / length(preacts);
                    end
                otherwise
                    preParam = 1.0;
            end
            switch posttype
                case ActivityPrecedenceType.POST_OR
                    for posta = 1:length(postacts)
                        postaidx = findstring(lsn.hashnames, ['A:',tasks{t}.precedences(ap).postActs{posta}]);
                        probs = tasks{t}.precedences(ap).postParams;
                        postParam = probs(posta);
                        lsn.graph(preaidx, postaidx) = preParam * postParam;
                        lsn.actpretype(preaidx) = sparse(tasks{t}.precedences(ap).preType);
                        lsn.actposttype(postaidx) = sparse(tasks{t}.precedences(ap).postType);
                    end
                case ActivityPrecedenceType.POST_AND
                    for posta = 1:length(postacts)
                        postaidx = findstring(lsn.hashnames, ['A:',tasks{t}.precedences(ap).postActs{posta}]);
                        lsn.graph(preaidx, postaidx) = 1;
                        lsn.actpretype(preaidx) = sparse(tasks{t}.precedences(ap).preType);
                        lsn.actposttype(postaidx) = sparse(tasks{t}.precedences(ap).postType);
                    end
                case ActivityPrecedenceType.POST_LOOP
                    counts = tasks{t}.precedences(ap).postParams;
                    % add the end activity
                    enda = length(postacts);
                    loopentryaidx = findstring(lsn.hashnames, ['A:',tasks{t}.precedences(ap).preActs{1}]);
                    loopstartaidx = findstring(lsn.hashnames, ['A:',tasks{t}.precedences(ap).postActs{1}]);
                    loopendaidx = findstring(lsn.hashnames, ['A:',tasks{t}.precedences(ap).postActs{enda}]);

                    if counts < 1
                        % When expected iterations < 1, we may skip loop entirely
                        % E[iterations] = counts means P(enter loop) = counts
                        lsn.graph(loopentryaidx, loopstartaidx) = counts;
                        lsn.graph(loopentryaidx, loopendaidx) = 1.0 - counts;
                        % Process activities inside the loop as serial chain
                        curaidx = loopstartaidx;
                        for posta = 2:(length(postacts)-1)
                            postaidx = findstring(lsn.hashnames, ['A:',tasks{t}.precedences(ap).postActs{posta}]);
                            lsn.graph(curaidx, postaidx) = 1.0;
                            lsn.actposttype(postaidx) = sparse((tasks{t}.precedences(ap).postType));
                            curaidx = postaidx;
                        end
                        % After loop body, always exit to end (no looping back)
                        lsn.graph(curaidx, loopendaidx) = 1.0;
                        lsn.actposttype(loopstartaidx) = sparse((tasks{t}.precedences(ap).postType));
                    else
                        % When expected iterations >= 1, always enter loop
                        % E[iterations] = 1/(1-p) = counts => p = 1 - 1/counts
                        curaidx = loopentryaidx;
                        for posta = 1:(length(postacts)-1) % last one is 'end' of loop activity
                            postaidx = findstring(lsn.hashnames, ['A:',tasks{t}.precedences(ap).postActs{posta}]);
                            lsn.graph(curaidx, postaidx) = 1.0;
                            lsn.actposttype(postaidx) = sparse((tasks{t}.precedences(ap).postType));
                            curaidx = postaidx;
                        end
                        loop_back_edges(curaidx, loopstartaidx) = true;
                        lsn.graph(curaidx, loopstartaidx) = 1.0 - 1.0 / counts;
                        lsn.graph(curaidx, loopendaidx) = 1.0 / counts;
                    end
                    lsn.actposttype(loopendaidx) = sparse((tasks{t}.precedences(ap).postType));
                otherwise
                    for posta = 1:length(postacts)
                        postaidx = findstring(lsn.hashnames, ['A:',tasks{t}.precedences(ap).postActs{posta}]);
                        postParam = 1.0;
                        lsn.graph(preaidx, postaidx) = preParam * postParam;
                        lsn.actpretype(preaidx) = sparse(tasks{t}.precedences(ap).preType);
                        lsn.actposttype(postaidx) = sparse(tasks{t}.precedences(ap).postType);
                    end
            end
        end
    end
end

%% Process forwarding calls from entries
for e = 1:length(self.entries)
    entry = self.entries{e};
    eidx = findstring(lsn.hashnames, ['E:', entry.name]);
    if eidx <= 0
        eidx = findstring(lsn.hashnames, ['I:', entry.name]);
    end
    if eidx <= 0
        continue;
    end
    source_tidx = lsn.parent(eidx);

    for fw = 1:length(entry.forwardingDests)
        target_entry_name = entry.forwardingDests{fw};
        target_eidx = findstring(lsn.hashnames, ['E:', target_entry_name]);
        if target_eidx <= 0
            target_eidx = findstring(lsn.hashnames, ['I:', target_entry_name]);
        end
        if target_eidx <= 0
            line_error(mfilename, sprintf('Entry "%s" forwards to non-existent entry "%s".', entry.name, target_entry_name));
        end
        target_tidx = lsn.parent(target_eidx);

        % Validate: cannot forward to same task
        if source_tidx == target_tidx
            line_error(mfilename, sprintf('Entry "%s" cannot forward to entry "%s" on the same task.', entry.name, target_entry_name));
        end

        cidx = cidx + 1;
        lsn.calltype(cidx, 1) = CallType.FWD;
        lsn.callpair(cidx, 1:2) = [eidx, target_eidx];
        lsn.callnames{cidx, 1} = [lsn.names{eidx}, '~>', lsn.names{target_eidx}];
        lsn.callhashnames{cidx, 1} = [lsn.hashnames{eidx}, '~>', lsn.hashnames{target_eidx}];

        % Forwarding probability (stored as mean calls)
        fwdProb = entry.forwardingProbs(fw);
        callDist = Geometric(1.0 / fwdProb);
        lsn.callproc{cidx, 1} = callDist;
        [lsn.callproc_type(cidx), lsn.callproc_params{cidx}, lsn.callproc_mean(cidx), lsn.callproc_scv(cidx), lsn.callproc_proc{cidx}] = extractDistParams(callDist);

        % Update task graph to reflect forwarding relationship
        lsn.taskgraph(source_tidx, target_tidx) = 1;
        lsn.graph(eidx, target_eidx) = 1;

        % NOTE: We do NOT set issynccaller for forwarding calls because:
        % - Forwarding is NOT a sync call - the forwarder does NOT wait for the target
        % - The forwarder sends the request to the target and the target replies
        %   directly to the original caller (not back to the forwarder)
        % - The forwarding target's workload comes from the forwarding probability,
        %   not from a separate closed class in the queueing network decomposition
        %
        % Setting issynccaller here was causing incorrect throughput calculations
        % because it created separate closed classes for forwarding targets.
    end
end

% Check for entries without boundTo activities
unbound = find(all(~lsn.graph(lsn.eshift+1 : lsn.eshift+lsn.nentries, ...
    lsn.ashift+1 : lsn.ashift+lsn.nacts), 2)); %#ok<EFIND>
if ~isempty(unbound)
    line_error(mfilename, 'An entry does not have any boundTo activity.');
end

%lsn.replies = false(1,lsn.nacts);
%lsn.replygraph = 0*lsn.graph;
% Snapshot which entries have explicit replies (from repliesTo calls)
% before adding implicit ones. This way OrFork branches all get marked.
hasExplicitReply = any(lsn.replygraph, 1);
for t = 1:lsn.ntasks
    tidx = lsn.tshift+t;
    for aidx = lsn.actsof{tidx}
        postaidxs = find(lsn.graph(aidx, :));
        isreply = true;
        % if no successor is an action of tidx
        for postaidx = postaidxs
            if any(lsn.actsof{tidx} == postaidx)
                isreply = false;
            end
        end
        if isreply
            % this is a leaf node, search backward for the parent entry,
            % which is assumed to be unique
            %lsn.replies(aidx-lsn.nacts) = true;
            parentidx = aidx;
            while lsn.type(parentidx) ~= LayeredNetworkElement.ENTRY
                ancestors = find(lsn.graph(:,parentidx));
                parentidx = at(ancestors,1); % only choose first ancestor
            end
            if lsn.type(parentidx) == LayeredNetworkElement.ENTRY
                eidx = parentidx - lsn.eshift;
                % Only add implicit reply if no explicit reply exists for this entry
                % This supports Phase-2 activities (activities after an explicit reply)
                if ~hasExplicitReply(eidx)
                    lsn.replygraph(aidx-lsn.ashift, eidx) = true;
                end
            end
        end
    end
end
lsn.ncalls = size(lsn.calltype,1);

% correct multiplicity for infinite server stations
for tidx = find(lsn.sched == SchedStrategy.INF)
    if lsn.type(tidx) == LayeredNetworkElement.TASK
        callers = find(lsn.taskgraph(:, tidx));
        callers_inf = strcmp(lsn.mult(callers), SchedStrategy.INF);
        if any(callers_inf)
            % if a caller is also inf, then we would need to recursively
            % determine the maximum multiplicity, we instead use a
            % heuristic value
            lsn.mult(tidx) = sum(lsn.mult(~callers_inf)) + sum(callers_inf)*max(lsn.mult);
        else
            lsn.mult(tidx) = sum(lsn.mult(callers));
        end
    end
end

lsn.refset = zeros(lsn.nidx,1);
[conncomps, roots]=graph_connected_components(lsn.taskgraph(lsn.nhosts+1:end, lsn.nhosts+1:end));
lsn.conntasks = conncomps;
for r=1:length(roots)
    lsn.conntasks(find(lsn.conntasks == r)) = lsn.tshift+roots(r);
end

lsn.isref = lsn.sched == SchedStrategy.REF;
lsn.iscache(1:(lsn.tshift+lsn.ntasks)) = lsn.nitems(1:(lsn.tshift+lsn.ntasks))>0;
lsn.isfunction(1:(lsn.tshift+lsn.ntasks)) = ~cellfun(@isempty, lsn.setuptime(1:(lsn.tshift+lsn.ntasks)));

% the dag differs from the graph:
% - dag swaps the direction of entry-task edges
% - dag removes loop edges
dag = lsn.graph;
n = size(dag, 1);
% Reverse edges from TASK to ENTRY if not a reference
for i = 1:n
    if lsn.type(i) == LayeredNetworkElement.TASK && ~lsn.isref(i)
        for j = 1:n
            if lsn.type(j) == LayeredNetworkElement.ENTRY && dag(i,j)
                dag(i,j) = 0;
                dag(j,i) = 1;
            end
        end
    end
end

% Compute entry-to-activity reachability within the same task
for eoff = 1:lsn.nentries
    eidx = lsn.eshift + eoff;          % global entry index
    tidx = lsn.parent(eidx);           % parent task index
    visited = false(1,nidx);           % global visit mask
    stack   = eidx;
    visited(eidx) = true;
    while ~isempty(stack)
        v        = stack(end);
        stack(end) = [];
        nbrs     = find(lsn.graph(v,:));
        nbrs     = nbrs(~visited(nbrs));
        visited(nbrs) = true;
        stack    = [stack nbrs];       %#ok<AGROW>
    end
    acts = find(visited & ...
        lsn.type'  == LayeredNetworkElement.ACTIVITY & ...
        lsn.parent == tidx);
    lsn.actsof{lsn.eshift+eoff} = acts;
end

%% check for errors
dag(loop_back_edges(:)) = 0;
lsn.dag = dag;
% Compute bounds on multiplicies for host processors and non-ref tasks
if is_dag(dag)
    lsn.maxmult = lsn_max_multiplicity(lsn);
    lsn.maxmult = lsn.maxmult(1:(lsn.tshift+lsn.ntasks));
else
    line_error(mfilename, 'A cycle exists in an activity graph.');
end
% Check for non-terminal reply activities
% An activity that replies to an entry should not have Phase 1 successor activities
% Phase 2 successors are allowed (post-reply processing)
for a = 1:lsn.nacts
    if any(lsn.replygraph(a, :))  % activity 'a' replies to some entry
        aidx = lsn.ashift + a;     % global activity index
        successors = find(lsn.graph(aidx, :));
        for succ = successors
            if succ > lsn.eshift + lsn.nentries  % successor is an activity
                succ_act_idx = succ - lsn.ashift;  % convert to activity array index
                if lsn.actphase(succ_act_idx) == 1
                    line_error(mfilename, 'Unsupported replyTo in non-terminal activity.');
                end
            end
        end
    end
end

% Check calls
if ~isempty(lsn.callpair)
    target_eidxs = unique(lsn.callpair(:,2));
    for eidx=target_eidxs(:)'
        call_types_to_eidx = lsn.calltype(find(lsn.callpair(:,2) == eidx),1);
        if ~all(call_types_to_eidx == call_types_to_eidx(1))
            line_error(mfilename, 'An entry is called both synchronously and asynchronously.');
        end
    end
end

self.lsn = lsn;
end

function [dtype, params, mean_val, scv_val, proc] = extractDistParams(dist)
% EXTRACTDISTPARAMS Extract primitive parameters from a Distribution object
%
% Returns:
%   dtype    - ProcessType enum value
%   params   - Vector/cell of primitive parameters
%   mean_val - Precomputed mean (NaN if unavailable)
%   scv_val  - Precomputed SCV (NaN if unavailable)
%   proc     - Process representation {D0, D1} or [] for non-Markovian

if isempty(dist)
    dtype = ProcessType.DISABLED;
    params = [];
    mean_val = NaN;
    scv_val = NaN;
    proc = {};
    return
end

% Get mean and SCV
try
    mean_val = dist.getMean();
catch
    mean_val = NaN;
end

try
    scv_val = dist.getSCV();
catch
    scv_val = NaN;
end

% Get process representation if available
try
    proc = dist.getRepres();
    if isempty(proc)
        proc = {};
    end
catch
    proc = {};
end

% Extract parameters based on distribution type
distClass = class(dist);
switch distClass
    case 'Disabled'
        dtype = ProcessType.DISABLED;
        params = [];

    case 'Immediate'
        dtype = ProcessType.IMMEDIATE;
        params = [];
        mean_val = 0;
        scv_val = 0;

    case 'Exp'
        dtype = ProcessType.EXP;
        params = dist.getParam(1).paramValue;  % lambda

    case 'Erlang'
        dtype = ProcessType.ERLANG;
        params = [dist.getParam(1).paramValue; dist.getParam(2).paramValue];  % [k; mu]

    case 'HyperExp'
        dtype = ProcessType.HYPEREXP;
        % p1, lambda1, lambda2
        params = [dist.getParam(1).paramValue; dist.getParam(2).paramValue; dist.getParam(3).paramValue];

    case 'Coxian'
        dtype = ProcessType.COXIAN;
        mu = dist.getParam(1).paramValue;
        phi = dist.getParam(2).paramValue;
        params = [length(mu); mu(:); phi(:)];

    case 'Cox2'
        dtype = ProcessType.COX2;
        mu1 = dist.getParam(1).paramValue;
        mu2 = dist.getParam(2).paramValue;
        phi = dist.getParam(3).paramValue;
        params = [mu1; mu2; phi];

    case 'APH'
        dtype = ProcessType.APH;
        alpha = dist.getParam(1).paramValue;
        T = dist.getParam(2).paramValue;
        n = length(alpha);
        params = [n; alpha(:); T(:)];

    case 'PH'
        dtype = ProcessType.PH;
        alpha = dist.getParam(1).paramValue;
        T = dist.getParam(2).paramValue;
        n = length(alpha);
        params = [n; alpha(:); T(:)];

    case 'Det'
        dtype = ProcessType.DET;
        params = dist.getParam(1).paramValue;

    case 'Uniform'
        dtype = ProcessType.UNIFORM;
        params = [dist.getParam(1).paramValue; dist.getParam(2).paramValue];  % [a; b]

    case 'Gamma'
        dtype = ProcessType.GAMMA;
        params = [dist.getParam(1).paramValue; dist.getParam(2).paramValue];  % [shape; scale]

    case 'Pareto'
        dtype = ProcessType.PARETO;
        params = [dist.getParam(1).paramValue; dist.getParam(2).paramValue];  % [shape; scale]

    case 'Weibull'
        dtype = ProcessType.WEIBULL;
        params = [dist.getParam(1).paramValue; dist.getParam(2).paramValue];  % [shape; scale]

    case 'Lognormal'
        dtype = ProcessType.LOGNORMAL;
        params = [dist.getParam(1).paramValue; dist.getParam(2).paramValue];  % [mu; sigma]

    case 'MAP'
        dtype = ProcessType.MAP;
        D0 = dist.getParam(1).paramValue;
        D1 = dist.getParam(2).paramValue;
        n = size(D0,1);
        params = [n; D0(:); D1(:)];

    case 'MMPP2'
        dtype = ProcessType.MMPP2;
        lambda0 = dist.getParam(1).paramValue;
        lambda1 = dist.getParam(2).paramValue;
        sigma0 = dist.getParam(3).paramValue;
        sigma1 = dist.getParam(4).paramValue;
        params = [lambda0; lambda1; sigma0; sigma1];

    case 'Geometric'
        dtype = ProcessType.GEOMETRIC;
        params = dist.getParam(1).paramValue;  % p (success probability)

    case 'Poisson'
        dtype = ProcessType.POISSON;
        params = dist.getParam(1).paramValue;  % lambda

    case 'Binomial'
        dtype = ProcessType.BINOMIAL;
        params = [dist.getParam(1).paramValue; dist.getParam(2).paramValue];  % [n; p]

    case 'Bernoulli'
        dtype = ProcessType.BERNOULLI;
        params = dist.getParam(1).paramValue;  % p

    case {'Replayer', 'Trace'}
        dtype = ProcessType.REPLAYER;
        params = [];  % Trace data not stored as primitive

    case 'DiscreteSampler'
        dtype = ProcessType.GEOMETRIC;  % Approximation
        params = [];

    otherwise
        % Fallback for unhandled types
        dtype = ProcessType.DISABLED;
        params = [];
        proc = {};
end
end