function stationStateAggr = sampleAggr(self, node, numEvents, markActivePassive)
% STATIONSTATEAGGR = SAMPLEAGGR(NODE, NUMEVENTS)

if GlobalConstants.DummyMode
    stationStateAggr = NaN;
    return
end

if nargin<2 %~exist('node','var')
    line_error(mfilename,'sampleAggr requires to specify a node.');
end

if nargin<4
    markActivePassive = false;
end

if nargin<3 %~exist('numEvents','var')
    %numEvents = -1;
else
    line_warning(mfilename,'JMT does not allow to fix the number of events for individual nodes. The number of returned events may be inaccurate.\n');
    %numEvents = numEvents - 1; % we include the initialization as an event
end
sn = self.getStruct;

Q = getAvgQLenHandles(self);
% create a temp model
modelCopy = self.model.copy;
modelCopy.resetNetwork;

% determine the nodes to logs
isNodeClassLogged = false(modelCopy.getNumberOfNodes, modelCopy.getNumberOfClasses);
ind = self.model.getNodeIndex(node.name);
for r=1:modelCopy.getNumberOfClasses
    isNodeClassLogged(ind,r) = true;
end
% apply logging to the copied model
Plinked = sn.rtorig;
isNodeLogged = max(isNodeClassLogged,[],2);
logpath = lineTempDir;
modelCopy.linkAndLog(Plinked, isNodeLogged, logpath);
% simulate the model copy and retrieve log data
solverjmt = SolverJMT(modelCopy, self.getOptions);
if nargin>=3 && numEvents > 0
    solverjmt.maxEvents = numEvents*sn.nnodes*sn.nclasses;
else
    solverjmt.maxEvents = -1;
    numEvents = self.getOptions.samples;
end
solverjmt.getAvg(); % log data
logData = SolverJMT.parseLogs(modelCopy, isNodeLogged, MetricType.toText(MetricType.QLen));

% from here convert from nodes in logData to stations
sn = modelCopy.getStruct;
ind = self.model.getNodeIndex(node.getName);
isf = sn.nodeToStateful(ind);
t = [];
nir = cell(1,sn.nclasses);
event = cell(1,sn.nclasses);
%ids = cell(1,sn.nclasses);

for r=1:sn.nclasses
    if isempty(logData{ind,r})
        nir{r} = [];
    else
        [~,uniqTS] = unique(logData{ind,r}.t);
        if isNodeClassLogged(isf,r)
            if ~isempty(logData{ind,r})
                t = logData{ind,r}.t(uniqTS);
                t = [t(2:end);t(end)];
                nir{r} = logData{ind,r}.QLen(uniqTS);
                event{r} = logData{ind,r}.event;
                %ids{r} = logData{ind,r}.
            end
        end
    end
end
if isfinite(self.options.timespan(2))
    stopAt = find(t>self.options.timespan(2),1,'first');
    if ~isempty(stopAt) && stopAt>1
        t = t(1:(stopAt-1));
        for r=1:length(nir)
            nir{r} = nir{r}(1:(stopAt-1));
        end
    end
end

if length(t) < 1+numEvents
    line_warning(mfilename,'LINE could not estimate correctly the JMT simulation length to return the desired number of events at the specified node. Try to re-run increasing the number of events.\n');
end

stationStateAggr = struct();
stationStateAggr.handle = node;
stationStateAggr.t = t;
stationStateAggr.t = stationStateAggr.t(1:min(length(t),1+numEvents),:);
stationStateAggr.t = [0; stationStateAggr.t(1:end-1)];
stationStateAggr.state = cell2mat(nir);
stationStateAggr.state = stationStateAggr.state(1:min(length(t),1+numEvents),:);
%stationStateAggr.job_id =

event = cellmerge(event);
event = {event{cellisa(event,'Event')}}';
event_t = cellfun(@(c) c.t, event);
event_t = event_t(event_t <= max(stationStateAggr.t));
[~,I]=sort(event_t);
stationStateAggr.event = {event{I}};
stationStateAggr.event = stationStateAggr.event';
stationStateAggr.isaggregate = true;

if markActivePassive
    apevent = cell(1,length(stationStateAggr.t)-1);
    for ti = 1:length(apevent)
        apevent{ti} = struct('active',[],'passive',[]);
    end
    for e=1:length(stationStateAggr.event)
        ti = find(stationStateAggr.event{e}.t == stationStateAggr.t);
        if ~isempty(ti)
            switch stationStateAggr.event{e}.event
                case EventType.ARV
                    apevent{ti-1}.passive = stationStateAggr.event{e};
                otherwise
                    apevent{ti-1}.active = stationStateAggr.event{e};
            end
        end
    end
    stationStateAggr.event = apevent';
end
end