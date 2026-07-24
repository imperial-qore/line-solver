function [state, evtype, evclass, evjob] = parseTranState(fileArv, fileDep, nodePreload)
% [STATE, EVTYPE, EVCLASS, EVJOB] = PARSETRANSTATE(FILEARV, FILEDEP, NODEPRELOAD)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

load(fileArv,'jobArvTS','jobArvClassID','jobArvID');
load(fileDep,'jobDepTS','jobDepClassID','jobDepID');

%% compute joint state at station
nClasses = length(nodePreload);
state = [jobArvTS,zeros(length(jobArvTS),nClasses);
    jobDepTS,zeros(length(jobDepTS),nClasses)];

evtype = nan(size(state,1),1);
evclass = zeros(size(state,1),1);
evjob = zeros(size(state,1),1);

for i=1:size(jobArvTS,1)
    state(i,1+jobArvClassID(i)) = +1;
    evtype(i) = EventType.ARV;
    evclass(i) = jobArvClassID(i);
    evjob(i) = jobArvID(i);
end

for i=1:size(jobDepTS,1)
    state(length(jobArvTS)+i,1+jobDepClassID(i)) = -1;
    evtype(length(jobArvTS)+i) = EventType.DEP;
    evclass(length(jobArvTS)+i) = jobDepClassID(i);
    evjob(length(jobArvTS)+i) = jobDepID(i);
end
[state,I] = sortrows(state,1); % sort on timestamps
state = [0,nodePreload;state];
evtype = [EventType.INIT; evtype(I)];
evclass = [NaN;evclass(I)];
evjob = [NaN;evjob(I)];

ev_inst = find(diff(state(:,1)) == 0); % instantaneous events
%evjobs_inst = unique(evjob(ev_inst)); % determine jobs involved in instantaneous events
toSwap = [];
for ev=ev_inst(:)'
    ej = evjob(ev); % find job involved in the instantaneous event;
    prev_ev_ej = find(evjob(1:(ev-1))==ej,1,'last'); % find last event in which this job was involved
    next_ev_ej = ev+find(evjob((ev+1):end)==ej,1); % find next event in which this job was involved
    et = evtype(ev); % find type of instantaneous event
    if ~isempty(prev_ev_ej)
        switch et
            case EventType.ARV
                if evtype(prev_ev_ej) == EventType.ARV
                    % causality error, two successive arrivals, swap position
                    % with next instantaneous event for this job
                    %% swap state
                    tmp = state(ev,:);
                    state(ev,:) = state(next_ev_ej,:);
                    state(next_ev_ej,:) = tmp;
                    %% swap evjob
                    tmp = evjob(ev,:);
                    evjob(ev,:) = evjob(next_ev_ej,:);
                    evjob(next_ev_ej,:) = tmp;
                    %% swap evtype
                    tmp = evtype(ev,:);
                    evtype(ev,:) = evtype(next_ev_ej,:);
                    evtype(next_ev_ej,:) = tmp;
                end
            case EventType.DEP
                if evtype(prev_ev_ej) == EventType.DEP
                    % causality error, two successive departures, swap position
                    % with next instantaneous event for this job
                    %% swap state
                    tmp = state(ev,:);
                    state(ev,:) = state(next_ev_ej,:);
                    state(next_ev_ej,:) = tmp;
                    %% swap evjob
                    tmp = evjob(ev,:);
                    evjob(ev,:) = evjob(next_ev_ej,:);
                    evjob(next_ev_ej,:) = tmp;
                    %% swap evtype
                    tmp = evtype(ev,:);
                    evtype(ev,:) = evtype(next_ev_ej,:);
                    evtype(next_ev_ej,:) = tmp;
                end
        end
    end
end

for j=2:(nClasses+1)
    state(:,j) = cumsum(state(:,j)); %+nodePreload(j-1);
end

end