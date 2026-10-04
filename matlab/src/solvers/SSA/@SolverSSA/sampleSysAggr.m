function tranSysState = sampleSysAggr(self, numEvents, markActivePassive)
% TRANSYSSTATEAGGR = sampleSysAggr(NUMSAMPLES)

% The trajectory is the C++ engine's; markActivePassive below is a rearrangement
% of the event cell and applies to it unchanged.
useCpp = isfield(self.options,'lang') && strcmp(self.options.lang,'cpp');

options = self.getOptions;

if GlobalConstants.DummyMode
    tranSysState = NaN;
    return
end

if nargin<2 %~exist('numEvents','var')
    numEvents = options.samples;
end

if nargin<3
    markActivePassive = false;
end

if useCpp
    tranSysState = CPPLINE.sysSamplePath(self.name, self.model, self.options, numEvents, true);
else
switch options.method
    case {'default','serial'}
        options.samples = numEvents;
        options.force = true;
        sn = self.getStruct;
        options.method = 'serial'; % nrm does not support tran*

        [~, tranSystemState, tranSync] = self.runAnalyzer(options);
        tranSysState = struct();
        tranSysState.handle = self.model.getStatefulNodes';
        tranSysState.t = tranSystemState{1};
        tranSysState.state = {tranSystemState{2:end}};
        tranSysState.event = tranSync;
        event = tranSync;

        for isf=1:sn.nstateful
            if size(tranSysState.state{isf},1) > numEvents
                tranSysState.t = tranSystemState(1:numEvents);
                tranSysState.state{isf} = tranSysState.state{isf}(1:numEvents,:);
            end
            [~,tranSysState.state{isf}] = State.toMarginal(sn,sn.statefulToNode(isf),tranSysState.state{isf});
        end

        sn = self.getStruct;
        % Derived tags of each step (see solver_ssa): sn.sync carries none.
        startTag = {};
        preemptTag = {};
        if isfield(self.result,'startTag')
            startTag = self.result.startTag;
            preemptTag = self.result.preemptTag;
        end
        tranSysState.event = {};
        for e = 1:length(event)
            for a=1:length(sn.sync{event(e)}.active)
                tranSysState.event{end+1} = sn.sync{event(e)}.active{a};
                tranSysState.event{end}.t = tranSysState.t(e);
            end
            for p=1:length(sn.sync{event(e)}.passive)
                tranSysState.event{end+1} = sn.sync{event(e)}.passive{p};
                tranSysState.event{end}.t = tranSysState.t(e);
            end
            % PREEMPT before START at the same instant: the victim leaves the
            % server before the job that displaced it takes it.
            if e <= numel(preemptTag) && ~isempty(preemptTag{e})
                for j = 1:size(preemptTag{e},1)
                    tranSysState.event{end+1} = Event(EventType.PREEMPT, ...
                        sn.statefulToNode(preemptTag{e}(j,1)), preemptTag{e}(j,2));
                    tranSysState.event{end}.t = tranSysState.t(e);
                end
            end
            if e <= numel(startTag) && ~isempty(startTag{e})
                for j = 1:size(startTag{e},1)
                    tranSysState.event{end+1} = Event(EventType.START, ...
                        sn.statefulToNode(startTag{e}(j,1)), startTag{e}(j,2));
                    tranSysState.event{end}.t = tranSysState.t(e);
                end
            end
        end
        tranSysState.isaggregate = true;
    otherwise
        line_error(mfilename,'sampleSys is not available in SolverSSA with the chosen method.');
end
end

if markActivePassive
    apevent = cell(1,length(tranSysState.t)-1);
    for ti = 1:length(apevent)
        % START/PREEMPT get their own fields, as in @SolverSSA/sample.m
        apevent{ti} = struct('active',[],'passive',[],'start',{{}},'preempt',{{}});
    end
    for e=1:length(tranSysState.event)
        ti = find(tranSysState.event{e}.t == tranSysState.t);
        if ~isempty(ti) && ti<length(tranSysState.t)
            switch tranSysState.event{e}.event
                case EventType.ARV
                    apevent{ti}.passive = tranSysState.event{e};
                case EventType.START
                    apevent{ti}.start{end+1} = tranSysState.event{e};
                case EventType.PREEMPT
                    apevent{ti}.preempt{end+1} = tranSysState.event{e};
                otherwise
                    apevent{ti}.active = tranSysState.event{e};
            end
        end
    end
    tranSysState.event = apevent';
end
end