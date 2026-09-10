function [QN,UN,RN,TN,AN,WN] = getEnsembleAvg(self)
% [QN,UN,RN,TN,AN,WN] = GETENSEMBLEAVG(SELF)

% Check if solver was properly constructed (may have returned early due to unsupported features)
if isempty(self.ensemble)
    QN = []; UN = []; RN = []; TN = []; AN = []; WN = [];
    return;
end

% Show library attribution if verbose and not yet shown
if self.options.verbose ~= VerboseLevel.SILENT && ~GlobalConstants.isLibraryAttributionShown()
    libs = SolverLN.getLibrariesUsed([], self.options);
    if ~isempty(libs)
        line_printf('The solver will leverage %s.\n', strjoin(libs, ', '));
        GlobalConstants.setLibraryAttributionShown(true);
    end
end

if self.isPHEncoding()
    % the layers of a PH encoding carry one class per caller task, so the
    % per-element results are rebuilt analytically -- see getEnsembleAvgPH
    [QN,UN,RN,TN,AN,WN] = getEnsembleAvgPH(self);
    return
end

% Solver console: SolverLN is an ensemble solver and does not pass through
% NetworkSolver.runAnalyzerChecks, so it opens its own run here. The guard
% must live until this function returns.
consoleGuard = LineConsole.beginRun(self, self.options); %#ok<NASGU>
LineConsole.loop('solving the layered fixed point over %d layers', self.nlayers);
lnRuntime = iterate(self); % run iterations
QN  = nan(self.lqn.nidx,1);
UN  = nan(self.lqn.nidx,1);
RN  = nan(self.lqn.nidx,1);
TN  = nan(self.lqn.nidx,1);
PN  = nan(self.lqn.nidx,1); % utilization will be first stored here
SN  = nan(self.lqn.nidx,1); % response time will be first stored here
WN  = nan(self.lqn.nidx,1); % residence time
AN  = nan(self.lqn.nidx,1); % not available yet
WN_processed = false(self.lqn.nidx,1); % track activities already accumulated into task WN
E = self.nlayers;
for e=1:E
    clientIdx = self.ensemble{e}.attribute.clientIdx;
    sourceIdx = self.ensemble{e}.attribute.sourceIdx;
    hostStations = self.serverStationsOf(e, true);
    hasHostServer = ~isempty(hostStations);
    % determine processor metrics, one processor at a time under flat layering
    for hs = hostStations
        hidx = self.ensemble{e}.stations{hs}.attribute.idx;
        TN(hidx) = 0;
        PN(hidx) = 0;
        for c=1:self.ensemble{e}.getNumberOfClasses
            if self.ensemble{e}.classes{c}.completes
                t = 0;
                if ~isnan(clientIdx)
                    t = max(t, self.results{end,e}.TN(clientIdx,c));
                end
                if ~isnan(sourceIdx)
                    t = max(t, self.results{end,e}.TN(sourceIdx,c));
                end
                TN(hidx) = TN(hidx) + max(t,self.results{end,e}.TN(hs,c));
            end
            type = self.ensemble{e}.classes{c}.attribute(1);
            switch type
                case LayeredNetworkElement.ACTIVITY
                    if self.stationIdxOfClass(e,c) ~= hs
                        continue % the activity does not run on this processor
                    end
                    aidx = self.ensemble{e}.classes{c}.attribute(2);
                    tidx = self.lqn.parent(aidx);
                    if isnan(PN(aidx)), PN(aidx)=0; end
                    if isnan(PN(tidx)), PN(tidx)=0; end
                    PN(aidx) = PN(aidx) + self.results{end,e}.UN(hs,c);
                    PN(tidx) = PN(tidx) + self.results{end,e}.UN(hs,c);
                    PN(hidx) = PN(hidx) + self.results{end,e}.UN(hs,c);
            end
        end
        TN(hidx) = NaN; % added for consistency with LQNS
    end

    % determine remaining metrics
    for c=1:self.ensemble{e}.getNumberOfClasses
        type = self.ensemble{e}.classes{c}.attribute(1);
        serverIdx = self.stationIdxOfClass(e,c);
        switch type
            case LayeredNetworkElement.TASK
                tidx = self.ensemble{e}.classes{c}.attribute(2);
                if hasHostServer
                    if isnan(TN(tidx))
                        % store the result in the processor
                        % model
                        TN(tidx) = self.results{end,e}.TN(clientIdx,c);
                    end
                else
                    % nop
                end
            case LayeredNetworkElement.ENTRY
                eidx = self.ensemble{e}.classes{c}.attribute(2);
                tidx = self.lqn.parent(eidx);
                % For phase-2 models, use residt (caller's view with overtaking)
                % Otherwise use servt (total service time = response time)
                if self.hasPhase2 && self.servt_ph2(eidx) > GlobalConstants.FineTol
                    SN(eidx) = self.residt(eidx);  % Phase-1 + overtaking correction
                else
                    SN(eidx) = self.servt(eidx);
                end
                if hasHostServer
                    if isnan(TN(eidx))
                        % store the result in the processor model
                        if isnan(TN(eidx)), TN(eidx)=0; end
                        TN(eidx) = self.results{end,e}.TN(clientIdx,c);
                    end
                else
                    % nop
                end
            case LayeredNetworkElement.CALL
                cidx = self.ensemble{e}.classes{c}.attribute(2);
                aidx = self.lqn.callpair(cidx,1);
                % Only sync calls contribute to caller's response time
                if self.lqn.calltype(cidx) == CallType.SYNC
                    SN(aidx) = SN(aidx) + self.results{end,e}.RN(serverIdx,c) * self.lqn.callproc{cidx}.getMean();
                end
                if isnan(QN(aidx)), QN(aidx)=0; end
                QN(aidx) = QN(aidx) + self.results{end,e}.QN(serverIdx,c);
            case LayeredNetworkElement.ACTIVITY
                aidx = self.ensemble{e}.classes{c}.attribute(2);
                tidx = self.lqn.parent(aidx);
                if isnan(QN(tidx)), QN(tidx)=0; end
                QN(tidx) = QN(tidx) + self.results{end,e}.QN(serverIdx,c);
                if isnan(TN(aidx)), TN(aidx)=0; end
                if isnan(QN(aidx)), QN(aidx)=0; end

                % For forwarding targets: propagate activity metrics to entry and task
                % Check if task has its own class (non-forwarding targets do)
                hasTaskClass = any(self.ensemble{e}.attribute.tasks(:,2) == tidx);
                if ~hasTaskClass
                    % Propagate activity throughput to task for forwarding targets
                    if isnan(TN(tidx)), TN(tidx)=0; end
                    switch self.ensemble{e}.classes{c}.type
                        case JobClassType.CLOSED
                            TN(tidx) = TN(tidx) + self.results{end,e}.TN(serverIdx,c);
                        case JobClassType.OPEN
                            TN(tidx) = TN(tidx) + self.results{end,e}.TN(sourceIdx,c);
                    end
                end

                % Find the entry this activity is bound to
                for eidx_check = self.lqn.entriesof{tidx}
                    if self.lqn.graph(eidx_check, aidx) > 0  % Activity is bound to this entry
                        if isnan(TN(eidx_check)), TN(eidx_check)=0; end
                        if isnan(QN(eidx_check)), QN(eidx_check)=0; end
                        if isnan(SN(eidx_check)), SN(eidx_check)=0; end
                        % Add activity metrics to entry (for entries without their own classes)
                        actTput = 0;
                        switch self.ensemble{e}.classes{c}.type
                            case JobClassType.CLOSED
                                actTput = self.results{end,e}.TN(serverIdx,c);
                            case JobClassType.OPEN
                                actTput = self.results{end,e}.TN(sourceIdx,c);
                        end
                        % Only add if entry doesn't have its own class (forwarding target case)
                        hasEntryClass = any(self.ensemble{e}.attribute.entries(:,2) == eidx_check);
                        if ~hasEntryClass
                            TN(eidx_check) = TN(eidx_check) + actTput;
                            QN(eidx_check) = QN(eidx_check) + self.results{end,e}.QN(serverIdx,c);
                            SN(eidx_check) = SN(eidx_check) + self.results{end,e}.RN(serverIdx,c);
                        end
                        break;
                    end
                end
                switch self.ensemble{e}.classes{c}.type
                    case JobClassType.CLOSED
                        TN(aidx) = TN(aidx) + self.results{end,e}.TN(serverIdx,c);
                    case JobClassType.OPEN
                        TN(aidx) = TN(aidx) + self.results{end,e}.TN(sourceIdx,c);
                end
                %                            SN(aidx) = self.servt(aidx);
                if isnan(SN(aidx)), SN(aidx)=0; end
                SN(aidx) = SN(aidx) + self.results{end,e}.RN(serverIdx,c);
                if isnan(RN(aidx)), RN(aidx)=0; end
                RN(aidx) = RN(aidx) + self.results{end,e}.RN(serverIdx,c);
                if isnan(WN(aidx)), WN(aidx)=0; end
                if isnan(WN(tidx)), WN(tidx)=0; end
                % Use self.residt (computed via QN/TN_ref in updateMetricsDefault)
                % instead of layer WN to avoid fork+loop visit distortion
                WN(aidx) = self.residt(aidx);
                if ~WN_processed(aidx)
                    WN(tidx) = WN(tidx) + self.residt(aidx);
                    WN_processed(aidx) = true;
                end
                if isnan(QN(aidx)), QN(aidx)=0; end
                QN(aidx) = QN(aidx) + self.results{end,e}.QN(serverIdx,c);
        end
    end
end

% A replicated element is solved as ONE representative replica, so the rates and
% the busy time read off its layer are that replica's. What the element itself
% delivers is REPL times as much, which is the convention LDES reports and the
% one flow balance across a call needs: on a two-replica probe the caller runs at
% 8.30737 and calls the replicated task once per invocation, so the task's rate
% is 8.30737 and not the 4.14871 of either copy. Response times are per request
% and are left alone. Inert wherever nothing is replicated.
for tidx = (self.lqn.tshift+1):(self.lqn.tshift+self.lqn.ntasks)
    nrep = self.lqn.repl(tidx);
    if ~(nrep > 1)
        continue
    end
    scaleIdx = tidx;
    for eidx = self.lqn.entriesof{tidx}
        scaleIdx(end+1) = eidx; %#ok<AGROW>
    end
    for aidx = self.lqn.actsof{tidx}
        scaleIdx(end+1) = aidx; %#ok<AGROW>
    end
    for idx = scaleIdx
        if ~isnan(TN(idx)), TN(idx) = nrep * TN(idx); end
        if ~isnan(PN(idx)), PN(idx) = nrep * PN(idx); end
    end
end
for hidx = 1:self.lqn.nhosts
    if self.lqn.repl(hidx) > 1 && ~isnan(PN(hidx))
        PN(hidx) = self.lqn.repl(hidx) * PN(hidx);
    end
end

for e=1:self.lqn.nentries
    eidx = self.lqn.eshift + e;
    tidx = self.lqn.parent(eidx);
    if isnan(UN(tidx)), UN(tidx)=0; end

    % Phase-2 support: utilization includes both phases
    if self.hasPhase2 && self.servt_ph2(eidx) > GlobalConstants.FineTol
        % Phase-1 utilization
        self.util_ph1(eidx) = TN(eidx) * self.servt_ph1(eidx);
        % Phase-2 utilization
        self.util_ph2(eidx) = TN(eidx) * self.servt_ph2(eidx);
        % Total utilization = both phases (server is busy during both)
        UN(eidx) = self.util_ph1(eidx) + self.util_ph2(eidx);
    else
        % Standard calculation for entries without phase-2
        UN(eidx) = TN(eidx)*SN(eidx);
    end

    % Entry utilization = sum of activity processor utilizations for that entry
    entryActs = self.lqn.actsof{eidx};
    if ~isempty(entryActs)
        PN(eidx) = sum(PN(entryActs(~isnan(PN(entryActs)))));
    end

    for aidx=self.lqn.actsof{tidx}
        UN(aidx) = TN(aidx)*SN(aidx);
    end
    UN(tidx) = UN(tidx) + UN(eidx);
end

% AN IGNORED ELEMENT IS IDLE, NOT UNDEFINED, and the two are different cells.
% Its component holds no reference task, so nothing reaches it and every
% measure it HAS is zero -- but the measures its kind never has stay NaN,
% exactly as they do for a reachable element. A flat zero over all six columns
% broke the table's NaN mask (a processor with a queue length of 0, an arrival
% rate reported where no solver reports one), and the mask is part of the
% answer: see _kb/06-solver-catalog.md. Reported columns are QLen=UN, Util=PN,
% RespT=SN, ResidT=WN, ArvR=AN, Tput=TN, so the pre-swap QN and RN are
% discarded below and are not written here.
for idx=find(self.ignore)'
    PN(idx) = 0;   % every kind reports a utilization
    AN(idx) = NaN; % nothing reports an arrival rate on an LQN
    switch self.lqn.type(idx)
        case LayeredNetworkElement.PROCESSOR
            UN(idx) = NaN; SN(idx) = NaN; WN(idx) = NaN; TN(idx) = NaN;
        case LayeredNetworkElement.TASK
            UN(idx) = 0;   SN(idx) = NaN; WN(idx) = 0;   TN(idx) = 0;
        case LayeredNetworkElement.ENTRY
            UN(idx) = 0;   SN(idx) = 0;   WN(idx) = NaN; TN(idx) = 0;
        case LayeredNetworkElement.ACTIVITY
            UN(idx) = 0;   SN(idx) = 0;   WN(idx) = 0;   TN(idx) = 0;
    end
end

QN = UN;
UN = PN;
RN = SN;

% Closing banner, the line every NetworkSolver prints. SolverLN is an
% EnsembleSolver and never reaches NetworkSolver.setAvgResults, so it had none.
self.reportCompletion(lnRuntime);
end
