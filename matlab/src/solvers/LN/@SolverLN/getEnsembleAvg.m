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

iterate(self); % run iterations
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
    serverIdx = self.ensemble{e}.attribute.serverIdx;
    sourceIdx = self.ensemble{e}.attribute.sourceIdx;
    % determine processor metrics
    if self.ensemble{e}.stations{serverIdx}.attribute.ishost
        hidx = self.ensemble{e}.stations{serverIdx}.attribute.idx;
        TN(hidx) = 0;
        PN(hidx) = 0;
        for c=1:self.ensemble{e}.getNumberOfClasses
            if self.ensemble{e}.classes{c}.completes
                t = 0;
                u = 0;
                if ~isnan(clientIdx)
                    t = max(t, self.results{end,e}.TN(clientIdx,c));
                end
                if ~isnan(sourceIdx)
                    t = max(t, self.results{end,e}.TN(sourceIdx,c));
                end
                TN(hidx) = TN(hidx) + max(t,self.results{end,e}.TN(serverIdx,c));
            end
            type = self.ensemble{e}.classes{c}.attribute(1);
            switch type
                case LayeredNetworkElement.ACTIVITY
                    aidx = self.ensemble{e}.classes{c}.attribute(2);
                    tidx = self.lqn.parent(aidx);
                    if isnan(PN(aidx)), PN(aidx)=0; end
                    if isnan(PN(tidx)), PN(tidx)=0; end
                    PN(aidx) = PN(aidx) + self.results{end,e}.UN(serverIdx,c);
                    PN(tidx) = PN(tidx) + self.results{end,e}.UN(serverIdx,c);
                    PN(hidx) = PN(hidx) + self.results{end,e}.UN(serverIdx,c);
            end
        end
        TN(hidx) = NaN; % added for consistency with LQNS
    end

    % determine remaining metrics
    for c=1:self.ensemble{e}.getNumberOfClasses
        type = self.ensemble{e}.classes{c}.attribute(1);
        switch type
            case LayeredNetworkElement.TASK
                tidx = self.ensemble{e}.classes{c}.attribute(2);
                if self.ensemble{e}.stations{serverIdx}.attribute.ishost
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
                if self.ensemble{e}.stations{serverIdx}.attribute.ishost
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

    for aidx=self.lqn.actsof{tidx}
        UN(aidx) = TN(aidx)*SN(aidx);
    end
    UN(tidx) = UN(tidx) + UN(eidx);
end

for idx=find(self.ignore)
    QN(idx)=0;
    UN(idx)=0;
    RN(idx)=0;
    TN(idx)=0;
    PN(idx)=0;
    SN(idx)=0;
    WN(idx)=0;
    AN(idx)=0;
end

QN = UN;
UN = PN;
RN = SN;
end
