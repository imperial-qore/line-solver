function [QN,UN,RN,TN,AN,WN] = getEnsembleAvg(self,~,~,~,~, useLQNSnaming)
% [QN,UN,RN,TN,AN,WN] = GETENSEMBLEAVG(SELF,~,~,~,~,USELQNSNAMING)

if nargin < 5
    useLQNSnaming = false;
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
                SN(eidx) = self.servt(eidx);
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
                SN(aidx) = SN(aidx) + self.results{end,e}.RN(serverIdx,c) * self.lqn.callproc{cidx}.getMean();
                if isnan(QN(aidx)), QN(aidx)=0; end
                QN(aidx) = QN(aidx) + self.results{end,e}.QN(serverIdx,c);
            case LayeredNetworkElement.ACTIVITY
                aidx = self.ensemble{e}.classes{c}.attribute(2);
                tidx = self.lqn.parent(aidx);
                QN(tidx) = QN(tidx) + self.results{end,e}.QN(serverIdx,c);
                if isnan(TN(aidx)), TN(aidx)=0; end
                if isnan(QN(aidx)), QN(aidx)=0; end
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
                WN(aidx) = WN(aidx) + self.results{end,e}.WN(serverIdx,c);
                WN(tidx) = WN(tidx) + self.results{end,e}.WN(serverIdx,c);
                if isnan(QN(aidx)), QN(aidx)=0; end
                QN(aidx) = QN(aidx) + self.results{end,e}.QN(serverIdx,c);
        end
    end
end

for e=1:self.lqn.nentries
    eidx = self.lqn.eshift + e;
    tidx = self.lqn.parent(eidx);
    if isnan(UN(tidx)), UN(tidx)=0; end
    UN(eidx) = TN(eidx)*SN(eidx);
    for aidx=self.lqn.actsof{tidx}
        UN(aidx) = TN(aidx)*SN(aidx);
    end
    UN(tidx) = UN(tidx) + UN(eidx);
end

if ~useLQNSnaming % if LN standard naming
    QN = UN;
    UN = PN;
    RN = SN;
end
end
