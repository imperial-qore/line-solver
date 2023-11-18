function [QN,UN,RN,TN] = getStageAvg(self,~,~,~,~)
%             % [QN,UN,RN,TN] = GETSTAGEAVG(SELF,~,~,~,~)
%
%             runAnalyzer(self); % run iterations
%             E = self.nlayers;
%             QN  = zeros(E,self.lqn.nidx+self.lqn.ncalls);
%             UN  = zeros(E,self.lqn.nidx+self.lqn.ncalls);
%             RN  = zeros(E,self.lqn.nidx+self.lqn.ncalls);
%             TN  = zeros(E,self.lqn.nidx+self.lqn.ncalls);
%             for e=1:E
%                 clientIdx = self.ensemble{e}.attribute.clientIdx;
%                 serverIdx = self.ensemble{e}.attribute.serverIdx;
%                 sourceIdx = self.ensemble{e}.attribute.sourceIdx;
%                 for c=1:self.ensemble{e}.getNumberOfClasses
%                     type = self.ensemble{e}.classes{c}.attribute(1);
%                     switch type
%                         case LayeredNetworkElement.TASK
%                             tidx = self.ensemble{e}.classes{c}.attribute(2);
%                             if strcmp(self.lqn.sched(tidx), SchedStrategy.REF)
%                                 TN(e,tidx) = self.results{end,e}.TN(clientIdx,c);
%                             end
%                             UN(e,tidx) = NaN;
%                             RN(e,tidx) = NaN;
%                             QN(e,tidx) = NaN;
%                             %if ~isnan(clientIdx)
%                             %    TN(e,tidx) = max(self.results{end,e}.TN(clientIdx,c), TN(e,tidx));
%                             %end
%                             %UN(e,tidx) = self.util(tidx);
%                             %RN(e,tidx) = NaN;
%                             %QN(e,tidx) = self.results{end,e}.QN(serverIdx,c);
%                         case LayeredNetworkElement.ENTRY
%                             %eidx = self.ensemble{e}.classes{c}.attribute(2);
%                             %RN(e,eidx) = NaN;
%                             %UN(e,eidx) = self.results{end,e}.UN(serverIdx,c);
%                             %QN(e,eidx) = self.results{end,e}.QN(serverIdx,c);
%                             %TN(e,eidx) = self.results{end,e}.TN(clientIdx,c);
%                             %if ~isnan(clientIdx)
%                             %    TN(e,eidx) = max(self.results{end,e}.TN(clientIdx,c), TN(e,eidx));
%                             %end
%                         case LayeredNetworkElement.CALL
%                             cidx = self.ensemble{e}.classes{c}.attribute(2);
%                             idx = self.lqn.nidx + cidx;
%                             aidx = self.lqn.callpair(cidx,1);
%                             tidx = self.lqn.parent(aidx);
%                             % contribution of call to task
%                             %UN(e,tidx) = UN(e,tidx) + self.results{end,e}.UN(serverIdx,c);
%                             %QN(e,tidx) = QN(e,tidx) + self.results{end,e}.QN(serverIdx,c);
%                             switch self.ensemble{e}.classes{c}.type
%                                 case JobClassType.CLOSED
%                                     if self.ensemble{e}.classes{c}.completes
%                                         %     TN(e,tidx) = TN(e,tidx) + max(self.results{end,e}.TN(serverIdx,c), self.results{end,e}.TN(clientIdx,c));
%                                     end
%                                     %TN(e,aidx) = TN(e,aidx) + self.results{end,e}.TN(serverIdx,c);
%                                     TN(e,idx) = TN(e,idx) + self.results{end,e}.TN(serverIdx,c);
%                                     %TN(e,idx) = max(self.results{end,e}.TN(serverIdx,c), self.results{end,e}.TN(clientIdx,c));
%                                 case JobClassType.OPEN
%                                     if self.ensemble{e}.classes{c}.completes
%                                         %    TN(e,tidx) = TN(e,tidx) + self.results{end,e}.TN(sourceIdx,c);
%                                     end
%                                     %TN(e,aidx) = TN(e,aidx) + self.results{end,e}.TN(sourceIdx,c);
%                                     TN(e,idx) = TN(e,idx) + self.results{end,e}.TN(sourceIdx,c);
%                                     %TN(e,idx) = max(self.results{end,e}.TN(serverIdx,c), self.results{end,e}.TN(sourceIdx,c));
%                             end
%                             %RN(e,aidx) = RN(e,aidx) + self.results{end,e}.RN(serverIdx,c);
%                             %UN(e,aidx) = UN(e,aidx) + self.results{end,e}.UN(serverIdx,c);
%                             %QN(e,aidx) = QN(e,aidx) + self.results{end,e}.QN(serverIdx,c);
%                             RN(e,idx) = RN(e,idx) + self.results{end,e}.RN(serverIdx,c);
%                             UN(e,idx) = UN(e,idx) + self.results{end,e}.UN(serverIdx,c);
%                             QN(e,idx) = QN(e,idx) + self.results{end,e}.QN(serverIdx,c);
%                         case LayeredNetworkElement.ACTIVITY
%                             aidx = self.ensemble{e}.classes{c}.attribute(2);
%                             tidx = self.lqn.parent(aidx);
%                             % contribution of activity to task
%                             %UN(e,tidx) = UN(e,tidx) + self.results{end,e}.UN(serverIdx,c);
%                             %QN(e,tidx) = QN(e,tidx) + self.results{end,e}.QN(serverIdx,c);
%                             switch self.ensemble{e}.classes{c}.type
%                                 case JobClassType.CLOSED
%                                     if self.ensemble{e}.classes{c}.completes
%                                         %    TN(e,tidx) = TN(e,tidx) + max(self.results{end,e}.TN(serverIdx,c), self.results{end,e}.TN(clientIdx,c));
%                                     end
%                                     TN(e,aidx) = TN(e,aidx) + self.results{end,e}.TN(serverIdx,c);
%                                 case JobClassType.OPEN
%                                     if self.ensemble{e}.classes{c}.completes
%                                         %    TN(e,tidx) = TN(e,tidx) + self.results{end,e}.TN(sourceIdx,c);
%                                     end
%                                     %TN(e,aidx) = TN(e,aidx) + self.results{end,e}.TN(sourceIdx,c);
%                                     TN(e,aidx) = TN(e,aidx) + self.results{end,e}.TN(sourceIdx,c);
%                             end
%                             RN(e,aidx) = RN(e,aidx) + self.results{end,e}.RN(serverIdx,c);
%                             UN(e,aidx) = UN(e,aidx) + self.results{end,e}.UN(serverIdx,c);
%                             QN(e,aidx) = QN(e,aidx) + self.results{end,e}.QN(serverIdx,c);
%                     end
%                 end
%             end
end
