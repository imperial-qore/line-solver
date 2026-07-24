function [NodeAvgTable,CallAvgTable] = getRawAvgTables(self)
% [QN,UN,RN,TN] = GETRAWAVGTABLES(SELF)

runAnalyzer(self);

lqn = self.getStruct;
Node = label(lqn.names);
O = length(Node);
NodeType = label(O,1);
for o = 1:O
    switch lqn.type(o)
        case LayeredNetworkElement.PROCESSOR
            NodeType(o,1) = label({'Processor'});
        case LayeredNetworkElement.TASK
            NodeType(o,1) = label({'Task'});
        case LayeredNetworkElement.ENTRY
            NodeType(o,1) = label({'Entry'});
        case LayeredNetworkElement.ACTIVITY
            NodeType(o,1) = label({'Activity'});
        case LayeredNetworkElement.CALL
            NodeType(o,1) = label({'Call'});
    end
end
Utilization = self.result.RawAvg.Nodes.Utilization;
Phase1Utilization = self.result.RawAvg.Nodes.Phase1Utilization;
Phase2Utilization = self.result.RawAvg.Nodes.Phase2Utilization;
Phase1ServiceTime = self.result.RawAvg.Nodes.Phase1ServiceTime;
Phase2ServiceTime = self.result.RawAvg.Nodes.Phase2ServiceTime;
Throughput = self.result.RawAvg.Nodes.Throughput;
ProcWaiting = self.result.RawAvg.Nodes.ProcWaiting;
ProcUtilization = self.result.RawAvg.Nodes.ProcUtilization;
NodeAvgTable = Table(Node, NodeType, Utilization, Phase1Utilization,...
    Phase2Utilization, Phase1ServiceTime, Phase2ServiceTime, Throughput,...
    ProcWaiting, ProcUtilization);

lqn=self.getStruct;
if lqn.ncalls == 0
    CallAvgTable = Table();
else
    SourceNode = label({lqn.names{lqn.callpair(:,1)}})';
    TargetNode = label({lqn.names{lqn.callpair(:,2)}})';
    Type = lqn.calltype;
    Waiting = self.result.RawAvg.Edges.Waiting(1:lqn.ncalls);
    CallAvgTable = Table(SourceNode, TargetNode, Type, Waiting);
end
end
