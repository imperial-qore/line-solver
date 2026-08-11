function [AvgTable,QT,UT,RT,WT,AT,TT] = getAvgTableLayered(self)
% [AVGTABLE,QT,UT,RT,WT,AT,TT] = GETAVGTABLELAYERED()
%
% Per-LQN-element metrics table for a Java-backed solver run on a
% LayeredNetwork model (SolverLDES). Formatting mirrors
% SolverLN/SolverLQNS getAvgTable.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

self.obj.getAvg(); % runs the LN LDES analyzer
avgTable = self.obj.getLNAvgTable();
[QN,UN,RN,WN,AN,TN] = JLINE.arrayListToResults(avgTable);

lqn = self.model.getStruct;
Node = label(lqn.names);
O = length(Node);
NodeType = label(O,1);
for o = 1:O
    switch lqn.type(o)
        case LayeredNetworkElement.PROCESSOR
            NodeType(o,1) = label({'Processor'});
        case LayeredNetworkElement.TASK
            if lqn.isref(o)
                NodeType(o,1) = label({'RefTask'});
            else
                NodeType(o,1) = label({'Task'});
            end
        case LayeredNetworkElement.ENTRY
            NodeType(o,1) = label({'Entry'});
        case LayeredNetworkElement.ACTIVITY
            NodeType(o,1) = label({'Activity'});
        case LayeredNetworkElement.CALL
            NodeType(o,1) = label({'Call'});
    end
end
QLen = QN;
QT = Table(Node,QLen);
Util = UN;
UT = Table(Node,Util);
RespT = RN;
RT = Table(Node,RespT);
Tput = TN;
TT = Table(Node,Tput);
ResidT = WN;
WT = Table(Node,ResidT);
ArvR = AN;
AT = Table(Node,ArvR);
AvgTable = Table(Node, NodeType, QLen, Util, RespT, ResidT, ArvR, Tput);
end
