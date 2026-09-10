function varargout = getAvgTableLayered(self,varargin)
% [AVGTABLE,QT,UT,RT,WT,AT,TT] = GETAVGTABLELAYERED()
%
% Per-LQN-element metrics table for a Java-backed solver run on a
% LayeredNetwork model (SolverLDES). Formatting mirrors
% SolverLN/SolverLQNS getAvgTable.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
% The result recorder captures the returned table together with the solver
% that produced it, so cross-codebase parity is asserted against the values a
% solver RETURNED rather than the text it printed. Off unless a run asked for
% it (LineResultRecorder.enable), and then it costs one appdata lookup here.
% The wrapper exists so that recording happens on EVERY exit path, including
% the early returns inside the implementation below.
[scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
[varargout{1:max(nargout,1)}] = getAvgTableLayered_impl(self,varargin{:});
LineResultRecorder.capture(scope, self, 'avg', varargout{1});
end

function [AvgTable,QT,UT,RT,WT,AT,TT] = getAvgTableLayered_impl(self)
% GETAVGTABLELAYERED_IMPL Implementation of GETAVGTABLELAYERED; see the wrapper above.

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
