function [QN,UN,RN,TN,AN,WN] = getAvgLayered(self)
% [QN,UN,RN,TN,AN,WN] = GETAVGLAYERED()
%
% Per-LQN-element average metrics of a LayeredNetwork model, the branch getAvg
% takes for one (getAvg itself is sealed). This default serves a Java-backed
% LQN simulation (SolverLDES); SolverLQNS overrides it with its own ensemble.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

self.obj.getAvg(); % runs the LN LDES analyzer
avgTable = self.obj.getLNAvgTable();
[QN,UN,RN,WN,AN,TN] = JLINE.arrayListToResults(avgTable);
end
