function ProbSysAggr = getProbSysAggr(self)
% PROBSYSAGGR = GETPROBSYSAGGR()
% Aggregated joint system-state probability. For LDES trajectories are already
% per-class, so this equals getProbSys(). Fully JSON-mediated. Mirrors the
% Python-native getProbSysAggr().
ProbSysAggr = self.getProbSys();
end
