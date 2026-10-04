function [QN,UN,RN,TN,AN,WN] = getEnsembleAvg(self)
% [QN,UN,RN,TN,AN,WN] = GETENSEMBLEAVG(SELF)

runAnalyzer(self);
%QN = self.result.Avg.QLen;
UN = self.result.Avg.Util;
RN = self.result.Avg.RespT;
TN = self.result.Avg.Tput;
PN = self.result.Avg.ProcUtil;
SN = self.result.Avg.SvcT;
AN = TN*NaN;
WN = RN*NaN;
QN = UN;
UN = PN;
RN = SN;

% UN is lqns' proc-utilization, verbatim for hosts, tasks and activities and
% aggregated over the activity graph for entries, which lqns itself reports as
% 0 in the activity-graph form. Both lqns and SolverLN report the processor
% utilization summed over the host's servers, so no rescaling by the host
% multiplicity applies; see _kb/06-solver-catalog.md (LQNS wrapper).
end
