function [QN,UN,RN,TN,AN,WN] = getEnsembleAvg(self)
% [QN,UN,RN,TN,AN,WN] = GETENSEMBLEAVG(SELF)

runAnalyzer(self);
%QN = self.result.Avg.QLen;
UN = self.result.Avg.Util;
RN = self.result.Avg.RespT;
TN = self.result.Avg.Tput;
PN = self.result.Avg.ProcUtil;
SN = self.result.Avg.SvcT;
AN = TN*0;
WN = RN*0;
QN = UN;
UN = PN;
RN = SN;
end
