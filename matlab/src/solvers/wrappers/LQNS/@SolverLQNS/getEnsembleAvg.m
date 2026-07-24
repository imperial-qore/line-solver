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

% LQNS reports proc-utilization summed over all instances of the host
% processor; LN reports the per-server fraction. Rescale UN to match LN.
lqn = self.getStruct;
for idx = 1:lqn.nidx
    cur = idx;
    hostMult = 1;
    for hops = 0:lqn.nidx
        if cur < 1 || cur > lqn.nidx
            break
        end
        if lqn.type(cur) == LayeredNetworkElement.PROCESSOR
            m = lqn.mult(cur);
            if m > 0 && ~isinf(m)
                hostMult = m;
            end
            break
        end
        p = lqn.parent(cur);
        if p <= 0 || p == cur
            break
        end
        cur = p;
    end
    if hostMult > 1 && ~isnan(UN(idx))
        UN(idx) = UN(idx) / hostMult;
    end
end
end
