function [logNormConst] = getProbNormConstAggr(self)
% [LOGNORMCONST] = GETPROBNORMCONST()

if GlobalConstants.DummyMode
    logNormConst = NaN;
    return
end

runAnalyzer(self);
logNormConst = self.result.Prob.logNormConstAggr;
end