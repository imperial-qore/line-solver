function [logNormConst] = getProbNormConstAggr(self)
% [LOGNORMCONST] = GETPROBNORMCONST()


if GlobalConstants.DummyMode
    logNormConst = NaN;
    return
end

if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    logNormConst = CPPLINE.normConstAggr(self.name, self.model, self.options);
    return
end

runAnalyzer(self);
logNormConst = self.result.Prob.logNormConstAggr;
end