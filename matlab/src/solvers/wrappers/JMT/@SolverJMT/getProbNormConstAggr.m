function lNormConst = getProbNormConstAggr(self)
% LNORMCONST = GETPROBNORMCONST()

if GlobalConstants.DummyMode
    lNormConst = NaN;
    return
end

switch self.options.method
    case {'jmva','jmva.recal','jmva.comom'}%,'jmva.ls'}
        runAnalyzer(self);
        lNormConst = self.result.Prob.logNormConstAggr;
    otherwise
        lNormConst = NaN; %#ok<NASGU>
        % 'jmva.ls' was named here too until it was withdrawn from
        % listValidMethods, so the advice sent a caller to a method the name
        % gate refuses; SolverNC 'ls' is the reachable logistic sampler.
        line_error(mfilename,'Selected solver method does not compute normalizing constants. Choose either jmva, jmva.recal or jmva.comom.');
end
end