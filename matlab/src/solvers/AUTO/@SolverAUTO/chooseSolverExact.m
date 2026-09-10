function solver = chooseSolverExact(self, method)
% SOLVER = CHOOSESOLVEREXACT(METHOD)
%
% Ranked choice restricted to solvers that return an exact answer for the
% requested metric. Exactness overrides the global MVA > NC order: NC is the
% exact normalizing-constant solution and CTMC is exact by construction. The
% call ERRORS rather than silently returning an approximation.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

switch class(self.model)
    case 'Network'
        traits = self.solverTraits();
        switch method
            case {'getTranProb','getTranProbSys','getTranProbAggr','getTranProbSysAggr', ...
                    'getCdfRespT','getCdfPassT','getPerctRespT','getTranCdfPassT','getTranCdfRespT'}
                order = self.CANDIDATE_CTMC;
            case {'getTranAvg'}
                order = self.CANDIDATE_CTMC;
            case {'sample','sampleSys','sampleAggr','sampleSysAggr'}
                % A sample path is exact in distribution, not in the mean.
                order = [self.CANDIDATE_SSA, self.CANDIDATE_LDES];
            case {'getProb','getProbAggr','getProbSys','getProbSysAggr','getProbMarg','getProbNormConstAggr'}
                if traits.isProductForm
                    order = [self.CANDIDATE_NC, self.CANDIDATE_CTMC];
                else
                    order = self.CANDIDATE_CTMC;
                end
            otherwise
                if traits.isProductForm && ~traits.hasMultiServer
                    order = [self.CANDIDATE_NC, self.CANDIDATE_CTMC];
                else
                    order = [self.CANDIDATE_CTMC, self.CANDIDATE_NC];
                end
        end
        % Gate on the method-level rule, not just the feature set: NC and MVA
        % reject 'exact' on a non-product-form model, which a flat feature set
        % cannot express.
        solver = self.chooseSolverRanked(order, 'exact');
        if isempty(solver)
            line_error(mfilename, ['No exact solver supports this model for method ''%s''. ' ...
                'Use method ''default'' for the approximate heuristic.'], method);
        end
    case 'LayeredNetwork'
        % No layered solver is exact; NC layers are the closest available.
        solver = self.chooseSolverRanked([self.CANDIDATE_LN_NC, self.CANDIDATE_LN_MVA]);
        if isempty(solver)
            line_error(mfilename,'No LayeredNetwork solver available');
        end
    case 'Environment'
        solver = self.chooseSolverRanked([self.CANDIDATE_ENV_NC, self.CANDIDATE_ENV_MVA]);
        if isempty(solver)
            line_error(mfilename,'No Environment solver available');
        end
end
end
