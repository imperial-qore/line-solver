function solver = chooseSolverRanked(self, order, methodToken)
% SOLVER = CHOOSESOLVERRANKED(ORDER, METHODTOKEN)
%
% First solver in ORDER (a vector of CANDIDATE_* slot ids) that exists and
% whose feature set accepts the model. Returns [] when none qualifies, so the
% caller can fall back rather than hand delegate() an infeasible solver.
%
% METHODTOKEN, when given, tightens the gate from supports(model) to
% supportsModelMethod(METHODTOKEN): that is the method-level gate, and it is
% where the rules a flat feature set cannot express live (product form for
% 'exact', finite capacity, NC 'mem' applicability).
%
% Feature support is necessary but not sufficient for CTMC: the chain must
% also fit memory, so the slot is screened by SolverCTMC.isStateSpaceTractable.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3
    methodToken = '';
end

solver = [];
for k = 1:numel(order)
    s = order(k);
    if s < 1 || s > numel(self.solvers) || isempty(self.solvers{s})
        continue
    end
    candidate = self.solvers{s};
    try
        if ~isempty(methodToken) && ismethod(candidate, 'supportsModelMethod')
            supported = candidate.supportsModelMethod(methodToken);
        else
            % LayeredNetwork solvers return [bool, reason]; take the first output.
            supported = candidate.supports(self.model);
        end
    catch
        supported = true;
    end
    if supported && s == self.CANDIDATE_CTMC
        supported = SolverCTMC.isStateSpaceTractable(self.model, candidate.getOptions());
    end
    if supported
        % Pin the delegate's method: the gate above only ASKED whether the
        % candidate can run METHODTOKEN, and delegate() hands the solver its
        % own options, so without this an exactness-gated choice would still
        % run the solver's default (approximate) method.
        opts = candidate.getOptions();
        if isempty(methodToken)
            opts.method = 'default';
        else
            opts.method = methodToken;
        end
        candidate.setOptions(opts);
        solver = candidate;
        return
    end
end
end
