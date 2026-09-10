function solver = chooseSolverSim(self, method)
% SOLVER = CHOOSESOLVERSIM(METHOD)
%
% Ranked choice restricted to simulators. LDES leads everywhere except
% event-level sampling, where SSA is the native sample-path engine.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

switch class(self.model)
    case 'Network'
        switch method
            case {'sample','sampleSys'}
                order = [self.CANDIDATE_SSA, self.CANDIDATE_LDES];
            otherwise
                order = [self.CANDIDATE_LDES, self.CANDIDATE_SSA];
        end
        solver = self.chooseSolverRanked(order);
        if isempty(solver)
            line_error(mfilename,'No simulator supports this model');
        end
    case 'LayeredNetwork'
        % lqsim is the layered simulator; the LN solvers are the fallback.
        solver = self.chooseSolverRanked([self.CANDIDATE_LQNS, self.CANDIDATE_LN_MVA, ...
            self.CANDIDATE_LN_NC]);
        if isempty(solver)
            line_error(mfilename,'No LayeredNetwork simulator available');
        end
    case 'Environment'
        solver = self.chooseSolverRanked([self.CANDIDATE_ENV_MVA, self.CANDIDATE_ENV_NC, ...
            self.CANDIDATE_ENV_FLUID]);
        if isempty(solver)
            line_error(mfilename,'No Environment solver available');
        end
end
end
