function [Pi_t, SSsys] = getTranProbSys(self)
% [PI_T, SSSYS] = GETTRANPROBSYSSTATE()


if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    % A CHAIN HAS NO MODEL TO SERIALIZE. isChainSolver means the caller handed
    % SolverCTMC a bare generator rather than a Network, and linemodel_save has
    % nothing to write for it, so this one keeps a refusal while the network
    % case routes.
    if self.isChainSolver()
        CPPLINE.cppUnsupported(self.name, 'getTranProbSys', ...
            ['a chain given directly to SolverCTMC is not a Network, so there is no ' ...
            'model.json for line-cli to read']);
    end
    [Pi_t, SSsys] = CPPLINE.tranProb(self.name, self.model, self.options, [], false);
    return
end

if self.isChainSolver()
    % Chain mode: integrate (CTMC) or iterate (DTMC) from options.init_sol,
    % defaulting to the uniform distribution.
    options = self.getOptions;
    if ~isfield(options,'timespan') || ~isfinite(options.timespan(2))
        line_error(mfilename,'getTranProbSys requires a finite timespan T, e.g., SolverCTMC(chain,''timespan'',[0,T]).');
    end
    [pi_t, t] = solver_ctmc_chain_transient(self.chainModel, options.init_sol, options.timespan);
    Pi_t = [t, pi_t];
    SSsys = self.getStateSpace();
    return
end

self.assertPhaseTypeStates('getTranProbSys');

options = self.getOptions;
if isfield(options,'timespan')  && isfinite(options.timespan(2))
    sn = self.getStruct;
    % Output 10 is StateSpace; 11 is StateSpaceAggr, which is what this used to
    % return under the name SSsys. See getTranProb.m for the analyzer's list.
    [t,pi_t,~,~,~,~,~,~,~,SSsys]  = solver_ctmc_transient_analyzer(sn, options);
    Pi_t = [t, pi_t];
else
    line_error(mfilename,'getTranProbSys in SolverCTMC requires to specify a finite timespan T, e.g., SolverCTMC(model,''timespan'',[0,T]).');
end
end