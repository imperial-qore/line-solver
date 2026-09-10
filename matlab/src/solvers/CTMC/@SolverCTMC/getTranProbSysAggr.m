function [Pi_t, SSsysa] = getTranProbSysAggr(self)
% [PI_T, SSSYSA] = GETTRANPROBSYSSTATEAGGR()


if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    [Pi_t, SSsysa] = CPPLINE.tranProb(self.name, self.model, self.options, [], true);
    return
end

self.assertPhaseTypeStates('getTranProbSysAggr');

options = self.getOptions;
if isfield(options,'timespan')  && isfinite(options.timespan(2))
    sn = self.getStruct;
    % Output 11 is StateSpaceAggr; 12 is EventFiltration, which is what this
    % used to return under the name SSsysa. See getTranProb.m for the list.
    [t,pi_t,~,~,~,~,~,~,~,~,SSsysa] = solver_ctmc_transient_analyzer(sn, options);
    Pi_t = [t, pi_t];
else
    line_error(mfilename,'getTranProbSysAggr in SolverCTMC requires to specify a finite timespan T, e.g., SolverCTMC(model,''timespan'',[0,T]).');
end
end