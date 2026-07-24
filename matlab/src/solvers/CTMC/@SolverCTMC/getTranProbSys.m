function [Pi_t, SSsys] = getTranProbSys(self)
% [PI_T, SSSYS] = GETTRANPROBSYSSTATE()

options = self.getOptions;
if isfield(options,'timespan')  && isfinite(options.timespan(2))
    sn = self.getStruct;
    [t,pi_t,~,~,~,~,~,~,~,~,SSsys]  = solver_ctmc_transient_analyzer(sn, options);
    Pi_t = [t, pi_t];
else
    line_error(mfilename,'getTranProbSys in SolverCTMC requires to specify a finite timespan T, e.g., SolverCTMC(model,''timespan'',[0,T]).');
end
end