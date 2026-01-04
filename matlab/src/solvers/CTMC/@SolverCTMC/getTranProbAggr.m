function [Pi_t, SSnode_a] = getTranProbAggr(self, node)
% [PI_T, SSNODE_A] = GETTRANPROBSTATEAGGR(NODE)

options = self.getOptions;
if isfield(options,'timespan')  && isfinite(options.timespan(2))
    sn = self.getStruct;
    [t,pi_t,~,~,~,~,~,~,~,~,~,SSa] = solver_ctmc_transient_analyzer(sn, options);
    jnd = node.index;
    SSnode_a = SSa(:,(jnd-1)*sn.nclasses+1:jnd*sn.nclasses);
    Pi_t = [t, pi_t];
else
    line_error(mfilename,'getTranProbAggr in SolverCTMC requires to specify a finite timespan T, e.g., SolverCTMC(model,''timespan'',[0,T]).');
end
end