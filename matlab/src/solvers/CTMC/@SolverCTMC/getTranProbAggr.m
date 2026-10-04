function [Pi_t, SSnode_a] = getTranProbAggr(self, node)
% [PI_T, SSNODE_A] = GETTRANPROBSTATEAGGR(NODE)


if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    [Pi_t, SSnode_a] = CPPLINE.tranProb(self.name, self.model, self.options, node.index, true);
    return
end

self.assertPhaseTypeStates('getTranProbAggr');

options = self.getOptions;
if isfield(options,'timespan')  && isfinite(options.timespan(2))
    sn = self.getStruct;
    % OUTPUT 11 IS THE AGGREGATE STATE SPACE. This read output 12, which is
    % EventFiltration -- a CELL of per-event matrices -- so the slice below
    % returned a 1x1 cell rather than the node's per-class counts, whatever the
    % model. See getTranProb.m for the analyzer's output list.
    [t,pi_t,~,~,~,~,~,~,~,~,SSa] = solver_ctmc_transient_analyzer(sn, options);
    jnd = node.index;
    SSnode_a = SSa(:,(jnd-1)*sn.nclasses+1:jnd*sn.nclasses);
    Pi_t = [t, pi_t];
else
    line_error(mfilename,'getTranProbAggr in SolverCTMC requires to specify a finite timespan T, e.g., SolverCTMC(model,''timespan'',[0,T]).');
end
end