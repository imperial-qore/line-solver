function [Pi_t, SSnode] = getTranProb(self, node)
% [PI_T, SSNODE] = GETTRANPROBSTATE(NODE)

options = self.getOptions;
if isfield(options,'timespan')  && isfinite(options.timespan(2))
    sn = self.getStruct;
    [t,pi_t,~,~,~,~,~,~,~,~,SS] = solver_ctmc_transient_analyzer(sn, options);
    jnd = node.index;
    shift = 1;
    for isf = 1:sn.nstateful
        len = length(sn.state{isf});
        if sn.statefulToNode(isf) == jnd
            SSnode = SS(:,shift:shift+len-1);
            break;
        end
        shift = shift+len;
    end
    Pi_t = [t, pi_t];
else
    line_error(mfilename,'getTranProb in SolverCTMC requires to specify a finite timespan T, e.g., SolverCTMC(model,''timespan'',[0,T]).');
end
end