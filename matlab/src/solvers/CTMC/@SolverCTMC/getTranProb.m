function [Pi_t, SSnode] = getTranProb(self, node)
% [PI_T, SSNODE] = GETTRANPROBSTATE(NODE)


if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    [Pi_t, SSnode] = CPPLINE.tranProb(self.name, self.model, self.options, node.index, false);
    return
end

self.assertPhaseTypeStates('getTranProb');

options = self.getOptions;
if isfield(options,'timespan')  && isfinite(options.timespan(2))
    sn = self.getStruct;
    % OUTPUT 10 IS THE DETAILED STATE SPACE. This read output 11, which is
    % StateSpaceAggr, and then sliced it with the DETAILED block widths, so the
    % columns returned as "this node's state" were per-class job counts cut at
    % encoding offsets. The analyzer's list is
    % [t,pit,QNt,UNt,RNt,TNt,CNt,XNt,InfGen,StateSpace,StateSpaceAggr,EventFiltration,...].
    [t,pi_t,~,~,~,~,~,~,~,SS] = solver_ctmc_transient_analyzer(sn, options);
    % THE BLOCK WIDTHS ARE THE LOCAL SPACE'S, not length(sn.state{isf}): every
    % row of a node's block is stored at that node's WIDEST encoding, so a
    % declared state narrower than it shifts every later node's slice left.
    % See _kb/11-conventions-and-gotchas.md (CTMC state-vector padding).
    [~, localStateSpace] = getStateSpace(self);
    jnd = node.index;
    shift = 1;
    for isf = 1:sn.nstateful
        len = size(localStateSpace{isf}, 2);
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