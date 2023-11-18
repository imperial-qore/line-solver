function [Pi_t, SSnode_a] = getTranProbAggr(self, node)
% [PI_T, SSNODE_A] = GETTRANPROBSTATEAGGR(NODE)
if GlobalConstants.DummyMode
    Pi_t = NaN;
    SSnode_a = NaN;
    return
end
line_error(mfilename,'Method not implemented yet.')
options = self.getOptions;
initSeed = self.options.seed;
if isfield(options,'timespan')  && isfinite(options.timespan(2))
    tu = [];
    TranSysStateAggr = cell(options.iter_max,1);
    for it=1:options.iter_max
        self.options.seed = initSeed + it -1;
        TranSysStateAggr{it} = self.sampleSysAggr;
        tu = union(tu, TranSysStateAggr{it}.t);
    end
    for j=1:length(TranSysStateAggr{it}.state)
        for it=1:options.iter_max
            TranSysStateAggr{it}.state{j} = interp1(TranSysStateAggr{it}.t, TranSysStateAggr{it}.state{j}, tu, 'previous');
            
        end
    end
    
    %                self.options.seed = initSeed; % in case of interruption
    %                jnd = node.index;
    %                SSnode_a = SSa(:,(jnd-1)*sn.nclasses+1:jnd*sn.nclasses);
    %                Pi_t = [t, pi_t];
    SSnode_a = [];
    Pi_t = [];
else
    line_error(mfilename,'getTranProbAggr in SolverJMT requires to specify a finite timespan T, e.g., SolverJMT(model,''timespan'',[0,T]).');
end
end