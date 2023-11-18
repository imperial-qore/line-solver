function solver = chooseSolverHeur(self, method)
% SOLVER = CHOOSESOLVERHEUR(METHOD)
%
switch method
    case {'getAvgChainTable', 'getAvgTputTable', 'getAvgRespTTable', ...
            'getAvgUtilTable',  'getAvgSysTable', 'getAvgNodeTable', ...
            'getAvgTable', 'getAvg', 'getAvgChain', 'getAvgSys', ...
            'getAvgNode', 'getAvgArvRChain', 'getAvgQLenChain', ...
            'getAvgUtilChain', 'getAvgRespTChain', 'getAvgTputChain', ...
            'getAvgSysRespT', 'getAvgSysTput'}
        solver = chooseAvgSolverHeur(self);
    case {'getTranAvg'}
        this_model = self.model;
        switch class(this_model)
            case 'Network'
                if SolverFluid.supports(this_model)
                    solver = self.solvers{self.CANDIDATE_FLUID};
                else
                    solver = self.solvers{self.CANDIDATE_JMT};
                end
            case 'LayeredNetwork'
                line_error(mfilename,'Method not yet supported with LayeredNetworks');
        end
    case {'getCdfRespT'}
        this_model = self.model;
        switch class(this_model)
            case 'Network'
                if this_model.hasHomogeneousScheduling(SchedStrategy.FCFS) && SolverNC.supports(this_model) && this_model.hasProductFormSolution()
                    solver = self.solvers{self.CANDIDATE_NC};
                elseif SolverFluid.supports(this_model)
                    solver = self.solvers{self.CANDIDATE_FLUID};
                else
                    solver = self.solvers{self.CANDIDATE_JMT};
                end
            case 'LayeredNetwork'
                line_error(mfilename,'Method not yet supported with LayeredNetworks');
        end
    case {'getTranCdfPassT','getTranCdfRespT'}
        this_model = self.model;
        switch class(this_model)
            case 'Network'
                if SolverFluid.supports(this_model)
                    solver = self.solvers{self.CANDIDATE_FLUID};
                else
                    solver = self.solvers{self.CANDIDATE_JMT};
                end
            case 'LayeredNetwork'
                line_error(mfilename,'Method not yet supported with LayeredNetworks');
        end
    case {'getTranProb','getTranProbSys','getTranProbAggr','getTranProbSysAggr'}
        this_model = self.model;
        switch class(this_model)
            case 'Network'
                solver = self.solvers{self.CANDIDATE_CTMC};
            case 'LayeredNetwork'
                line_error(mfilename,'Method not yet supported with LayeredNetworks');
        end
    case {'sample','sampleSys'}
        this_model = self.model;
        switch class(this_model)
            case 'Network'
                solver = self.solvers{self.CANDIDATE_SSA};
            case 'LayeredNetwork'
                line_error(mfilename,'Method not yet supported with LayeredNetworks');
        end
    case {'sampleAggr','sampleSysAggr'}
        this_model = self.model;
        switch class(this_model)
            case 'Network'
                solver = self.solvers{self.CANDIDATE_JMT};
            case 'LayeredNetwork'
                line_error(mfilename,'Method not yet supported with LayeredNetworks');
        end        
    case {'getProb','getProbAggr','getProbSys','getProbSysAggr','getProbNormConstAggr'}
        this_model = self.model;
        switch class(this_model)
            case 'Network'
                if SolverNC.supports(this_model) && this_model.hasProductFormSolution()
                    solver = self.solvers{self.CANDIDATE_NC};
                else
                    solver = self.solvers{self.CANDIDATE_JMT};
                end
            case 'LayeredNetwork'
                line_error(mfilename,'Method not yet supported with LayeredNetworks');
        end
end
end