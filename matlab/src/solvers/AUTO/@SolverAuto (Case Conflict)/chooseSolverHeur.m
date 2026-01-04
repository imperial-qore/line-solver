function solver = chooseSolverHeur(self, method)
% SOLVER = CHOOSESOLVERHEUR(METHOD)
%
switch method
    case {'getAvgChainTable', 'getAvgTputTable', 'getAvgRespTTable', ...
            'getAvgUtilTable',  'getAvgSysTable', 'getAvgNodeTable', ...
            'getAvgTable', 'getAvg', 'getAvgChain', 'getAvgSys', ...
            'getAvgNode', 'getAvgArvRChain', 'getAvgQLenChain', ...
            'getAvgUtilChain', 'getAvgRespTChain', 'getAvgTputChain', ...
            'getAvgSysRespT', 'getAvgSysTput', ...
            'getAvgQLen', 'getAvgUtil', 'getAvgRespT', 'getAvgResidT', ...
            'getAvgWaitT', 'getAvgTput', 'getAvgArvR', 'getAvgQLenTable', ...
            'getAvgResidTChain', 'getAvgNodeQLenChain', 'getAvgNodeUtilChain', ...
            'getAvgNodeRespTChain', 'getAvgNodeResidTChain', 'getAvgNodeTputChain', ...
            'getAvgNodeArvRChain', 'getAvgNodeChainTable', 'getPerctRespT', ...
            'getResults', 'hasResults', 'getAvgHandles', 'getTranHandles', ...
            'getAvgQLenHandles', 'getAvgUtilHandles', 'getAvgRespTHandles', ...
            'getAvgTputHandles', 'getAvgArvRHandles', 'getAvgResidTHandles'}
        solver = chooseAvgSolverHeur(self);
    case {'getEnsembleAvg', 'getEnsembleAvgTables', 'getSolver', 'setSolver', ...
            'getNumberOfModels', 'getIteration', 'get_state', 'set_state', ...
            'update_solver'}
        % LayeredNetwork/EnsembleSolver methods - use first available LayeredNetwork solver
        this_model = self.model;
        switch class(this_model)
            case 'LayeredNetwork'
                % Try SolverLN with MVA first (fastest), then others
                if ~isempty(self.solvers{self.CANDIDATE_LN_MVA})
                    solver = self.solvers{self.CANDIDATE_LN_MVA};
                elseif ~isempty(self.solvers{self.CANDIDATE_LN_NC})
                    solver = self.solvers{self.CANDIDATE_LN_NC};
                elseif ~isempty(self.solvers{self.CANDIDATE_LN_MAM})
                    solver = self.solvers{self.CANDIDATE_LN_MAM};
                elseif ~isempty(self.solvers{self.CANDIDATE_LN_FLUID})
                    solver = self.solvers{self.CANDIDATE_LN_FLUID};
                elseif ~isempty(self.solvers{self.CANDIDATE_LQNS})
                    solver = self.solvers{self.CANDIDATE_LQNS};
                else
                    line_error(mfilename,'No LayeredNetwork solver available');
                end
            case 'Network'
                line_error(mfilename,'Method %s is only available for LayeredNetwork models', method);
        end
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