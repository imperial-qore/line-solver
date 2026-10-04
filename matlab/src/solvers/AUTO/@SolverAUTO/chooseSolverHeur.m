function solver = chooseSolverHeur(self, method)
% SOLVER = CHOOSESOLVERHEUR(METHOD)
%
% Metric-aware ranked choice. Each accessor group maps to an ordered list of
% candidate slots; the first feasible one wins. LDES is the only simulation
% candidate and MVA > NC > MAM the analytical order, inverted to NC > MVA on
% cache models, with Fluid promoted where a smooth or high-load answer is
% wanted (sensitivity, distributions, transients, large populations).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

switch class(self.model)
    case 'Network'
        solver = chooseNetworkSolver(self, method);
    case 'LayeredNetwork'
        solver = chooseLayeredSolver(self, method);
    case 'Environment'
        % Fluid leads: the blending method restarts each stage from the mean
        % state left by the previous one, so an inner solver without transient
        % analysis returns zeros for every blended metric.
        solver = self.chooseSolverRanked([self.CANDIDATE_ENV_FLUID, ...
            self.CANDIDATE_ENV_MVA, self.CANDIDATE_ENV_NC]);
        if isempty(solver)
            line_error(mfilename,'No Environment solver available');
        end
    otherwise
        solver = self.chooseAvgSolverHeur();
end
end

function solver = chooseNetworkSolver(self, method)
% SOLVER = CHOOSENETWORKSOLVER(METHOD)

order = [];
switch method
    case {'getAvgChainTable', 'getAvgTputTable', 'getAvgRespTTable', ...
            'getAvgUtilTable',  'getAvgSysTable', 'getAvgNodeTable', ...
            'getAvgTable', 'getAvgTableLayered', 'getAvg', 'getAvgChain', 'getAvgSys', ...
            'getAvgNode', 'getAvgNodeChain', 'getAvgArvRChain', 'getAvgQLenChain', ...
            'getAvgUtilChain', 'getAvgRespTChain', 'getAvgTputChain', ...
            'getAvgSysRespT', 'getAvgSysTput', ...
            'getAvgQLen', 'getAvgUtil', 'getAvgRespT', 'getAvgResidT', ...
            'getAvgWaitT', 'getAvgTput', 'getAvgArvR', 'getAvgQLenTable', ...
            'getAvgResidTChain', 'getAvgNodeQLenChain', 'getAvgNodeUtilChain', ...
            'getAvgNodeRespTChain', 'getAvgNodeResidTChain', 'getAvgNodeTputChain', ...
            'getAvgNodeArvRChain', 'getAvgNodeChainTable', ...
            'getResults', 'hasResults', 'getAvgHandles', 'getTranHandles', ...
            'getAvgQLenHandles', 'getAvgUtilHandles', 'getAvgRespTHandles', ...
            'getAvgTputHandles', 'getAvgArvRHandles', 'getAvgResidTHandles', ...
            'getMethodFeatureSet', 'supportsModelMethod', 'isStochasticMethod', ...
            'libraries', 'showLibraryAttribution', 'citations'}
        solver = self.chooseAvgSolverHeur();
        return
    case {'getEnsembleAvg', 'getEnsembleAvgTables', 'getSolver', 'setSolver', ...
            'getNumberOfModels', 'getIteration', 'get_state', 'set_state', ...
            'update_solver'}
        line_error(mfilename,'Method %s is only available for LayeredNetwork models', method);
    case {'getTranAvg'}
        order = [self.CANDIDATE_FLUID, self.CANDIDATE_LDES];
    case {'getCdfRespT','getCdfPassT','getPerctRespT'}
        % NC gives the exact passage-time distribution on FCFS product form;
        % otherwise Fluid is the smooth approximation, then the simulators.
        if self.model.hasHomogeneousScheduling(SchedStrategy.FCFS) && self.model.hasProductFormSolution()
            order = [self.CANDIDATE_NC, self.CANDIDATE_FLUID, self.CANDIDATE_LDES];
        else
            order = [self.CANDIDATE_FLUID, self.CANDIDATE_LDES];
        end
    case {'getTranCdfPassT','getTranCdfRespT'}
        order = [self.CANDIDATE_FLUID, self.CANDIDATE_LDES];
    case {'getTranProb','getTranProbSys','getTranProbAggr','getTranProbSysAggr'}
        order = self.CANDIDATE_CTMC;
    case {'sample','sampleSys'}
        % Event-level trajectories: SSA is the native sample-path engine.
        order = [self.CANDIDATE_SSA, self.CANDIDATE_LDES];
    case {'sampleAggr','sampleSysAggr'}
        order = [self.CANDIDATE_LDES, self.CANDIDATE_SSA];
    case {'getProb','getProbAggr','getProbSys','getProbSysAggr','getProbMarg','getProbNormConstAggr'}
        if self.model.hasProductFormSolution()
            order = [self.CANDIDATE_NC, self.CANDIDATE_CTMC, self.CANDIDATE_LDES];
        else
            order = [self.CANDIDATE_CTMC, self.CANDIDATE_LDES];
        end
    case {'getAvgCacheTable','getAvgCacheT','getAvgItemTable','getAvgItemT', ...
            'cacheAvgT','itemAvgT','aCaT','aIT'}
        % Cache metrics invert the analytical order: NC leads MVA here.
        order = [self.CANDIDATE_NC, self.CANDIDATE_MVA, self.CANDIDATE_FLUID, ...
            self.CANDIDATE_CTMC, self.CANDIDATE_LDES];
    case {'getAvgLossTable','getAvgLossT','getAvgRegionLossTable','getAvgRegionLossT', ...
            'lossAvgT','regionLossAvgT','aLT','aRLT'}
        % Drop counts come from a sample path or an exact chain, not from AMVA.
        order = [self.CANDIDATE_LDES, self.CANDIDATE_CTMC, self.CANDIDATE_SSA];
    case {'getAvgOrbitTable','getAvgOrbitT','getAvgOrbit','orbitAvgT','aOT'}
        order = [self.CANDIDATE_MVA, self.CANDIDATE_CTMC, self.CANDIDATE_LDES];
    case {'getMomentTable','getMomentChainTable','getMomentStationTable', ...
            'getMomentT','getMomentChainT','getMomentStationT', ...
            'momentT','momentChainT','momentStationT','mT','mCT','mST'}
        order = [self.CANDIDATE_MVA, self.CANDIDATE_NC, self.CANDIDATE_CTMC, self.CANDIDATE_LDES];
    case {'getSensitivityTable','getSensitivityT','sensitivityT','sT', ...
            'supportsExactSensitivity'}
        % Sensitivity wants a differentiable model, which is what Fluid gives.
        order = [self.CANDIDATE_FLUID, self.CANDIDATE_MVA, self.CANDIDATE_NC];
    otherwise
        solver = self.chooseAvgSolverHeur();
        return
end

solver = self.chooseSolverRanked(order);
if isempty(solver)
    % No solver in the metric's ranking supports the model: the average
    % heuristic is the floor, and delegate() still retries the other
    % candidates if that one cannot serve the method either.
    solver = self.chooseAvgSolverHeur();
end
end

function solver = chooseLayeredSolver(self, method)
% SOLVER = CHOOSELAYEREDSOLVER(METHOD)

switch method
    case {'getTranAvg','getCdfRespT','getCdfPassT','getPerctRespT', ...
            'getTranCdfPassT','getTranCdfRespT'}
        order = [self.CANDIDATE_LN_FLUID, self.CANDIDATE_LN_MVA, self.CANDIDATE_LN_NC];
    case {'sample','sampleSys','sampleAggr','sampleSysAggr'}
        order = [self.CANDIDATE_LQNS, self.CANDIDATE_LN_MVA, self.CANDIDATE_LN_NC];
    case {'getProb','getProbAggr','getProbSys','getProbSysAggr','getProbMarg', ...
            'getProbNormConstAggr','getTranProb','getTranProbSys', ...
            'getTranProbAggr','getTranProbSysAggr'}
        order = [self.CANDIDATE_LN_NC, self.CANDIDATE_LN_MVA];
    case {'getEnsembleAvg', 'getEnsembleAvgTables', 'getSolver', 'setSolver', ...
            'getNumberOfModels', 'getIteration', 'get_state', 'set_state', ...
            'update_solver'}
        order = [self.CANDIDATE_LN_MVA, self.CANDIDATE_LN_NC, self.CANDIDATE_LN_MAM, ...
            self.CANDIDATE_LN_FLUID, self.CANDIDATE_LQNS];
    otherwise
        solver = self.chooseAvgSolverHeur();
        return
end

solver = self.chooseSolverRanked(order);
if isempty(solver)
    solver = self.chooseAvgSolverHeur();
end
end
