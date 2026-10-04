function solver = chooseAvgSolverHeur(self)
% SOLVER = CHOOSEAVGSOLVERHEUR()
%
% Ranked choice for mean-value metrics. The global order is MVA > NC > MAM,
% inverted to NC > MVA on cache models, with Fluid promoted for large
% populations and MAM promoted when the traffic is autocorrelated.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

switch class(self.model)
    case 'Network'
        traits = self.solverTraits();
        % Small populations: an approximation buys nothing there, so take an
        % exact solver whenever one is available, preferring MVA over NC over
        % CTMC. The 'exact' method name is what makes this a claim rather than a
        % preference: MVA and NC reject it without a product-form solution and
        % CTMC rejects it when the chain does not fit memory.
        solver = [];
        if traits.totalJobs > 0 && traits.totalJobs <= self.EXACT_POPULATION_MAX
            if traits.hasCache
                % Cache models keep the inverted analytical order.
                exactOrder = [self.CANDIDATE_NC, self.CANDIDATE_MVA, self.CANDIDATE_CTMC];
            else
                exactOrder = [self.CANDIDATE_MVA, self.CANDIDATE_NC, self.CANDIDATE_CTMC];
            end
            solver = self.chooseSolverRanked(exactOrder, 'exact');
        end
        if isempty(solver)
            order = avgOrder(self, traits);
            solver = self.chooseSolverRanked(order);
        end
        if isempty(solver)
            % Nothing in the ranked list is feasible: fall back to the whole
            % pool in the global order rather than returning nothing.
            solver = self.chooseSolverRanked([self.CANDIDATE_MVA, self.CANDIDATE_NC, ...
                self.CANDIDATE_MAM, self.CANDIDATE_FLUID, self.CANDIDATE_LDES, ...
                self.CANDIDATE_CTMC, self.CANDIDATE_SSA]);
        end
        if isempty(solver)
            line_error(mfilename,'No solver supports this model');
        end
    case 'LayeredNetwork'
        % A layered cache model needs the NC layer solver: the cache layer is
        % where NC beats MVA.
        for t = 1:length(self.model.tasks)
            if isa(self.model.tasks{t},'CacheTask')
                solver = self.solvers{self.CANDIDATE_LN_NC};
                return
            end
        end
        if ~isempty(self.solvers{self.CANDIDATE_LQNS}) && SolverLQNS.isAvailable()
            solver = self.solvers{self.CANDIDATE_LQNS};
        else
            solver = self.chooseSolverRanked([self.CANDIDATE_LN_MVA, self.CANDIDATE_LN_NC, ...
                self.CANDIDATE_LN_MAM, self.CANDIDATE_LN_FLUID]);
        end
        if isempty(solver)
            line_error(mfilename,'No LayeredNetwork solver available');
        end
    case 'Environment'
        % The blending method is transient-based: it restarts each stage from
        % the average state left by the previous one. An inner solver without
        % a transient analysis contributes nothing to that recursion and the
        % blended metrics come back as zeros, so the fluid inner solver leads
        % the ranking and the mean-value ones are kept only as a fallback.
        solver = self.chooseSolverRanked([self.CANDIDATE_ENV_FLUID, ...
            self.CANDIDATE_ENV_MVA, self.CANDIDATE_ENV_NC]);
        if isempty(solver)
            line_error(mfilename,'No Environment solver available');
        end
end
end

function order = avgOrder(self, traits)
% ORDER = AVGORDER(TRAITS)

if traits.hasCache
    order = [self.CANDIDATE_NC, self.CANDIDATE_MVA, self.CANDIDATE_FLUID, ...
        self.CANDIDATE_CTMC, self.CANDIDATE_LDES];
elseif traits.hasFCR
    if traits.totalJobs <= 10
        order = [self.CANDIDATE_NC, self.CANDIDATE_CTMC, self.CANDIDATE_LDES];
    else
        order = [self.CANDIDATE_NC, self.CANDIDATE_LDES, self.CANDIDATE_CTMC];
    end
elseif strcmp(traits.prio,'preempt')
    order = [self.CANDIDATE_MAM, self.CANDIDATE_LDES, self.CANDIDATE_CTMC, self.CANDIDATE_SSA];
elseif strcmp(traits.prio,'ps')
    order = [self.CANDIDATE_CTMC, self.CANDIDATE_LDES, self.CANDIDATE_SSA];
elseif traits.hasMAP
    order = [self.CANDIDATE_MAM, self.CANDIDATE_MVA, self.CANDIDATE_FLUID, ...
        self.CANDIDATE_LDES];
elseif strcmp(traits.prio,'hol')
    order = [self.CANDIDATE_MVA, self.CANDIDATE_MAM, self.CANDIDATE_FLUID, ...
        self.CANDIDATE_CTMC, self.CANDIDATE_LDES];
elseif traits.popPerChain > 30
    order = [self.CANDIDATE_FLUID, self.CANDIDATE_MVA, self.CANDIDATE_NC];
elseif traits.totalJobs > 0 && traits.totalJobs <= self.EXACT_POPULATION_MAX
    % No exact solver was available at this population (chooseAvgSolverHeur
    % tried first), so keep the exact-leaning approximate order.
    order = [self.CANDIDATE_NC, self.CANDIDATE_MVA, self.CANDIDATE_MAM];
elseif self.model.hasHomogeneousScheduling(SchedStrategy.INF)
    order = [self.CANDIDATE_MVA, self.CANDIDATE_NC, self.CANDIDATE_FLUID];
else
    order = [self.CANDIDATE_MVA, self.CANDIDATE_NC, self.CANDIDATE_MAM, ...
        self.CANDIDATE_FLUID, self.CANDIDATE_LDES];
end
end
