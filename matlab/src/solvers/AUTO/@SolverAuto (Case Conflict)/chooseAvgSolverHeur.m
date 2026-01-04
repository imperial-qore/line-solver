function solver = chooseAvgSolverHeur(self)
% SOLVER = CHOOSEAVGSOLVER()

this_model = self.model;
switch class(this_model)
    case 'Network'
        if this_model.hasSingleChain
            solver = self.solvers{self.CANDIDATE_NC}; % exact O(1) solution
        else % MultiChain
            if this_model.hasHomogeneousScheduling(SchedStrategy.INF)
                solver = self.solvers{self.CANDIDATE_MVA}; % all infinite servers
            elseif this_model.hasProductFormSolution() && ~this_model.hasMultiServer
                solver = self.solvers{self.CANDIDATE_NC};
            elseif this_model.hasHomogeneousScheduling(SchedStrategy.FCFS) && ~this_model.hasMultiServer
                if sum(this_model.getNumberOfJobs) / sum(this_model.getNumberOfChains) > 30 % likely fluid regime  
                    solver = self.solvers{self.CANDIDATE_FLUID};
                elseif sum(this_model.getNumberOfJobs) / sum(this_model.getNumberOfChains) > 10 % mid/heavy load
                    solver = self.solvers{self.CANDIDATE_MVA};
                elseif sum(this_model.getNumberOfJobs) < 5 % light load, avoid errors of AMVA in low populations
                    solver = self.solvers{self.CANDIDATE_NC};
                end
            elseif (this_model.hasHomogeneousScheduling(SchedStrategy.PS)||this_model.hasHomogeneousScheduling(SchedStrategy.PSPRIO)) && this_model.hasMultiServer
                solver = self.solvers{self.CANDIDATE_MVA};
            elseif this_model.hasHomogeneousScheduling(SchedStrategy.FCFS) && this_model.hasMultiServer
                if sum(this_model.getNumberOfJobs) < 5 % light load, avoid errors of AMVA in low populations
                    solver = self.solvers{self.CANDIDATE_NC};
                else
                    solver = self.solvers{self.CANDIDATE_MVA};
                end
            else
                solver = self.solvers{self.CANDIDATE_MVA};
            end
        end
    case 'LayeredNetwork'
        % if the LQN has caches use LN
        for t=1:length(self.model.tasks)
            if isa(self.model.tasks{t},'CacheTask')
                solver = self.solvers{self.CANDIDATE_LN_NC};
                return
            end
        end
        % Try LQNS first if available, otherwise fall back to SolverLN variants
        if ~isempty(self.solvers{self.CANDIDATE_LQNS}) && SolverLQNS.isAvailable()
            solver = self.solvers{self.CANDIDATE_LQNS};
        elseif ~isempty(self.solvers{self.CANDIDATE_LN_MVA})
            solver = self.solvers{self.CANDIDATE_LN_MVA};
        elseif ~isempty(self.solvers{self.CANDIDATE_LN_NC})
            solver = self.solvers{self.CANDIDATE_LN_NC};
        elseif ~isempty(self.solvers{self.CANDIDATE_LN_MAM})
            solver = self.solvers{self.CANDIDATE_LN_MAM};
        elseif ~isempty(self.solvers{self.CANDIDATE_LN_FLUID})
            solver = self.solvers{self.CANDIDATE_LN_FLUID};
        else
            line_error(mfilename,'No LayeredNetwork solver available');
        end
end
end