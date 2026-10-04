function solver = chooseSolver(self, method)
% SOLVER = CHOOSESOLVER(METHOD)

line_debug('AUTO chooseSolver: method=%s, selectionMode=%s', method, self.selectionMode);

switch self.selectionMode
    case {'default', 'heur'}
        line_debug('Using heuristic solver selection');
        solver = chooseSolverHeur(self, method);
    case 'exact'
        line_debug('Using exact solver selection');
        solver = chooseSolverExact(self, method);
    case 'sim'
        line_debug('Using simulation solver selection');
        solver = chooseSolverSim(self, method);
    case 'fast'
        % Cheapest analytical answer; the heuristic is the floor for metrics
        % no mean-value solver can serve. The candidate slot ids are per model
        % class, so only Network models take the ranked shortcut.
        solver = [];
        if isa(self.model,'Network')
            solver = chooseSolverRanked(self, [self.CANDIDATE_MVA, self.CANDIDATE_NC, ...
                self.CANDIDATE_FLUID, self.CANDIDATE_MAM]);
        end
        if isempty(solver)
            solver = chooseSolverHeur(self, method);
        end
    case 'accurate'
        % Smooth or matrix-analytic answer preferred over the fastest one.
        solver = [];
        if isa(self.model,'Network')
            solver = chooseSolverRanked(self, [self.CANDIDATE_FLUID, self.CANDIDATE_MAM, ...
                self.CANDIDATE_CTMC, self.CANDIDATE_LDES]);
        end
        if isempty(solver)
            solver = chooseSolverHeur(self, method);
        end
    otherwise
        line_debug('Using heuristic solver selection (fallback)');
        solver = chooseSolverHeur(self, method);
end

line_debug('AUTO chooseSolver selected: %s', solver.getName());
end
