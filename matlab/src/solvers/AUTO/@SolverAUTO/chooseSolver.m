function solver = chooseSolver(self, method)
% SOLVER = CHOOSESOLVER(METHOD)

line_debug('AUTO chooseSolver: method=%s, options.method=%s', method, self.options.method);

switch self.options.method
    case 'default'
        line_debug('Default method: using heuristic solver selection\n');
        line_debug('Using heuristic solver selection');
        solver = chooseSolverHeur(self, method);
    % case 'ai'  % AI method not yet available
    %     solver = chooseSolverAI(self, method);
    otherwise
        line_debug('Using heuristic solver selection (fallback)');
        solver = chooseSolverHeur(self, method);
end

line_debug('AUTO chooseSolver selected: %s', solver.getName());
end
