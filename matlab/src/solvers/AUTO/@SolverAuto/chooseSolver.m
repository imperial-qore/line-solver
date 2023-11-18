function solver = chooseSolver(self, method)
% SOLVER = CHOOSESOLVER(METHOD)

switch self.options.method
    case 'default'
        solver = chooseSolverHeur(self, method);
    case 'ai'
        solver = chooseSolverAI(self, method);
    case 'nn'
        solver = chooseSolverNN(self, method);
end
end
