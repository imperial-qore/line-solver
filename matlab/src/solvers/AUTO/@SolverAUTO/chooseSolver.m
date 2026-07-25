function solver = chooseSolver(self, method)
% SOLVER = CHOOSESOLVER(METHOD)

line_debug('AUTO chooseSolver: method=%s, selectionMode=%s', method, self.selectionMode);

switch self.selectionMode
    case 'default'
        line_debug('Default method: using heuristic solver selection\n');
        line_debug('Using heuristic solver selection');
        solver = chooseSolverHeur(self, method);
    case 'tree'
        % Learned selection refines the heuristic, it does not replace it: the
        % tree abstains outside its training support and the heuristic remains
        % the floor. See chooseSolverTree.m for the abstention thresholds.
        solver = chooseSolverLearned(self, method);
    otherwise
        line_debug('Using heuristic solver selection (fallback)');
        solver = chooseSolverHeur(self, method);
end

line_debug('AUTO chooseSolver selected: %s', solver.getName());
end

function solver = chooseSolverLearned(self, method)
% SOLVER = CHOOSESOLVERLEARNED(SELF, METHOD)

solver = [];
if isa(self.model, 'Network')
    choice = chooseSolverTree(self.model.getStruct());
    if ~isempty(choice)
        slot = candidateOfName(self, choice);
        % Feasibility gates the learned choice: a mispredicting selector must
        % never reach a solver whose feature set rejects the model.
        if ~isempty(slot) && ~isempty(self.solvers{slot}) && ...
                self.solvers{slot}.supports(self.model)
            line_debug('AUTO learned selection: %s', choice);
            solver = self.solvers{slot};
        else
            line_debug('AUTO learned selection %s rejected by feature set', choice);
        end
    else
        line_debug('AUTO learned selection abstained');
    end
end
if isempty(solver)
    solver = chooseSolverHeur(self, method);
end
end

function slot = candidateOfName(self, name)
% SLOT = CANDIDATEOFNAME(SELF, NAME)
%
% Map an algsel solver key onto a SolverAUTO candidate slot. The exact and
% approximate variants of MVA and NC share a slot: the harness distinguishes
% them by method, the candidate list does not.

switch name
    case {'mva', 'mva_exact'}
        slot = self.CANDIDATE_MVA;
    case {'nc', 'nc_exact'}
        slot = self.CANDIDATE_NC;
    case 'ctmc'
        slot = self.CANDIDATE_CTMC;
    case 'fld'
        slot = self.CANDIDATE_FLUID;
    case 'mam'
        slot = self.CANDIDATE_MAM;
    case 'ssa'
        slot = self.CANDIDATE_SSA;
    case 'jmt'
        slot = self.CANDIDATE_JMT;
    case 'ldes'
        slot = self.CANDIDATE_LDES;
    otherwise
        slot = [];
end
end
