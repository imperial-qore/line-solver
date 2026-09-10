function QN = solveMeansForStruct(self, sn)
% QN = SOLVEMEANSFORSTRUCT(SN) Mean queue lengths of a perturbed structure,
% under this solver's own method.
%
% This is the mean-value oracle behind the numerical-derivative path of
% getMomentTable, getMomentChainTable and getMomentStationTable. The moment
% identity Cov[n,n] = L dQ/dL does not care HOW the mean queue lengths were
% obtained, only that they are the means of a product-form model as a function
% of its demands. So rather than hand-differentiating each algorithm, the solver
% is re-run on a perturbed structure and differentiated numerically.
%
% Running THIS solver, rather than one chosen algorithm, is what makes the path
% general: it covers every method of every solver, including the
% normalizing-constant methods of SolverNC (comom, ca, le, mom, ...) and the
% summation methods of SolverMVA (sum, esum), none of which any
% hand-differentiated implementation reaches. Restricting the oracle to one
% analyzer would restrict the moments to that analyzer's methods for no
% mathematical reason.
%
% The perturbed SN is injected by copying the model and overwriting its cached
% structure, then constructing a fresh solver of this class with this solver's
% options. The copy matters: Network is a handle object (Model < Copyable <
% handle), so writing sn on the caller's model would corrupt it, and the fresh
% solver matters because a solver caches its results and would otherwise return
% the unperturbed answer.
%
% See also: getMomentStationTable, getMomentChainTable, getMomentTable.

model = self.model.copy;
model.sn = sn;
model.hasStruct = true;   % the copy must not regenerate and discard the perturbation

solver = feval(class(self), model, self.getOptions);
QN = solver.getAvgQLen();
end
