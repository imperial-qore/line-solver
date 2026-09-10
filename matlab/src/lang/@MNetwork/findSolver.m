function T = findSolver(self, metric, showAll)
% T = FINDSOLVER(METRIC, SHOWALL)
%
% Which solvers and solver methods can analyze THIS model.
%
%   model.findSolver()                  every (solver, method) pair that runs
%   model.findSolver('cdf')             ... that returns a passage-time law
%   model.findSolver('getCdfRespT')     the same question, asked by accessor
%   model.findSolver('', true)          also the pairs that are refused, and why
%
% The returned table has one row per pair, with columns Solver, Method,
% Runnable, Class ('exact', 'approx', 'bound' or 'simulation'), Metrics and
% Reason. Method is the method name to pass as a solver method, so a row can be
% acted on directly:
%
%   T = model.findSolver('cdf');
%   solver = LINE(model, T.Method{1});
%
% FINDMETHOD and HELP are aliases of this method.
%
% See also SolverAUTO.findSolver, SolverAUTO.listValidMethods,
% MNetwork.getUsedLangFeatures
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    metric = '';
end
if nargin < 3
    showAll = false;
end
% The gate lives in SolverAUTO, which is the class that already knows every
% family, how to build one and what each refuses. Asking it here rather than
% reimplementing the walk is what keeps the model's answer and AUTO's own
% dispatch from being two opinions.
% The guard covers the CONSTRUCTION as well as the walk: SolverAUTO probes
% every candidate with supports(model) as it builds them, which warns on a
% model one of them refuses, and a report must not print.
verboseGuard = GlobalConstants.pushVerbose(VerboseLevel.SILENT); %#ok<NASGU>
T = SolverAUTO(self, 'verbose', 0).findSolver(metric, showAll);
end
