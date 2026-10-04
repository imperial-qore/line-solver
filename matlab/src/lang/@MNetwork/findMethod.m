function T = findMethod(self, metric, showAll)
% T = FINDMETHOD(METRIC, SHOWALL)
%
% Alias of FINDSOLVER: which solvers and solver methods can analyze this
% model. The two names exist because the question is asked both ways round --
% "which solver do I use" and "which method do I pass" -- and the answer is
% the same table, whose Method column carries the method name either caller needs.
%
% See also MNetwork.findSolver
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    metric = '';
end
if nargin < 3
    showAll = false;
end
T = self.findSolver(metric, showAll);
end
