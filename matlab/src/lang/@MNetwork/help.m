function T = help(self, metric, showAll)
% T = HELP(METRIC, SHOWALL)
%
% Alias of FINDSOLVER: what can this model be solved with?
%
% It SHADOWS the builtin HELP for Network objects, which is the point: an
% object handed to HELP is a request about that model, not about its class,
% and MATLAB's own answer to help(model) -- the Network class documentation --
% is reached by name as `help Network`, which this does not affect.
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
