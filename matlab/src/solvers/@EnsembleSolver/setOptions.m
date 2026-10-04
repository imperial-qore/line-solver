function self = setOptions(self, eoptions)
% SELF = SETOPTIONS(EOPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

defaultOptions = self.defaultOptions;
if ~isfield(eoptions,'init_sol')
    eoptions.init_sol = defaultOptions.init_sol;
end
if ~isfield(eoptions,'iter_max')
    eoptions.iter_max = defaultOptions.iter_max;
end
if ~isfield(eoptions,'iter_tol')
    eoptions.iter_tol = defaultOptions.iter_tol;
end
if ~isfield(eoptions,'tol')
    eoptions.tol = defaultOptions.tol;
end
if ~isfield(eoptions,'verbose')
    eoptions.verbose = defaultOptions.verbose;
end
self.options = eoptions;
end
