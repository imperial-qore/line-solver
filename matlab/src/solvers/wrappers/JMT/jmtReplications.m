function nreps = jmtReplications(options)
% NREPS = JMTREPLICATIONS(OPTIONS)
%
% The size of SolverJMT's transient ensemble: the number of independent JSIM
% runs, seeded seed, seed+1, ..., that runAnalyzer (method 'default' on a
% finite timespan) and getTranProbAggr average on a common time grid.
%
% Read as SolverLDES reads the same knob: a flat options.replications when set,
% else options.config.replications, else 10. Solver.parseOptions
% replaces the defaults wholesale when the caller passes its own struct, so
% neither spelling can be assumed present.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

nreps = 10;
if isfield(options, 'replications') && ~isempty(options.replications)
    nreps = options.replications;
elseif isfield(options, 'config') && isstruct(options.config) ...
        && isfield(options.config, 'replications') && ~isempty(options.config.replications)
    nreps = options.config.replications;
end
if ~isscalar(nreps) || ~isnumeric(nreps) || ~isfinite(nreps) || nreps < 1 || nreps ~= round(nreps)
    line_error(mfilename, sprintf(['options.config.replications is the number of independent ' ...
        'JSIM replications and must be a positive integer (got %s).'], mat2str(nreps)));
end
end
