function tf = lineTimeoutExceeded(options)
% TF = LINETIMEOUTEXCEEDED(OPTIONS)
% Cooperative wall-clock checkpoint for iterative solvers. Returns true if the
% elapsed time since the solver launch (options.timeout_tic, a tic handle set at
% the top of the per-solver runAnalyzer) exceeds options.timeout (seconds).
% A missing/non-finite/non-positive budget, or a missing tic, means no budget.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

tf = false;
if nargin == 0
    % Global-deadline form: deep utility functions (e.g. multichoose) that
    % have no options argument consult the session-level deadline set by the
    % per-solver runAnalyzer (see SolverCTMC). Empty appdata means no budget.
    d = getappdata(0, 'LINEtimeoutDeadline');
    if isempty(d)
        return;
    end
    tf = toc(d.tic) > d.budget;
    return;
end
if ~isstruct(options) || ~isfield(options, 'timeout') || ~isfield(options, 'timeout_tic')
    return;
end
budget = options.timeout;
if isempty(budget) || ~isfinite(budget) || budget <= 0
    return;
end
tf = toc(options.timeout_tic) > budget;
end
