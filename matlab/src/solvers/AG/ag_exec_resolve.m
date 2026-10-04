function exec = ag_exec_resolve(options)
% EXEC = AG_EXEC_RESOLVE(OPTIONS)
%
% Resolve the execution backend of the RCAT fixed point into a struct the
% sweep can act on without re-reading options.
%
% THE BACKEND CHANGES WHO EVALUATES AN AGENT, NEVER WHAT IT EVALUATES TO.
% Agent k's generator is Q_k(x) = L_k + sum_{c passive at k} x_c Pb_c, which
% depends on the other agents only through the scalar reversed rates x. The
% sweep is Jacobi -- every x_a is read off the PREVIOUS sweep's stationary
% vectors, then all agents re-solve -- so the agent order is immaterial and a
% parallel or distributed sweep produces the SAME iterates as the serial one.
%
% HOW FAR THAT SURVIVES FLOATING POINT DEPENDS ON WHO RUNS THE AGENT SOLVE.
% 'parallel' is bit-identical to 'serial': the same code in the same process,
% differing only in an order that does not matter. 'cluster' is bit-identical
% only when the worker runs the same implementation as the coordinator -- the
% wire is exact, since JSON round-trips a double without loss, but the
% stationary vector comes back from the WORKER's solve, so a MATLAB coordinator
% driving a Java ag-worker agrees to a few ulp rather than bit for bit. That is
% the ordinary cross-codebase difference, not a protocol defect, but a cluster
% run must not regenerate a golden that a same-language run will later read.
%
% Fields:
%   mode       'serial' | 'parallel' (alias 'para') | 'cluster'
%   nworkers   'parallel': pool size, 0 = MATLAB's default pool
%   endpoints  'cluster': cell array of 'host:port' worker addresses
%   timeout    'cluster': seconds to wait on a worker before solving its
%              agents locally instead
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

exec = struct('mode', 'serial', 'nworkers', 0, 'endpoints', {{}}, 'timeout', 30);

if nargin < 1 || isempty(options) || ~isfield(options, 'config') || isempty(options.config)
    return;
end
cfg = options.config;

if isfield(cfg, 'exec') && ~isempty(cfg.exec)
    mode = lower(char(cfg.exec));
    % 'para' is the accepted alias of 'parallel', the same pair SolverSSA answers
    % to (solver_ssa_analyzer.m: case {'para','parallel'}). It is NORMALISED here
    % so everything downstream -- the switch in solver_ag.m, the pool check below
    % -- compares against one spelling only.
    if strcmp(mode, 'para')
        mode = 'parallel';
    end
    if ~any(strcmp(mode, {'serial','parallel','cluster'}))
        % 'threads' was this backend's name until 2026-08-19. Name the rename
        % rather than reporting a backend that still exists as unknown.
        if strcmp(mode, 'threads')
            line_error(mfilename, ['The ''threads'' execution backend was renamed to ' ...
                '''parallel'' (alias ''para'').']);
        end
        line_error(mfilename, ['Unknown AG execution backend ''%s''. Use ''serial'', ' ...
            '''parallel'' (alias ''para'') or ''cluster''.'], mode);
    end
    exec.mode = mode;
end
if isfield(cfg, 'nworkers') && ~isempty(cfg.nworkers)
    exec.nworkers = cfg.nworkers;
end
if isfield(cfg, 'endpoints') && ~isempty(cfg.endpoints)
    exec.endpoints = cfg.endpoints;
end
if isfield(cfg, 'worker_timeout') && ~isempty(cfg.worker_timeout)
    exec.timeout = cfg.worker_timeout;
end

if strcmp(exec.mode, 'cluster') && isempty(exec.endpoints)
    line_error(mfilename, ['The ''cluster'' execution backend needs worker endpoints: set ' ...
        'options.config.endpoints to a cell array of ''host:port'' strings, each one an ' ...
        'ag-worker started with ''java -cp jline.jar jline.cli.AgWorker -p <port>''.']);
end

if strcmp(exec.mode, 'parallel')
    % A pool is REQUIRED rather than optional: parfor without one silently
    % degrades to a serial loop, which would make a 'parallel' run report a
    % backend it did not use. Opening it here also keeps the cost out of the
    % sweep, which would otherwise pay it on the first iteration.
    if isempty(ver('parallel'))
        line_error(mfilename, ['The ''parallel'' execution backend needs the Parallel ' ...
            'Computing Toolbox. Use ''serial'', or ''cluster'' with ag-worker processes.']);
    end
    pool = gcp('nocreate');
    if isempty(pool)
        if exec.nworkers > 0
            parpool(exec.nworkers);
        else
            parpool;
        end
    end
end

end
