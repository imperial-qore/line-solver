function [T, Y] = lsoda_solve(odefun, tspan, y0, options)
% LSODA_SOLVE  Solve ODE system using LSODA (auto stiff/nonstiff switching)
%
%   [T, Y] = lsoda_solve(odefun, tspan, y0)
%   [T, Y] = lsoda_solve(odefun, tspan, y0, options)
%
%   Solves the ODE system dy/dt = odefun(t,y) from tspan(1) to tspan(end).
%   LSODA automatically switches between nonstiff (Adams) and stiff (BDF)
%   methods based on the problem behavior.
%
%   Inputs:
%     odefun  - Function handle @(t,y) returning column vector dy/dt
%     tspan   - Time span [t0 tf] or vector of output times [t0 t1 ... tf]
%     y0      - Initial conditions (column vector)
%     options - (optional) struct with fields:
%               .RelTol  - Relative tolerance (scalar or vector, default 1e-6)
%               .AbsTol  - Absolute tolerance (scalar or vector, default 1e-9)
%               .MaxStep - Maximum internal steps per interval (default 5000)
%               .MaxOrdNonStiff - Max order for Adams method (1-12, default 12)
%               .MaxOrdStiff    - Max order for BDF method (1-5, default 5)
%               .Backend - 'mex' or 'matlab'; by default the compiled MEX is
%                          used where present and lsoda_matlab.m otherwise
%
%   Outputs:
%     T - Column vector of output times
%     Y - Solution matrix (length(T) x length(y0))
%
%   Example:
%     f = @(t,y) [-0.04*y(1) + 1e4*y(2)*y(3);
%                  0.04*y(1) - 1e4*y(2)*y(3) - 3e7*y(2)^2;
%                  3e7*y(2)^2];
%     [T, Y] = lsoda_solve(f, [0 4e10], [1; 0; 0]);
%
%   Based on LSODA by L.R. Petzold and A.C. Hindmarsh (LLNL).
%   C implementation from liblsoda (MIT License).

    % Ensure y0 is a column vector
    y0 = y0(:);

    % [t0 tf] runs itask=2 (every internal step is collected), more than two
    % entries runs itask=1 (integrate to each prescribed output time).
    tspan = tspan(:)';

    % Default options
    rtol = 1e-6;
    atol = 1e-9;
    mxstep = 5000;
    mxordn = 0;  % 0 = use library default (12)
    mxords = 0;  % 0 = use library default (5)

    if nargin >= 4 && ~isempty(options)
        if isfield(options, 'RelTol'), rtol = options.RelTol; end
        if isfield(options, 'AbsTol'), atol = options.AbsTol; end
        if isfield(options, 'MaxStep'), mxstep = options.MaxStep; end
        if isfield(options, 'MaxOrdNonStiff'), mxordn = options.MaxOrdNonStiff; end
        if isfield(options, 'MaxOrdStiff'), mxords = options.MaxOrdStiff; end
    end

    % Select the backend: the compiled MEX when it is available, the pure
    % MATLAB port otherwise. options.Backend forces one of 'mex', 'matlab'.
    backend = '';
    if nargin >= 4 && ~isempty(options) && isfield(options, 'Backend')
        backend = lower(options.Backend);
    end

    thisDir = fileparts(mfilename('fullpath'));
    mexName = ['lsoda_mex.' mexext];
    haveMex = exist(fullfile(thisDir, mexName), 'file') == 3;

    if strcmp(backend, 'matlab')
        [T, Y] = lsoda_matlab(odefun, tspan, y0, rtol, atol, mxstep, mxordn, mxords);
        return
    end

    if ~haveMex
        oldDir = cd(thisDir);
        try
            build_lsoda_mex();
            haveMex = exist(fullfile(thisDir, mexName), 'file') == 3;
        catch ME
            if strcmp(backend, 'mex')
                cd(oldDir);
                rethrow(ME);
            end
        end
        cd(oldDir);
    end

    if haveMex
        [T, Y] = lsoda_mex(odefun, tspan, y0, rtol, atol, mxstep, mxordn, mxords);
    else
        [T, Y] = lsoda_matlab(odefun, tspan, y0, rtol, atol, mxstep, mxordn, mxords);
    end
end
