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

    % Handle tspan: if only [t0 tf], pass directly to MEX for adaptive
    % stepping (itask=2 mode collects all internal solver steps).
    % If >2 elements, integrate to each prescribed output time (itask=1).
    if length(tspan) == 2
        tspan = tspan(:)';
    else
        tspan = tspan(:)';
    end

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

    % Build MEX if needed
    thisDir = fileparts(mfilename('fullpath'));
    mexName = ['lsoda_mex.' mexext];
    if ~exist(fullfile(thisDir, mexName), 'file')
        fprintf('lsoda_mex not found, building...\n');
        oldDir = cd(thisDir);
        try
            build_lsoda_mex();
        catch ME
            cd(oldDir);
            rethrow(ME);
        end
        cd(oldDir);
    end

    % Call the MEX function
    [T, Y] = lsoda_mex(odefun, tspan, y0, rtol, atol, mxstep, mxordn, mxords);
end
