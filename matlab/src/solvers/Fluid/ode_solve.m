function [t, yt_e] = ode_solve(ode_h, trange, y0, ode_opt, options)
% [T, YT_E] = SOLVEODE(Y0)

if ode_opt.AbsTol <= GlobalConstants.CoarseTol
    if GlobalConstants.Verbose == VerboseLevel.DEBUG
        if trange(1)==0
            line_printf('ode_solve: running accurateOdeSolver');
        end
    end
    [t, yt_e] = feval(options.odesolvers.accurateOdeSolver, ode_h, trange, y0, ode_opt);
else
    if GlobalConstants.Verbose == VerboseLevel.DEBUG
        if trange(1)==0
            line_printf('ode_solve: running fastOdeSolver');
        end
    end
    [t, yt_e] = feval(options.odesolvers.fastOdeSolver, ode_h, trange, y0, ode_opt);
end
end

