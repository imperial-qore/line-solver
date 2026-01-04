function [t, yt_e] = ode_solve_stiff(ode_h, trange, y0, ode_opt, options)
% [T, YT_E] = SOLVEODESTIFF(Y0)

if options.stiff
    [t, yt_e] = feval(options.odesolvers.accurateStiffOdeSolver, ode_h, trange, y0, ode_opt);
else
    ode_opt.NonNegative = []; % not supported by ode23s
    [t, yt_e] = feval(options.odesolvers.fastStiffOdeSolver, ode_h, trange, y0, ode_opt);
end
end
