function tf = fluid_hide_immediate(sn, options)
% TF = FLUID_HIDE_IMMEDIATE(SN, OPTIONS)
%
% Whether the fluid drift of this model should be built on the stochastic
% complement of its INSTANTANEOUS coordinates (see ODE_ELIMINATE_IMMEDIATE).
%
% Every fluid route that builds its drift from the station/class/phase event set
% asks here rather than reading OPTIONS.CONFIG.HIDE_IMMEDIATE directly, so that
% the answer is the same one across 'matrix', 'closing', 'statedep', 'tbi',
% 'minnormal', 'refined' and 'dae'. The flag defaults to TRUE for SolverFLD: a
% coordinate whose exit rate is GlobalConstants.Immediate is LINE's stand-in for
% infinity, and integrating it is meaningless work no integrator does well.
%
% THE STOCHASTIC PETRI NET ROUTE IS THE ONE EXCEPTION, and it is not a refusal.
% SOLVER_FLUID_PETRI carries immediate firings as ALGEBRAIC unknowns of an
% index-1 DAE (FLUID_PETRI_IMMEDIATE), which is a stronger treatment than
% absorbing them: it keeps the firing flow itself as a solved quantity rather
% than folding it into the timed events. That route never builds the event set
% this reduction acts on, so the answer here is simply FALSE and the flag is
% left alone rather than being turned into an error the way it once was.
%
% Parameters:
%   sn      - NetworkStruct
%   options - solver options
%
% Returns:
%   tf - true when the immediate coordinates are to be complemented away
%
% See also ODE_ELIMINATE_IMMEDIATE, SOLVER_FLUID_ODES, FLUID_MOMENT_TERMS.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

tf = false;
if ~isfield(options,'config') || ~isfield(options.config,'hide_immediate') ...
        || ~options.config.hide_immediate
    return
end
if isfield(sn,'nodetype') && any(sn.nodetype == NodeType.Transition)
    return
end
tf = true;
end % fluid_hide_immediate
