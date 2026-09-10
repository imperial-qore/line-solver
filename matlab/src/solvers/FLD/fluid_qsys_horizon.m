function [bool, reason, t0, t1] = fluid_qsys_horizon(options)
% [BOOL, REASON, T0, T1] = FLUID_QSYS_HORIZON(OPTIONS)
%
% @brief The integration window of a time-varying single-station fluid limit.
%
% The time-varying limits ('mol', 'mtginf', 'tvms') report a TRAJECTORY, so an
% infinite or absent upper end of options.timespan leaves them nothing to
% report; the fluid analyzer applies the same rule to its own timespan.
%
% A horizon is a solver OPTION and not a model feature, so the feature
% registry has no name for it and SolverFLD.supportsModelMethod has to ask
% this predicate directly. SOLVER_FLUID_QSYS_ANALYZER asks the same one on the
% solve path, which is what keeps the report and the run from disagreeing
% about whether a method can be asked for.
%
% @param options solver options carrying the timespan
% @return bool true when the window is a finite non-empty interval
% @return reason the refusal, or '' when BOOL is true
% @return t0 lower end of the window (0 when unset)
% @return t1 upper end of the window (meaningful only when BOOL is true)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

t0 = 0; t1 = 1;
if isfield(options,'timespan') && numel(options.timespan) >= 2
    if isfinite(options.timespan(1))
        t0 = options.timespan(1);
    end
    t1 = options.timespan(2);
end
if ~isfinite(t1) || t1 <= t0
    bool = false;
    reason = ['A time-varying fluid method needs a finite horizon: set ' ...
        'options.timespan = [t0 t1].'];
    return
end
bool = true;
reason = '';
end
