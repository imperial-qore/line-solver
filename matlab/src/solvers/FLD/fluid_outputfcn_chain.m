function fcn = fluid_outputfcn_chain(fcns)
% FCN = FLUID_OUTPUTFCN_CHAIN(FCNS)
%
% Compose several ODE OutputFcns into the ONE slot odeset provides.
%
% WHY IT EXISTS. `odeset` carries a single 'OutputFcn', and the fluid solver has
% two independent reasons to halt a window: FLUID_CONSERVATION_GUARD, when the
% moment-closure drift has left the model, and FLUID_FIXED_POINT_GUARD, when the
% state has settled and the rest of the span is known in closed form. Neither is
% a special case of the other and both must run on the SAME integrator, so the
% slot holds a chain rather than a choice. Empty entries are dropped, so a
% caller may pass a slot that was never filled.
%
% MATLAB's protocol is honoured for every member: 'init' first, '' after each
% accepted step, 'done' last. A halt is the disjunction -- ANY member returning
% nonzero halts -- but every member is still called on that step, because a
% member may be keeping state and skipping it would corrupt it. 'done' always
% reports 0: it is a teardown notification, not a decision point, and MATLAB
% ignores the status there.
%
% Parameters:
%   fcns - cell array of OutputFcn handles, called as fcn(t, y, flag). Empty
%          elements are ignored, so {[], guard} is a chain of one
%
% Returns:
%   fcn - handle for odeset('OutputFcn', ...), or [] when nothing is left after
%         dropping the empties, which the caller may store back into the slot
%         unchanged
%
% See also SOLVER_FLUID_ITERATION, FLUID_CONSERVATION_GUARD,
% FLUID_FIXED_POINT_GUARD.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

keep = false(1, numel(fcns));
for i = 1:numel(fcns)
    keep(i) = ~isempty(fcns{i});
end
fcns = fcns(keep);
if isempty(fcns)
    fcn = [];
    return
end
if numel(fcns) == 1
    fcn = fcns{1};
    return
end

fcn = @chained;

    function status = chained(t, y, flag)
        status = 0;
        for i = 1:numel(fcns)
            s = fcns{i}(t, y, flag);
            if ~isempty(s) && s ~= 0
                status = 1;
            end
        end
        if strcmp(flag, 'done')
            status = 0;
        end
    end
end
