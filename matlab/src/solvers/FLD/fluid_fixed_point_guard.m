function fcn = fluid_fixed_point_guard(ode_h, fpRate)
% FCN = FLUID_FIXED_POINT_GUARD(ODE_H, FPRATE)
%
% An ODE OutputFcn that halts a window whose state has SETTLED MID-FLIGHT.
%
% WHY IT EXISTS. SOLVER_FLUID_ITERATION already short-circuits a window ENTERED
% at its fixed point, in closed form, because a stiff step controller handed a
% state it is already at cannot pick a step. That test is taken once, at the
% window's first instant, and a window that reaches the fixed point AFTER its
% first step is left to grind out the rest of its span at the reciprocal of the
% FASTEST rate while the span runs to 10*iter/min_rate, set by the SLOWEST.
%
% WHY IT IS AN OUTPUTFCN AND NOT AN EVENTS FUNCTION. The obvious slot is
% odeset('Events'), which costs nothing when it does not fire and is free of the
% conservation guard. It is the wrong one: the default accurateStiffOdeSolver is
% @lsoda_accurate_stiff, and LSODA_ODESOLVE implements the OutputFcn init/step/
% done protocol but has NO Events support -- so an Events-based check would work
% under @ode15s, do NOTHING on the default stiff path, and fail silently there.
% FLUID_OUTPUTFCN_CHAIN composes this guard with FLUID_CONSERVATION_GUARD so the
% one slot both integrators honour can carry both tests.
%
% THE THRESHOLD IS ROUND-OFF, NOT THE SOLVER TOLERANCE, and it is the same
% GlobalConstants.Zero the entry test uses, deliberately: a state that merely
% satisfies the window loop's driftTol (1e-4) is still moving, and cutting the
% window there was measured to shift results by 1.8e-5. Below Zero the drift is
% zero to double precision, so holding the state for the rest of the span is
% exact rather than approximate.
%
% THE CALLER MUST STILL PAD THE TAIL. Halting through an OutputFcn returns a
% SHORT final time, which SOLVER_FLUID_ITERATION reads as the conservation guard
% and turns into 'LINE:FluidNonHyperbolic'. LOCAL_SETTLE_TAIL re-tests the
% residual at the returned end state and, when it is settled, extends the
% trajectory to the requested instants at that state -- so a settled window is
% indistinguishable from one that ran to its endpoint, and only a genuine
% conservation excursion still reads as short.
%
% Parameters:
%   ode_h  - the drift, called as ode_h(t, y); must be the SAME handle the
%            window is integrating, and AUTONOMOUS, which is what the caller's
%            arming condition (no RT_BREAKS) already guarantees
%   fpRate - the slowest exit rate, normalising the residual to a rate ratio so
%            the test is scale-free in time as well as in population
%
% Returns:
%   fcn - handle for odeset('OutputFcn', ...). Returns status 1 to halt.
%
% See also SOLVER_FLUID_ITERATION, FLUID_CONSERVATION_GUARD,
% FLUID_OUTPUTFCN_CHAIN.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

fcn = @guard;

    function status = guard(t, y, flag)
        status = 0;
        % 'init' is the entry state, which the caller's own closed-form test
        % has already read at the same threshold, so re-testing it here would
        % halt a window before its first step for no gain. 'done' carries no
        % state at all.
        if isempty(y) || strcmp(flag, 'init') || strcmp(flag, 'done')
            return
        end
        for c = 1:size(y, 2)
            yc = y(:, c);
            ytot = sum(yc);
            if ytot <= 0
                continue
            end
            tc = t(min(c, numel(t)));
            if norm(ode_h(tc, yc), 1) / 2 / ytot / fpRate < GlobalConstants.Zero
                status = 1;
                return
            end
        end
    end
end
