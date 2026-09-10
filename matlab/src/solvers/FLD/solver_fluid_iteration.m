function [xvec_it, xvec_t, t, iter] = solver_fluid_iteration(sn, N, Mu, Phi, PH, P, S, xvec_it, ydefault, slowrate, Tstart, max_time, options)
% [XVEC_IT, XVEC_T, T, ITER] = SOLVER_FLUID_ITERATION(QN, N, MU, PHI, PH, P, S, YMEAN, YDEFAULT, SLOWRATE, TSTART, MAX_TIME, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

iter_max = options.iter_max;
verbose = options.verbose;
tol = options.tol;
iter_tol = options.iter_tol;
stiff = options.stiff;
timespan = options.timespan;

goon = true; % max stiff solver
iter = 0;
consoleSettled = false; % console: the mass stopped moving, stop reporting
LineConsole.loop('integrating the fluid ODEs over successive time windows');
lastmsg = '';
t=[];
xvec_t=[];
% per-pass chunks, concatenated once at the end: growing xvec_t and t in place
% recopies the whole history on every pass, which is quadratic in iter_max
tchunk = cell(1,iter_max);
xchunk = cell(1,iter_max);
% heuristic to select stiff or non-stiff ODE solver
nonZeroRates = slowrate(:);
nonZeroRates = nonZeroRates( nonZeroRates >tol );

nonZeroRates = nonZeroRates(isfinite(nonZeroRates));
if isempty(nonZeroRates)
    nonZeroRates = 1; % fallback when all rates are zero or infinite
end
%rategap = log10(max(nonZeroRates)/min(nonZeroRates)); % if the max rate is InfRate and the min is 1, then rategap = 6

% init ode
[ode_h, ~, rt_breaks, absorb] = solver_fluid_odes(sn, N, Mu, Phi, PH, P, S, sn.sched, sn.schedparam, options);

% THE INITIAL POINT HAS TO BE PROJECTED TOO. Once the instantaneous coordinates
% are complemented away no event moves them any more, so whatever mass the
% initial condition parked there -- SOLVER_FLUID_INITSOL starts everything in
% phase 1, but a warm start from an earlier LN iterate does not -- would be
% frozen for the whole integration and lost from its chain. ABSORB sends it
% where the eliminated coordinate would have sent it instantaneously.
if ~isempty(absorb)
    for ia = 1:numel(xvec_it)
        if ~isempty(xvec_it{ia})
            xvec_it{ia} = reshape(absorb' * xvec_it{ia}(:), size(xvec_it{ia}));
        end
    end
end

% Early stop on the DRIFT RESIDUAL (see the termination test below).
% Disabled with options.config.fluid_earlystop = false, which restores the
% unconditional iter_max windows.
earlystop = true;
if isfield(options, 'config') && isfield(options.config, 'fluid_earlystop')
    earlystop = options.config.fluid_earlystop;
end
% The residual cannot be driven below the error the integrator itself carries,
% so a caller asking for less than tol is asking for something unobservable.
driftTol = max(iter_tol, tol);
% A mode relaxing at the slowest exit rate needs about 10/minrate to show at
% all: the test is not consulted before the loop has integrated that far.
slowestRate = min(nonZeroRates);
minHorizon = 10 / slowestRate;
driftBelow = 0; % consecutive windows satisfying the residual test
movedMassRatioPrev = Inf; % previous window's moved mass, for the tail estimate
rhoHist = nan(1,3);   % recent contraction ratios; the worst one drives the tail
driftSafety = 0.01;   % headroom on the tail, since rho is estimated, not known

T0 = timespan(1);
odeopt = odeset('AbsTol', tol, 'RelTol', tol, 'NonNegative', 1:length(xvec_it{1}));
% ARM THE CONSERVATION GUARD ONLY WHERE THE MOMENT CLOSURE IS ACTIVE.
% SOLVER_FLUID_MOMENTS runs its first pass at sigma2 = 0 -- that IS the
% first-order solve -- and only the LATER passes integrate a drift that can
% leave the simplex. Testing sigma2 keeps the guard off every first-order call,
% so 'matrix' and 'closing' (which are also the ladder's own fallback) keep
% byte-identical behaviour and cannot be sent down a fallback by their own
% watchdog. See FLUID_CONSERVATION_GUARD.
closureActive = isfield(options, 'config') && isfield(options.config, 'moment_sigma2') ...
    && ~isempty(options.config.moment_sigma2) && any(options.config.moment_sigma2 ~= 0);
if closureActive
    % Same layout SOLVER_FLUID builds, recomputed from Mu rather than read off
    % sn.phases, because that is what fixes the coordinate blocks here.
    guardPhases = zeros(sn.nstations, sn.nclasses);
    for gi = 1:sn.nstations
        for gk = 1:sn.nclasses
            guardPhases(gi,gk) = length(Mu{gi}{gk});
        end
    end
    % ONLY WHEN THE STATE IS THE ONE THE PHASES DESCRIBE. options.config
    % .hide_immediate folds immediate coordinates out through
    % ODE_ELIMINATE_IMMEDIATE, which SHRINKS the vector, and the block offsets
    % read off phases would then address the wrong coordinates and trip on a
    % sum that was never that class's population. The lengths agreeing is the
    % exact test for that, so the guard stays off rather than guessing.
    if sum(guardPhases(:)) == length(xvec_it{1})
        if isfield(sn, 'chains') && ~isempty(sn.chains)
            guardChains = sn.chains;
        else
            guardChains = eye(sn.nclasses);
        end
        odeopt = odeset(odeopt, 'OutputFcn', ...
            fluid_conservation_guard(guardPhases, sn.njobs(:)', guardChains, 0.1));
    else
        closureActive = false;
    end
end
T = 0;
while (isfinite(timespan(2)) && T < timespan(2)) || (goon && iter < iter_max)
    iter = iter + 1;
    if toc(Tstart) > max_time
        goon = false;
        break;
    end
    
    % determine entry state vector in e
    y0 = xvec_it{iter-1 +1};
    
    if iter == 1 % first iteration
        T = min(timespan(2),abs(10/min(nonZeroRates))); % solve ode until T = 1 event with slowest exit rate
    else
        T = min(timespan(2),abs(10*iter/min(nonZeroRates)));
    end
    trange = [T0, T];
    % A CALLER MAY ASK FOR THE OUTPUT TIMES ITSELF. `options.tranpoints`, when
    % set, is an increasing vector of instants the trajectory is wanted at, and
    % the chunk integrates over exactly those (plus its own endpoints) instead
    % of over [T0, T]. MATLAB's ode solvers return precisely the vector they are
    % handed, so this costs no accuracy and no extra steps.
    %
    % It exists because SolverENV sums an exit average against a sojourn CDF:
    % over a horizon of 1e3 read through an Exp(1) clock, the integrator's own
    % grid puts almost every point where the weight is zero, and INTERPOLATING
    % it cannot recover resolution the trajectory never had. Unset, nothing
    % changes for any other caller.
    if isfield(options, 'tranpoints') && ~isempty(options.tranpoints)
        pts = options.tranpoints(:).';
        pts = pts(pts > T0 & pts < T);
        if ~isempty(pts)
            trange = unique([T0, pts, T]);
        end
    end

    % Armed only for an AUTONOMOUS drift: with a rate schedule (RT_BREAKS) a zero
    % residual now says nothing about the next segment, so the window is stepped.
    % It rides on EARLYSTOP because it reads the same residual that test does.
    % A FINITE timespan is a transient request, which must reach its end time
    % rather than stop at the fixed point -- the same gate the tail test carries.
    if earlystop && isempty(rt_breaks) && ~isfinite(timespan(2))
        fpRate = slowestRate;
    else
        fpRate = [];
    end
    try
        [t_iter, ymean_t_iter] = local_integrate(ode_h, trange, y0, odeopt, options, stiff, rt_breaks, fpRate);
    catch ME
        line_printf('\nThe initial point is invalid, Fluid solver switching to default initialization.');
        odeopt = odeset('AbsTol', tol, 'RelTol', tol, 'NonNegative', 1:length(ydefault));
        [t_iter, ymean_t_iter] = local_integrate(ode_h, trange, ydefault, odeopt, options, false, rt_breaks, fpRate);
    end
    % A SHORT FINAL TIME MEANS THE GUARD HALTED THE WINDOW. MATLAB's ode
    % solvers return the requested endpoint unless an OutputFcn stops them, so
    % this is the guard and nothing else. Raise the identifier the fallback
    % ladder in @SolverFLD/runAnalyzer.m already catches, so the model is
    % answered by 'dae' and then by 'matrix'/'closing' instead of hanging here.
    if closureActive && ~isempty(t_iter) && t_iter(end) < trange(end) - 1e-9*max(1,abs(trange(end)))
        throw(MException('LINE:FluidNonHyperbolic', ...
            ['[%s.m] The moment-closure drift left the model: a closed chain '...
             'lost more than 10%% of its population inside the window ending at '...
             't = %g, which the drift conserves exactly, so the excursion is the '...
             'NonNegative clamp injecting mass rather than a solution. Integrating '...
             'on would not return. Falling back to a first-order closure.'], ...
            mfilename, T));
    end
    xchunk{iter} = ymean_t_iter;
    tchunk{iter} = t_iter;
    xvec_it{iter +1} = ymean_t_iter(end,:);
    movedMassRatio = norm(xvec_it{iter +1} - xvec_it{iter-1 +1}, 1) / 2 / sum(xvec_it{iter-1 +1});
    % the loop deliberately runs to iter_max even after the state settles
    % (see below), so the console reports only the windows that still move
    if ~consoleSettled
        if movedMassRatio > 0
            LineConsole.iter(iter, ['window %d up to t = %g: moved mass %.3e, ' ...
                '%d ODE points'], iter, T, movedMassRatio, numel(t_iter));
        else
            consoleSettled = true;
            % through DETAIL, not STEP: the outer moment-closure loop restarts
            % this integration once per iterate and would repeat the line
            LineConsole.detail(sprintf(['fluid state settled at t = %g, integrating ' ...
                'on to the horizon without further change'], T));
        end
    end
    T0  = T; % for next iteration
    
    % check termination condition
    
    if verbose > 0
        llmsg = length(lastmsg);
        if llmsg>0
            for ib=1:llmsg
                line_printf('\b');
            end
        end
    end
    
    % MOVEDMASSRATIO IS NOT THE TERMINATION TEST. It is the mass moved over ONE
    % window, and it underestimates the distance still left to the fixed point
    % by exactly the geometric tail it drops: for a mode contracting by rho per
    % window the remaining distance is r*rho/(1-rho), which is what is tested
    % here. rho is read off the iteration itself, so no rate in the model has to
    % stand in for the slowest system mode -- on a slowly mixing model there is
    % none that can, which is why the drift norm alone stops 3% short.
    % The drift F(x), zero AT a fixed point, is kept as a second, independent
    % bound: BOTH must hold, on two consecutive windows, past the slowest
    % relaxation time. See _kb/06-solver-catalog.md.
    if earlystop && goon && T >= minHorizon && iter > 1
        % One ratio sees only the mode that dominates THESE two windows: with
        % several modes in play, and a window that lengthens as the loop runs,
        % it reads faster than the slowest one. Take the worst ratio still on
        % record, and hold the tail to a fraction of the tolerance, so what is
        % left when the loop stops sits below what the caller asked for.
        rhoHist(1 + mod(iter, numel(rhoHist))) = ...
            movedMassRatio / max(movedMassRatioPrev, GlobalConstants.Zero);
        rho = max(rhoHist(isfinite(rhoHist)));
        xend = xvec_it{iter +1}(:);
        driftDispl = norm(ode_h(T, xend), 1) / 2 / max(sum(xend), GlobalConstants.Zero) / slowestRate;
        % a non-contracting iteration has no tail to sum: it is not converging
        %
        % THE 1e-6 GATE IS NOT AN OVERSIGHT, EVEN THOUGH IT SITS BELOW THE
        % INTEGRATOR'S OWN tol. Relaxing it to "the moved mass reached the
        % integrator floor, so trust the drift residual alone" was TRIED and
        % REVERTED: it stops the M/M/1 rho = 0.9 minnormal solve at 7.018088
        % against the 7.021524680 all four codebases agree on
        % (MinNormalTest.testOpenMm1), and it ends the statedep trajectory of
        % test_exampleCdfRespT2StatedepMethod at t = 60 instead of its 2000.
        % The accuracy of this loop comes from running the windows, so the
        % stop has to stay conservative. It is also not what makes a solve
        % hang: see _kb/06-solver-catalog.md, where the minnormal closure
        % diverges outright on a bounded multiserver station.
        if rho < 1
            tailEstimate = movedMassRatio * rho / (1 - rho);
            if tailEstimate < driftSafety * driftTol && driftDispl < driftTol
                driftBelow = driftBelow + 1;
                if driftBelow >= 2
                    goon = false;
                    LineConsole.detail(sprintf(['fluid tail estimate %.3e and drift %.3e ' ...
                        'below %.3e at t = %g, stopping after %d windows'], ...
                        tailEstimate, driftDispl, driftTol, T, iter));
                end
            else
                driftBelow = 0;
            end
        else
            driftBelow = 0;
        end
    end
    movedMassRatioPrev = movedMassRatio;
    
    if T >= timespan(2)
        goon = false;
    end
end

if iter > 0
    xvec_t = vertcat(xchunk{1:iter});
    t = vertcat(tchunk{1:iter});
end

end

function [t_iter, y_iter] = local_integrate(ode_h, trange, y0, odeopt, options, stiff, rt_breaks, fpRate)
% [T_ITER, Y_ITER] = LOCAL_INTEGRATE(ODE_H, TRANGE, Y0, ODEOPT, OPTIONS, STIFF, RT_BREAKS, FPRATE)
%
% Integrate one window, RESTARTING THE INTEGRATOR AT EVERY JUMP OF THE DRIFT.
%
% FPRATE, when non-empty, is the slowest exit rate, and it arms the FIXED-POINT
% SHORT CIRCUIT below. See the block at the head of the function body.
%
% RT_BREAKS holds the instants where an NHPP schedule steps its intensity, so
% the right-hand side is piecewise constant with jumps there and nothing in it
% tells a step controller where they are. A step that straddles one integrates a
% rate the model never had over part of its span, and no error estimate catches
% that: the states either side are both smooth, so the controller reads the
% discontinuity as a well-resolved segment. Making each jump a boundary is the
% fix, and it costs nothing when there are none -- the common case, where this
% reduces to the single call it replaced.
%
% See also SOLVER_FLUID_RATEMULT, ODE_SOLVE_STIFF.

if nargin < 8
    fpRate = [];
end

% A FIXED POINT ENDS THE WINDOW IN CLOSED FORM, and this is what keeps a window
% that has already converged from becoming a window that never returns. The
% caller arms FPRATE only when RT_BREAKS is empty, i.e. when the drift is
% AUTONOMOUS, so F(y*) = 0 means y(t) = y* for every later t and the rest of the
% span is known exactly rather than integrated.
%
% THE THRESHOLD IS ROUND-OFF, NOT THE SOLVER TOLERANCE. The window loop's
% driftDispl < driftTol (1e-4 by default) says "converged to what the caller
% asked for", and a state that merely satisfies THAT is still moving -- cutting
% the window there was measured to shift results by 1.8e-5 in the Python twin. A
% normalized residual below GlobalConstants.Zero is the stronger claim that the
% drift is zero to double precision, and that is what makes skipping the rest of
% the span exact instead of approximate.
%
% WHY THE WINDOW DOES NOT END ON ITS OWN. A stiff step controller handed a state
% it is already at cannot pick a step: on the LN layer of test_LQN_13 the Python
% twin advanced t by 0.011 in 20000 steps from a state with |F| = 1.5e-16, and
% covered the whole 1000-unit span in 60 steps once that state was nudged 1e-6
% off the equilibrium. The layer carries an Immediate() coordinate --
% an eigenvalue of exactly -GlobalConstants.Immediate = -1e8 that
% ODE_ELIMINATE_IMMEDIATE did not fold out -- so the controller is pinned near
% 1/1e8 while the window runs to 10*iter/min(nonZeroRates). 272 windows of that
% layer took 3.0 s between them and the 273rd had not returned after 143 s.
if nargin >= 8 && ~isempty(fpRate) && fpRate > 0
    ytot = sum(y0);
    if ytot > 0 && norm(ode_h(trange(1), y0(:)), 1) / 2 / ytot / fpRate < GlobalConstants.Zero
        % TRANGE, not just its endpoints: a caller that asked for output
        % instants (options.tranpoints) gets them, at the settled state.
        t_iter = trange(:);
        y_iter = repmat(y0(:)', numel(t_iter), 1);
        return
    end
end

% THE ENTRY TEST ABOVE IS TAKEN ONCE, and a window that reaches the fixed point
% after its first step is left grinding out the rest of its span. The same test
% therefore also rides on the integrator, through the OutputFcn slot chained
% with whatever the caller already put there. Events would be the natural slot
% and is NOT usable: the default accurateStiffOdeSolver is LSODA, which honours
% OutputFcn and ignores Events, so an Events check would silently do nothing on
% the default stiff path. LOCAL_SETTLE_TAIL below turns the resulting short
% window back into a full one.
if nargin >= 8 && ~isempty(fpRate) && fpRate > 0
    odeopt = odeset(odeopt, 'OutputFcn', fluid_outputfcn_chain( ...
        {odeget(odeopt, 'OutputFcn', []), fluid_fixed_point_guard(ode_h, fpRate)}));
end

edges = [];
if ~isempty(rt_breaks)
    lo = min(trange(1), trange(end));
    hi = max(trange(1), trange(end));
    span = max(1, abs(hi - lo));
    edges = rt_breaks(rt_breaks > lo + 1e-12*span & rt_breaks < hi - 1e-12*span);
end
if isempty(edges)
    [t_iter, y_iter] = local_odecall(ode_h, trange, y0, odeopt, options, stiff);
    [t_iter, y_iter] = local_settle_tail(t_iter, y_iter, trange, ode_h, fpRate);
    return
end

% One sub-window per segment, each asked for the output instants the caller
% wanted inside it. The shared boundary row is emitted once, by the segment
% that starts there, so the concatenated trajectory carries no duplicate.
pts = unique([trange(:)', edges]);
t_iter = [];
y_iter = [];
ycur = y0;
for e = 1:numel(edges)+1
    if e == 1
        a = pts(1);
    else
        a = edges(e-1);
    end
    if e == numel(edges)+1
        b = pts(end);
    else
        b = edges(e);
    end
    sub = pts(pts > a & pts < b);
    if isempty(sub)
        subrange = [a, b];
    else
        subrange = [a, sub, b];
    end
    [tk, yk] = local_odecall(ode_h, subrange, ycur, odeopt, options, stiff);
    if isempty(tk)
        break
    end
    if isempty(t_iter)
        t_iter = tk(:);
        y_iter = yk;
    else
        t_iter = [t_iter; tk(2:end)]; %#ok<AGROW>
        y_iter = [y_iter; yk(2:end,:)]; %#ok<AGROW>
    end
    ycur = yk(end,:);
    % A SEGMENT THAT STOPPED SHORT IS THE CONSERVATION GUARD HALTING, and the
    % caller reads a short final time as exactly that. Integrating the later
    % segments would hide the halt behind a full-length window.
    if tk(end) < b - 1e-9*max(1,abs(b))
        break
    end
end
end

function [tk, yk] = local_odecall(ode_h, subrange, y0, odeopt, options, stiff)
if stiff
    [tk, yk] = ode_solve_stiff(ode_h, subrange, y0, odeopt, options);
else
    [tk, yk] = ode_solve(ode_h, subrange, y0, odeopt, options);
end
end

function [t_iter, y_iter] = local_settle_tail(t_iter, y_iter, trange, ode_h, fpRate)
% [T_ITER, Y_ITER] = LOCAL_SETTLE_TAIL(T_ITER, Y_ITER, TRANGE, ODE_H, FPRATE)
%
% Extend a window that FLUID_FIXED_POINT_GUARD halted out to the instants the
% caller asked for, holding the settled state.
%
% AN OUTPUTFCN HALT IS AMBIGUOUS and this is what disambiguates it. Two guards
% share the slot, and both report the same way: a final time short of the
% window's end. SOLVER_FLUID_ITERATION reads that as the CONSERVATION guard and
% raises 'LINE:FluidNonHyperbolic' so the fallback ladder answers the model. So
% the residual is re-tested at the returned end state, at the threshold the
% fixed-point guard used: below it the halt was a settled window and the rest of
% the span is x(t) = x*, which is emitted; at or above it the halt was the
% conservation guard and the short window is left exactly as it came back.
%
% See also FLUID_FIXED_POINT_GUARD, FLUID_CONSERVATION_GUARD.

if isempty(fpRate) || fpRate <= 0 || isempty(t_iter) || isempty(y_iter)
    return
end
tend = trange(end);
if t_iter(end) >= tend - 1e-9 * max(1, abs(tend))
    return % ran to its endpoint: nothing to extend
end
ylast = y_iter(end, :).';
ytot = sum(ylast);
if ytot <= 0
    return
end
if norm(ode_h(t_iter(end), ylast), 1) / 2 / ytot / fpRate >= GlobalConstants.Zero
    return % still moving: this was the conservation guard, leave it short
end
% TRANGE, not just its endpoint: a caller that asked for output instants
% (options.tranpoints) gets the ones it has not reached yet, at the settled
% state, exactly as the entry short circuit hands back the whole vector.
pts = trange(:).';
pts = pts(pts > t_iter(end));
if isempty(pts)
    return
end
t_iter = [t_iter(:); pts(:)];
y_iter = [y_iter; repmat(ylast.', numel(pts), 1)];
end
