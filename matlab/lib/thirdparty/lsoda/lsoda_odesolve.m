function [T, Y] = lsoda_odesolve(odefun, tspan, y0, odeopt, mxordn, mxords, forceStiff)
% LSODA_ODESOLVE  LSODA behind MATLAB's ODE-solver signature and odeset struct
%
%   [T, Y] = lsoda_odesolve(odefun, tspan, y0, odeopt, mxordn, mxords, forceStiff)
%
%   The shared body of lsoda_fast, lsoda_accurate, lsoda_fast_stiff and
%   lsoda_accurate_stiff, i.e. of the four hooks in options.odesolvers. It reads
%   the odeset struct the fluid solver builds and honours the fields LSODA has
%   a counterpart for:
%
%     RelTol, AbsTol - tolerances, scalar or one per equation
%     MaxStep        - max step SIZE, mapped to LSODA's hmax. Note that
%                      lsoda_solve reads its own options.MaxStep as a step
%                      COUNT; that is the lower-level contract and is left
%                      alone. Here the odeset meaning applies.
%     InitialStep    - first step size, else LSODA's own heuristic
%     NonNegative    - index vector of components held at or above zero
%     MaxSteps       - LSODA-specific step budget per output interval (10000)
%     OutputFcn      - called as fcn(tspan,y0,'init'), fcn(t,y,'') after every
%                      accepted output point and fcn([],[],'done') at the end;
%                      a nonzero return HALTS the integration and the caller is
%                      handed the trajectory up to that point, as MATLAB's own
%                      solvers do. The fluid solver's conservation guard is one
%                      of these, and it is what stands between a moment-closure
%                      excursion and a window that never returns, so this slot
%                      cannot be a no-op here (see FLUID_CONSERVATION_GUARD).
%
%   NONNEGATIVE is the field LSODA has nothing for, and dropping it is not an
%   option in the fluid solver, where a queue-length coordinate resting at zero
%   dips negative on nearly every step. The rule reproduced here is the one
%   ode15s.m applies and that native Python already mirrors in
%   solver_fld/methods/closing.py:
%
%     - the drift is projected on the floor, so a coordinate at zero is never
%       pushed below it;
%     - an excursion within the band max(0,-ynew) <= AbsTol*RelTol is roundoff
%       on a coordinate resting at zero: the recorded state is clipped and the
%       integrator is left running, because resetting on every such dip pins
%       the step size and the window never advances;
%     - a deeper excursion is charged as error, as ode15s.m charges
%       max(0,-ynew)/AbsTol against RelTol: the step cap is halved and the
%       integrator is RESTARTED from the clipped state, which is what keeps the
%       multistep history from straddling the kink at x = 0. Clipping without
%       the restart leaves the history inconsistent, which is why projection
%       alone does not fix a stiff layer;
%     - the restart carries the step just accepted as its initial step, since a
%       cold restart would re-enter the initial-step heuristic on a settled
%       trajectory and re-select a step orders of magnitude smaller.
%
%   See also: lsoda_matlab, lsoda_solve, lsoda_fast, lsoda_accurate_stiff

if nargin < 4, odeopt = struct(); end
if nargin < 5 || isempty(mxordn), mxordn = 0; end
if nargin < 6 || isempty(mxords), mxords = 0; end
if nargin < 7 || isempty(forceStiff), forceStiff = false; end
if isempty(odeopt), odeopt = struct(); end

rtol = odeField(odeopt, 'RelTol', 1e-6);
atol = odeField(odeopt, 'AbsTol', 1e-9);
hmax = odeField(odeopt, 'MaxStep', 0);
h0 = odeField(odeopt, 'InitialStep', 0);
mxstep = odeField(odeopt, 'MaxSteps', 10000);
nonneg = odeField(odeopt, 'NonNegative', []);
outfcn = odeField(odeopt, 'OutputFcn', []);
mxordn = odeField(odeopt, 'MaxOrdNonStiff', mxordn);
mxords = odeField(odeopt, 'MaxOrdStiff', mxords);
if ~isfinite(hmax), hmax = 0; end
if ~isfinite(h0), h0 = 0; end

y0 = y0(:).';
neq = numel(y0);
tspan = tspan(:).';
nn = false(1, neq);
if ~isempty(nonneg)
    idx = nonneg(:).';
    nn(idx(idx >= 1 & idx <= neq)) = true;
end

opts = struct('hmax', hmax, 'h0', h0, 'forceStiff', forceStiff);

if ~any(nn)
    if forceStiff || hmax > 0 || h0 ~= 0
        [T, Y] = lsoda_matlab(odefun, tspan, y0, rtol, atol, mxstep, mxordn, mxords, opts);
    else
        solveopt = struct('RelTol', rtol, 'AbsTol', atol, 'MaxStep', mxstep, ...
            'MaxOrdNonStiff', mxordn, 'MaxOrdStiff', mxords);
        [T, Y] = lsoda_solve(odefun, tspan, y0, solveopt);
    end
    % Without the NonNegative loop there is no per-step hook to call the output
    % function from, so it is applied to the trajectory that came back: a halt
    % still truncates what the caller sees, which is the contract the caller
    % reads (a short final time), it just costs the whole window first.
    [T, Y] = applyOutputFcn(outfcn, tspan, T, Y);
    return
end

% The band an excursion may sit in without being charged as error. AbsTol is
% the solver tolerance the caller asked for, as in ode15s.m's errNN test.
tolerated = max(rtol) * max(max(atol), max(rtol));
fproj = @(t, yv) projectFloor(odefun(t, yv), yv, nn);

outstate = outputInit(outfcn, tspan, y0);
[T, Y, outstate] = integrateNonneg(fproj, tspan, y0, nn, rtol, atol, ...
    mxstep, mxordn, mxords, hmax, h0, forceStiff, tolerated, outstate);
outputDone(outstate);
end

% ---------------------------------------------------------------------------

function st = outputInit(outfcn, tspan, y0)
% Open an OutputFcn session, in MATLAB's own protocol: 'init' first, then one
% call per accepted output point, then 'done'. An absent function leaves an
% inert state that every hook below short-circuits on.
st = struct('fcn', [], 'halted', false);
if isempty(outfcn)
    return
end
st.fcn = outfcn;
if feval(outfcn, [tspan(1) tspan(end)], y0(:), 'init') ~= 0
    st.halted = true;
end
end

function st = outputStep(st, t, y)
% Report one accepted point. Returns with HALTED set when the function asks to
% stop, which every loop below reads as "return what you have".
if isempty(st.fcn) || st.halted
    return
end
if feval(st.fcn, t, y(:), '') ~= 0
    st.halted = true;
end
end

function outputDone(st)
if ~isempty(st.fcn)
    feval(st.fcn, [], [], 'done');
end
end

function [T, Y] = applyOutputFcn(outfcn, tspan, T, Y)
% The after-the-fact form used on the path that has no per-step hook: replay the
% trajectory through the function and truncate at the first halt.
if isempty(outfcn)
    return
end
st = outputInit(outfcn, tspan, Y(1, :));
last = numel(T);
if ~st.halted
    for i = 1:numel(T)
        st = outputStep(st, T(i), Y(i, :));
        if st.halted
            last = i;
            break
        end
    end
else
    last = 1;
end
outputDone(st);
T = T(1:last);
Y = Y(1:last, :);
end

% ---------------------------------------------------------------------------

function [T, Y, outstate] = integrateNonneg(f, tspan, y0, nn, rtol, atol, mxstep, ...
    mxordn, mxords, hmaxIn, h0In, forceStiff, tolerated, outstate)
% Integrate TSPAN under odeset('NonNegative') semantics, restarting the
% integrator wherever a step lands below the floor by more than roundoff.
%
% A VECTOR TSPAN IS STILL ONE INTEGRATION. lsoda_matlab runs itask = 1 across a
% vector of output times, carrying its Nordsieck history from one instant to the
% next, so an extra output instant costs an interpolation and nothing more.
% Integrating each interval as its own initial-value problem instead COLD-STARTS
% the multistep method at every requested instant: the order falls back to one,
% the initial-step heuristic runs again, and the per-interval local errors
% accumulate along the grid rather than being controlled by a single continuous
% integration. SOLVER_FLUID_PASSAGE_TIME refines its response-time CDF onto up to
% 20001 instants and re-integrates the whole curve on them, which is where that
% showed: the delay row of cdf_respt_closed read an SCV of 1.0624 against the
% exact 1.0 of its Exp service, where one integration reads 1.0070, and the
% three-class model spent twenty minutes on the same call.
y0 = y0(:).';
neq = numel(y0);
touts = tspan(:).';
tdir = 1.0;
if touts(end) < touts(1), tdir = -1.0; end
% Two entries are lsoda_matlab's adaptive mode, where every internal step is an
% output. More than two ask for the solution AT those instants and nothing else.
dense = numel(touts) == 2;
capReset = hmaxIn;
if capReset <= 0, capReset = Inf; end
cap = capReset;
h0 = h0In;
tcur = touts(1);
ycur = y0;
if dense
    nalloc = 256;
else
    nalloc = numel(touts);
end
T = zeros(nalloc, 1);
Y = zeros(nalloc, neq);
T(1) = tcur;
Y(1, :) = ycur;
n = 1;
lastRestart = NaN;
stop = false;
while ~stop && (touts(end) - tcur) * tdir > 0.0
    rest = touts((touts - tcur) * tdir > 0.0);
    if isempty(rest)
        break
    end
    % ONE REMAINING INSTANT puts lsoda_matlab back in its adaptive mode, where
    % every internal step comes back. Only the last row is then an instant the
    % caller asked for; the others are still scanned for excursions, and dropped.
    call = [tcur, rest];
    wantAll = dense || numel(call) > 2;
    opts = struct('hmax', capHmax(cap), 'h0', h0, 'forceStiff', forceStiff);
    [Tk, Yk] = lsoda_matlab(f, call, ycur, rtol, atol, mxstep, mxordn, mxords, opts);
    restarted = false;
    for i = 2:numel(Tk)
        yn = Yk(i, :);
        neg = nn & (yn < 0.0);
        deep = false;
        if any(neg)
            deep = max(-yn(neg)) > tolerated;
            yn(neg) = 0.0;
        end
        if wantAll || i == numel(Tk)
            n = n + 1;
            if n > numel(T)
                T = [T; zeros(numel(T), 1)];  %#ok<AGROW>
                Y = [Y; zeros(size(Y))];      %#ok<AGROW>
            end
            T(n) = Tk(i);
            Y(n, :) = yn;
            outstate = outputStep(outstate, Tk(i), yn);
            if outstate.halted
                stop = true;
                break
            end
        end
        if ~deep
            cap = capReset;
            continue
        end
        span = cap;
        if ~isfinite(span), span = abs(touts(end) - touts(1)); end
        cap = max(span / 2.0, eps * max(1.0, abs(Tk(i))));
        tcur = Tk(i);
        ycur = yn;
        h0 = tdir * min([abs(Tk(i) - Tk(i-1)), cap, abs(touts(end) - tcur)]);
        restarted = true;
        break
    end
    if stop || ~restarted
        break
    end
    if tcur == lastRestart
        error('lsoda:nonnegative', ['lsoda_odesolve: the NonNegative restart at t=%g made no ' ...
            'progress, the drift pushes a floored component below zero faster than the step ' ...
            'size can be cut'], tcur);
    end
    lastRestart = tcur;
end
T = T(1:n);
Y = Y(1:n, :);
end

function hmax = capHmax(cap)
% LSODA reads hmax = 0 as unbounded, the loop tracks it as Inf.
if isfinite(cap)
    hmax = cap;
else
    hmax = 0.0;
end
end

function yp = projectFloor(yp, yv, nn)
% Hold the drift of a component resting on the floor at or above zero.
at = nn(:) & (yv(:) <= 0.0);
if any(at)
    yp(at) = max(yp(at), 0.0);
end
end

function v = odeField(s, name, defval)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = defval;
end
end
