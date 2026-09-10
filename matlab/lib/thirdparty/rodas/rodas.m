function [T, Y, TE, YE, IE, stats] = rodas(odefun, tspan, y0, odeopt)
% RODAS  Hairer-Wanner RODAS behind MATLAB's ODE-solver signature
%
%   [T, Y] = RODAS(ODEFUN, TSPAN, Y0, ODEOPT)
%   [T, Y, TE, YE, IE] = RODAS(ODEFUN, TSPAN, Y0, ODEOPT)  with EVENTS set
%
%   An index-1 DAE integrator for  M y' = f(t,y)  with a possibly SINGULAR M,
%   presented with the calling convention and the ODESET struct of @ode15s so
%   that it drops into OPTIONS.ODESOLVERS. It is the counterpart of
%   LSODA_ODESOLVE, which fills the same role for the plain ODE slots.
%
%   WHERE THIS IS SELECTED. SolverFLD's 'dae' method integrates its transient
%   through OPTIONS.ODESOLVERS.DAESOLVER, which defaults to @ode15s; setting it
%   to @rodas runs the step sequence that the C++, native Python and JAR ports
%   of the DAE route also run, so all four agree in the last digits. It is a
%   SEPARATE slot from ACCURATESTIFFODESOLVER on purpose -- see NONNEGATIVE
%   below.
%
%   ODESET fields honoured:
%     Mass          the mass matrix M. A constant matrix, or a function handle
%                   evaluated ONCE at TSPAN(1): rodas.f calls MAS a single time
%                   before the first step and has no channel for a matrix that
%                   varies, so a state- or time-dependent mass is refused by
%                   name rather than silently frozen. Banded M is detected and
%                   passed in band storage, which is the IJOB=3 path the fluid
%                   DAE uses; otherwise M goes in full (IJOB=5). Absent M means
%                   the identity (IJOB=1).
%     MassSingular  informational here: RODAS carries a singular M natively and
%                   needs no declaration. Accepted for ODE15S compatibility.
%     Jacobian      df/dy, a constant matrix or @(t,y). Absent, the Jacobian is
%                   formed by finite differences, as RODAS's own IJAC=0.
%     RelTol        scalar (default 1e-3)
%     AbsTol        scalar or one per equation (default 1e-6)
%     InitialStep   first step size (default 1e-6, RODAS's own)
%     MaxStep       largest step size (default TSPAN(end)-TSPAN(1))
%     Events        @(t,y) -> [value, isterminal, direction], ode15s's own
%                   convention. A sign change of VALUE across an accepted step is
%                   located ON THE STEP'S OWN INTERPOLANT by bisection -- RODAS
%                   carries a third-order one (CONTRO) for exactly this kind of
%                   question -- so locating a crossing costs no extra step and the
%                   state reported at it is the integrator's, not a re-solve's.
%                   The integration STOPS at a terminal event and the reported
%                   path ends there: no point beyond the crossing is returned,
%                   which is what makes the located time worth locating for the
%                   hybrid DAE that asks for it (SOLVER_FLUID_DAE, where the
%                   crossing is a capacity constraint starting or stopping to
%                   bind). A non-terminal event is recorded and the integration
%                   continues to the next accepted step, so it is located to the
%                   same tolerance but not stepped back to.
%
%   THE STATS ARE THE SIXTH OUTPUT, not the third: TE/YE/IE take positions 3 to 5
%   so that the signature is ode15s's, which is the point of this wrapper.
%
%   NONNEGATIVE IS REFUSED, NOT IGNORED. Holding a component at or above zero
%   is a decision taken INSIDE the step loop -- ode15s charges the excursion
%   against RelTol and reduces the step -- and rodas_core is a faithful
%   transliteration of vendored Fortran that has no such hook. Honouring the
%   field would mean editing the vendored numerics; ignoring it would return
%   negative queue lengths from a solver that was told not to. So it errors,
%   which is also why this is not wired into ACCURATESTIFFODESOLVER: the rest
%   of the fluid solver does set NonNegative, and the DAE route does not.
%
%   TSPAN with two entries reports at every accepted step, as ode15s does. With
%   more than two it reports at exactly those abscissae, read off CONTRO, the
%   third-order interpolant RODAS carries per step -- so an output grid costs no
%   extra step and no interpolation of the caller's own.
%
%   See also RODAS_CORE, RODAS_CONTRO, ODE15S, LSODA_ODESOLVE.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 4 || isempty(odeopt), odeopt = struct(); end
y0 = y0(:);
n = numel(y0);
tspan = tspan(:);
if numel(tspan) < 2
    error('rodas:tspan', 'TSPAN needs at least an initial and a final value.');
end
t0 = tspan(1);
tf = tspan(end);
if tf == t0
    error('rodas:tspan', 'TSPAN must span a nonzero interval.');
end

if local_has(odeopt, 'NonNegative') && ~isempty(odeopt.NonNegative)
    error('rodas:nonNegative', ...
        ['RODAS cannot hold a component at or above zero: that is decided ' ...
         'inside the step loop and rodas_core is a faithful transliteration ' ...
         'with no such hook. Use @ode15s for a path that sets NonNegative, ' ...
         'or drop the constraint.']);
end
if local_has(odeopt, 'MStateDependence') && ...
        ~isempty(odeopt.MStateDependence) && ...
        ~strcmpi(odeopt.MStateDependence, 'none')
    error('rodas:massState', ...
        ['RODAS reads the mass matrix ONCE, before the first step, so a ' ...
         'state-dependent mass has nowhere to be re-read. Set ' ...
         'MStateDependence to ''none'' or use @ode15s.']);
end

rtol = 1e-3;
if local_has(odeopt, 'RelTol') && ~isempty(odeopt.RelTol)
    rtol = odeopt.RelTol(1);
end
atol = 1e-6;
if local_has(odeopt, 'AbsTol') && ~isempty(odeopt.AbsTol)
    atol = odeopt.AbsTol;
end
if isscalar(atol)
    itol = 0;
    atolv = atol;
    rtolv = rtol;
else
    itol = 1;
    atolv = atol(:);
    if numel(atolv) ~= n
        error('rodas:absTol', ...
            'AbsTol has %d entries where the system has %d.', numel(atolv), n);
    end
    rtolv = repmat(rtol, n, 1);
end

opt = struct();
opt.ifcn = 1;    % f may depend on t; the extra df/dt difference is cheap and
                 % an autonomous system simply returns zero for it
opt.idfx = 0;
opt.uround = 1e-16;
opt.nmax = 100000;
opt.meth = 1;
opt.pred = true;
opt.safe = 0.9;
opt.fac1 = 5;
opt.fac2 = .16666666666666666;
opt.m1 = 0;
opt.m2 = 0;
opt.mljac = n;
opt.mujac = n;

h0 = 1e-6;
if local_has(odeopt, 'InitialStep') && ~isempty(odeopt.InitialStep)
    h0 = odeopt.InitialStep;
end
if local_has(odeopt, 'MaxStep') && ~isempty(odeopt.MaxStep)
    opt.hmax = odeopt.MaxStep;
else
    opt.hmax = tf - t0;
end

% ---- the Jacobian ----
if local_has(odeopt, 'Jacobian') && ~isempty(odeopt.Jacobian)
    Jspec = odeopt.Jacobian;
    if isa(Jspec, 'function_handle')
        opt.jac = @(x, y) full(Jspec(x, y));
    else
        Jc = full(Jspec);
        opt.jac = @(x, y) Jc;
    end
    opt.ijac = 1;
else
    opt.ijac = 0;
    opt.jac = [];
end

% ---- the mass matrix ----
if local_has(odeopt, 'Mass') && ~isempty(odeopt.Mass)
    Mspec = odeopt.Mass;
    if isa(Mspec, 'function_handle')
        % rodas.f calls MAS once, so this is the value it would have used; a
        % genuinely varying mass was refused above.
        try
            Mm = Mspec(t0, y0);
        catch
            Mm = Mspec(t0);
        end
    else
        Mm = Mspec;
    end
    Mm = full(Mm);
    if ~isequal(size(Mm), [n n])
        error('rodas:mass', 'Mass is %dx%d where the system has %d equations.', ...
            size(Mm,1), size(Mm,2), n);
    end
    [mlmas, mumas] = local_bandwidth(Mm);
    % Band storage pays only while the band is narrower than the matrix; at or
    % above that it IS the full case, which is also how rodas.f decides (it
    % compares MLMAS with NM1).
    if mlmas + mumas + 1 >= n
        opt.mlmas = n;
        opt.mumas = n;
        opt.mas = @() Mm;
    else
        opt.mlmas = mlmas;
        opt.mumas = mumas;
        band = zeros(mlmas + mumas + 1, n);
        for j = 1:n
            for i = max(1, j-mumas):min(n, j+mlmas)
                band(i - j + mumas + 1, j) = Mm(i, j);
            end
        end
        opt.mas = @() band;
    end
    opt.imas = 1;
else
    opt.imas = 0;
    opt.mas = [];
    opt.mlmas = 0;
    opt.mumas = 0;
end

fcn = @(x, y) local_rhs(odefun, x, y);

% ---- events ----
evfun = [];
if local_has(odeopt, 'Events') && ~isempty(odeopt.Events)
    evfun = odeopt.Events;
end
TE = zeros(0,1);
YE = zeros(0,n);
IE = zeros(0,1);
gprev = zeros(0,1);
stoppedByEvent = false;
if ~isempty(evfun)
    gprev = local_evvalue(evfun, t0, y0);
end

% ---- output plan ----
grid = tspan;
wantGrid = numel(tspan) > 2;
if wantGrid
    d = diff(grid);
    if ~(all(d > 0) || all(d < 0))
        error('rodas:tspan', 'TSPAN must be monotonic.');
    end
end

Tacc = zeros(0,1);
Yacc = zeros(0,n);
cursor = 1;            % next grid point still owed a value
% t0 is the initial condition, which no integrator has to be asked for.
if wantGrid
    Tacc(end+1,1) = grid(1);
    Yacc(end+1,:) = y0.';
    cursor = 2;
else
    Tacc(end+1,1) = t0;
    Yacc(end+1,:) = y0.';
end

    function irtrn = solout(nr, xold, x, y, dense) %#ok<INUSD>
        irtrn = 0;
        if ~isempty(evfun) && nr > 1
            gcur = local_evvalue(evfun, x, y);
            [hitIdx, hitTerm] = local_evcross(evfun, gprev, gcur, x, y);
            if ~isempty(hitIdx)
                % BISECT ON THE INTERPOLANT rather than on the integration: the
                % step that has just been accepted carries its own third-order
                % dense output, so the crossing is located inside it.
                lo = xold; hi = x;
                for bit = 1:60
                    mid = 0.5*(lo + hi);
                    ymid = zeros(n,1);
                    for ii = 1:n
                        ymid(ii) = rodas_contro(dense, ii, mid);
                    end
                    gmid = local_evvalue(evfun, mid, ymid);
                    if local_evsame(gprev(hitIdx(1)), gmid(hitIdx(1)))
                        lo = mid;
                    else
                        hi = mid;
                    end
                    if abs(hi - lo) <= 1e-12*max(1, abs(hi))
                        break
                    end
                end
                yhit = zeros(n,1);
                for ii = 1:n
                    yhit(ii) = rodas_contro(dense, ii, hi);
                end
                TE(end+1,1) = hi; %#ok<AGROW>
                YE(end+1,:) = yhit.'; %#ok<AGROW>
                IE(end+1,1) = hitIdx(1); %#ok<AGROW>
                if hitTerm
                    % THE OVERSHOOTING STEP IS NOT REPORTED: the path ends at the
                    % crossing, which is where the caller restarts it.
                    keep = Tacc <= hi + 1e-13*max(1, abs(hi));
                    Tacc = Tacc(keep);
                    Yacc = Yacc(keep,:);
                    Tacc(end+1,1) = hi;
                    Yacc(end+1,:) = yhit.';
                    stoppedByEvent = true;
                    irtrn = -1;
                    return
                end
            end
            gprev = gcur;
        end
        if wantGrid
            while cursor <= numel(grid) && ...
                    ((tf > t0 && grid(cursor) <= x + 1e-13) || ...
                     (tf < t0 && grid(cursor) >= x - 1e-13))
                tq = grid(cursor);
                if nr <= 1
                    % Before any step: CONT holds no coefficients yet and the
                    % state at that point is y itself.
                    xs = y(:).';
                else
                    xs = zeros(1, n);
                    for ii = 1:n
                        xs(ii) = rodas_contro(dense, ii, tq);
                    end
                end
                Tacc(end+1,1) = tq; %#ok<AGROW>
                Yacc(end+1,:) = xs; %#ok<AGROW>
                cursor = cursor + 1;
            end
        else
            if nr > 1
                Tacc(end+1,1) = x; %#ok<AGROW>
                Yacc(end+1,:) = y(:).'; %#ok<AGROW>
            end
        end
    end

opt.solout = @solout;
opt.iout = 1;

[yend, idid, stats, xend_reached] = rodas_core(n, fcn, t0, y0, tf, h0, ...
    rtolv, atolv, itol, opt);

if idid ~= 1 && ~(stoppedByEvent && idid == 2)
    error('rodas:failed', ...
        'RODAS returned idid=%d at t=%g (%s).', idid, xend_reached, ...
        local_idid(idid));
end

% RODAS lands exactly on TF, and the last accepted step's dense output is
% evaluated at its own right endpoint, where the two agree to rounding; take
% the integrator's own value there rather than the interpolant's. A run stopped
% by a terminal event already ends AT the crossing and must not be extended to
% the step that overshot it.
if ~stoppedByEvent
    if ~isempty(Tacc) && abs(Tacc(end) - tf) <= 1e-13*max(1, abs(tf))
        Yacc(end,:) = yend.';
    else
        Tacc(end+1,1) = xend_reached;
        Yacc(end+1,:) = yend.';
    end
end

T = Tacc;
Y = Yacc;
end

% ---------------------------------------------------------------------------
function g = local_evvalue(evfun, t, y)
g = evfun(t, y(:));
g = g(:);
end

% ---------------------------------------------------------------------------
function [idx, terminal] = local_evcross(evfun, gprev, gcur, t, y)
% Which event functions crossed zero across this step, in the direction they
% asked for. ODE15S's convention: DIRECTION -1 counts only a decreasing crossing,
% +1 only an increasing one, 0 either.
idx = zeros(0,1);
terminal = false;
[~, isterm, dirn] = evfun(t, y(:));
isterm = isterm(:);
dirn = dirn(:);
for k = 1:min(numel(gprev), numel(gcur))
    if gprev(k) == 0 || sign(gcur(k)) == sign(gprev(k))
        continue
    end
    dk = 0;
    if numel(dirn) >= k
        dk = dirn(k);
    end
    if dk < 0 && ~(gprev(k) > 0 && gcur(k) <= 0)
        continue
    end
    if dk > 0 && ~(gprev(k) < 0 && gcur(k) >= 0)
        continue
    end
    idx(end+1,1) = k; %#ok<AGROW>
    if numel(isterm) >= k && isterm(k)
        terminal = true;
    end
end
end

% ---------------------------------------------------------------------------
function tf = local_evsame(a, b)
% Same side of zero, which is what the bisection halves on.
tf = (a > 0 && b > 0) || (a < 0 && b < 0) || (a == 0 && b == 0);
end

% ---------------------------------------------------------------------------
function f = local_rhs(odefun, x, y)
f = odefun(x, y);
f = f(:);
end

% ---------------------------------------------------------------------------
function tf = local_has(s, f)
tf = isstruct(s) && isfield(s, f);
end

% ---------------------------------------------------------------------------
function [ml, mu] = local_bandwidth(M)
% The tightest band that still holds every nonzero of M.
[i, j] = find(M ~= 0);
if isempty(i)
    ml = 0; mu = 0;
    return
end
ml = max(0, max(i - j));
mu = max(0, max(j - i));
end

% ---------------------------------------------------------------------------
function s = local_idid(idid)
switch idid
    case 2,  s = 'stopped by the output routine';
    case -2, s = 'more steps than nmax';
    case -3, s = 'step size became too small';
    case -4, s = 'the matrix is repeatedly singular';
    otherwise, s = 'unknown status';
end
end
