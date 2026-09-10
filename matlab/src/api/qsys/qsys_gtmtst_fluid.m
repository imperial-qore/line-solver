function result = qsys_gtmtst_fluid(lambdaFun, sFun, muFun, patienceCcdf, T, varargin)
% QSYS_GTMTST_FLUID The Gt/Mt/st+GI many-server fluid queue.
%
% RESULT = QSYS_GTMTST_FLUID(LAMBDAFUN, SFUN, MUFUN, PATIENCECCDF, T) solves the
% time-varying many-server fluid queue on [0,T]: arrival rate LAMBDAFUN(t),
% staffing SFUN(t), exponential service at rate MUFUN(t), general patience with
% complementary cdf PATIENCECCDF, unlimited waiting room.
%
% THE MODEL ALTERNATES BETWEEN TWO REGIMES, and the whole algorithm is the
% bookkeeping of that alternation:
%
%   UNDERLOADED  the queue is empty and every arrival enters service at once, so
%                the system is the infinite-server fluid model and B obeys
%                B'(t) = lambda(t) - mu(t)B(t) (eq. 18 in its Mt form). It ends
%                when B reaches s while lambda exceeds the rate
%                Gamma(t) = s'(t) + s(t)mu(t) at which capacity frees up (15).
%   OVERLOADED   every server is busy, B(t) = s(t), fluid enters service at
%                exactly Gamma(t), and the queue is described by its BOUNDARY
%                WAITING TIME w(t), the age of the oldest fluid still waiting.
%                Content of age x is what arrived x ago and has not abandoned,
%                q(t,x) = lambda(t-x)F^c(x), and the boundary moves by the delay
%                differential equation (21)
%
%                    w'(t) = 1 - Gamma(t) / [lambda(t-w(t)) F^c(w(t))].
%
%                It ends when w returns to 0 with lambda no longer above
%                Gamma (14).
%
% WHY w AND NOT Q. The queue content is a functional of w, but not the other way
% round: two systems with the same Q and different age profiles abandon at
% different rates. Tracking the boundary keeps the age profile exact, which is
% what makes a general patience law admissible at all.
%
% Options:
%   'dt', DT          - grid step, default T/2000
%   'B0', B           - fluid in service at time 0, default 0
%   'w0', W           - boundary waiting time at time 0, default 0
%   'sPrime', SP      - s'(t); differentiated from SFUN numerically when absent
%   'patiencePdf', F  - the patience density, for the abandonment rate;
%                       differenced from the ccdf when absent
%   'lambdaPast', L   - the arrival rate before time 0, needed only when the
%                       queue starts non-empty
%
% Returns a struct on the grid with fields times, regime (1 overloaded,
% 0 underloaded), B, Q, X = B+Q, w, v (potential waiting time), sigma (service
% completion rate), alpha (abandonment rate), utilization, arrivalRate, staffing
% and capacityRate (Gamma).
%
% Example:
%   fc = @(x) exp(-0.5*x);
%   res = qsys_gtmtst_fluid(@(t) 100+30*sin(t), 100, 1, fc, 40);
%
% Reference: Y. Liu, W. Whitt (2012). The Gt/GI/st+GI many-server fluid queue.
% Queueing Systems 71, 405-444; Y. Liu, W. Whitt (2014). Algorithms for
% time-varying networks of many-server fluid queues. INFORMS J. on Computing
% 26(1), 59-73.
%
% See also QSYS_GGISGI_FLUID, QSYS_MTGINF, NPFQN_GTMTST_FLUID.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

options = struct('dt', [], 'b0', 0, 'w0', 0, 'sprime', [], 'patiencepdf', [], 'lambdapast', []);
for i = 1:2:numel(varargin)
    if i+1 > numel(varargin)
        line_error(mfilename, sprintf('option %s has no value', char(varargin{i})));
    end
    name = lower(char(varargin{i}));
    if ~isfield(options, name)
        line_error(mfilename, sprintf('unknown option %s', char(varargin{i})));
    end
    options.(name) = varargin{i+1};
end

if T <= 0
    line_error(mfilename, 'The horizon T must be positive.');
end
dt = options.dt;
if isempty(dt)
    dt = T/2000;
end
n = round(T/dt) + 1;
t = linspace(0, T, n);
dt = t(2) - t(1);

lam = qsys_gtmtst_fluid_eval(lambdaFun, t);
s = qsys_gtmtst_fluid_eval(sFun, t);
mu = qsys_gtmtst_fluid_eval(muFun, t);
if any(s <= 0)
    line_error(mfilename, 'The staffing function must be positive.');
end
if any(mu <= 0)
    line_error(mfilename, 'The service rate must be positive.');
end
if isempty(options.sprime)
    sp = gradient(s, dt);
else
    sp = qsys_gtmtst_fluid_eval(options.sprime, t);
end

if isempty(options.lambdapast)
    pastFun = lambdaFun;
else
    pastFun = options.lambdapast;
end
lamOf = @(u) qsys_gtmtst_fluid_at(lambdaFun, pastFun, u);
if isempty(options.patiencepdf)
    h = 1e-6;
    pdfOf = @(x) max(0, (patienceCcdf(max(0,x-h)) - patienceCcdf(x+h))/(2*h));
else
    pdfOf = options.patiencepdf;
end

B = zeros(1,n); Q = zeros(1,n); w = zeros(1,n); alpha = zeros(1,n);
regime = zeros(1,n);
B(1) = options.b0;
w(1) = options.w0;
gamma = sp + s.*mu;                        % Gamma(t), eq. (13)

over = w(1) > 0 || (B(1) >= s(1) - 1e-12 && lam(1) > gamma(1));
regime(1) = double(over);
if over
    B(1) = s(1);
end
Q(1) = qsys_gtmtst_fluid_queue(lamOf, patienceCcdf, t(1), w(1), dt);
alpha(1) = qsys_gtmtst_fluid_queue(lamOf, pdfOf, t(1), w(1), dt);

for i = 1:n-1
    if regime(i) == 0
        % Underloaded: B' = lambda - mu B, by RK4 on the grid step.
        f = @(tt,bb) interp1(t, lam, tt, 'linear', 'extrap') - ...
            interp1(t, mu, tt, 'linear', 'extrap')*bb;
        k1 = f(t(i), B(i));
        k2 = f(t(i)+dt/2, B(i)+dt*k1/2);
        k3 = f(t(i)+dt/2, B(i)+dt*k2/2);
        k4 = f(t(i)+dt, B(i)+dt*k3);
        Bnext = B(i) + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
        wnext = 0;
        if Bnext >= s(i+1) && lam(i+1) > gamma(i+1)
            % The servers just filled and the input outruns the freed capacity:
            % eq. (15), the underloaded interval ends here.
            Bnext = s(i+1);
            regime(i+1) = 1;
        else
            regime(i+1) = 0;
            Bnext = min(Bnext, s(i+1));
        end
    else
        % Overloaded: B = s and the boundary moves by eq. (21).
        g = @(tt,ww) qsys_gtmtst_fluid_wdot(tt, ww, lamOf, patienceCcdf, t, gamma);
        k1 = g(t(i), w(i));
        k2 = g(t(i)+dt/2, max(0, w(i)+dt*k1/2));
        k3 = g(t(i)+dt/2, max(0, w(i)+dt*k2/2));
        k4 = g(t(i)+dt, max(0, w(i)+dt*k3));
        wnext = w(i) + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
        Bnext = s(i+1);
        if wnext <= 0 && lam(i+1) <= gamma(i+1)
            % The queue has drained and the input no longer outruns the freed
            % capacity: eq. (14), the overloaded interval ends here.
            wnext = 0;
            regime(i+1) = 0;
        else
            wnext = max(wnext, 0);
            regime(i+1) = 1;
        end
    end
    B(i+1) = Bnext;
    w(i+1) = wnext;
    if regime(i+1) == 1
        Q(i+1) = qsys_gtmtst_fluid_queue(lamOf, patienceCcdf, t(i+1), wnext, dt);
        alpha(i+1) = qsys_gtmtst_fluid_queue(lamOf, pdfOf, t(i+1), wnext, dt);
    else
        Q(i+1) = 0;
        alpha(i+1) = 0;
    end
end

sigma = mu .* B;                            % service completion rate, eq. (3)
% The potential waiting time of an arrival at t is the u-t at which the boundary
% reaches it, i.e. the solution of u - w(u) = t. That map is non-decreasing, so
% one interpolation inverts it.
entry = t - w;
v = max(0, interp1(entry, t, t, 'linear', 'extrap') - t);

result.times = t;
result.regime = regime;
result.B = B;
result.Q = Q;
result.X = B + Q;
result.w = w;
result.v = v;
result.sigma = sigma;
result.alpha = alpha;
result.utilization = B ./ s;
result.arrivalRate = lam;
result.staffing = s;
result.capacityRate = gamma;
end

function y = qsys_gtmtst_fluid_eval(f, t)
% Evaluate a handle or a constant on the grid, accepting a scalar-only handle.
if ~isa(f, 'function_handle')
    y = f*ones(size(t));
    return
end
y = f(t);
if numel(y) == 1 && numel(t) > 1
    y = y*ones(size(t));
elseif numel(y) ~= numel(t)
    y = arrayfun(f, t);
end
y = reshape(y, size(t));
end

function v = qsys_gtmtst_fluid_at(lamFun, pastFun, u)
% The arrival rate at a possibly negative time.
if u < 0
    v = pastFun(u);
else
    v = lamFun(u);
end
v = v(1);
end

function val = qsys_gtmtst_fluid_queue(lamOf, weightFun, ti, wi, dt)
% int_0^w lambda(t-x) WEIGHT(x) dx by Simpson: the ccdf gives the queue content
% that has not abandoned, the density gives the abandonment rate.
if wi <= 0
    val = 0;
    return
end
m = max(8, ceil(wi/dt) + 1);
if mod(m,2) == 1
    m = m + 1;
end
x = linspace(0, wi, m+1);
vals = zeros(1, m+1);
for j = 1:m+1
    vals(j) = lamOf(ti - x(j)) * weightFun(x(j));
end
wgt = ones(1, m+1);
wgt(2:2:end-1) = 4;
wgt(3:2:end-2) = 2;
val = wi/(3*m) * sum(wgt .* vals);
end

function d = qsys_gtmtst_fluid_wdot(tt, ww, lamOf, ccdf, t, gamma)
% Eq. (21): w' = 1 - Gamma(t)/q~(t,w), q~(t,w) = lambda(t-w)F^c(w).
den = lamOf(tt - ww) * ccdf(ww);
if den <= 0
    % No fluid of that age survives, so the boundary can only advance with the
    % clock.
    d = 1;
else
    d = 1 - interp1(t, gamma, tt, 'linear', 'extrap')/den;
end
end
