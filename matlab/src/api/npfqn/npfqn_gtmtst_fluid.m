function result = npfqn_gtmtst_fluid(lambdaFuns, sFuns, muFuns, patienceCcdfs, P, T, varargin)
% NPFQN_GTMTST_FLUID Time-varying open network of many-server fluid queues.
%
% RESULT = NPFQN_GTMTST_FLUID(LAMBDAFUNS, SFUNS, MUFUNS, PATIENCECCDFS, P, T)
% solves an open network of M queues on [0,T]. Each queue is the Gt/Mt/st+GI
% fluid queue of QSYS_GTMTST_FLUID; the departure flow of queue i is routed to
% queue j with proportion P(i,j), and whatever is left leaves the network.
% LAMBDAFUNS, SFUNS, MUFUNS and PATIENCECCDFS are cell arrays with one entry per
% queue; P is an MxM substochastic matrix or a handle P(t) returning one.
%
% THE NETWORK IS A FIXED POINT. The total arrival rate of queue j is
%
%   lambda_j(t) = lambda_j^0(t) + sum_i sigma_i(t) P_ij(t),   sigma_i = mu_i B_i,
%
% (eqs. 23-24) and sigma_i itself depends on lambda_i, because B_i does. The
% iteration starts from the external rates alone and adds one more traversal of
% the network per round, so the nth iterate is the fluid that has made n
% transitions through the network; the map is a monotone contraction, so the
% rates increase to the fixed point rather than oscillating.
%
% Only the SERVICE COMPLETION flow is routed. Abandoning fluid leaves the
% network, which is what makes the traffic equations linear in sigma.
%
% Options:
%   'dt', DT        - grid step, default T/2000
%   'B0', VEC       - initial fluid in service at each queue
%   'w0', VEC       - initial boundary waiting time at each queue
%   'tol', TOL      - sup-norm tolerance on the arrival-rate iteration, 1e-6
%   'maxIter', K    - cap on the iterations, default 100
%
% Returns a struct with fields times, queues (a cell array of the per-queue
% structs of QSYS_GTMTST_FLUID), arrivalRates (the converged total rates, one
% row per queue), iterations and residual.
%
% Example:
%   fc = @(x) exp(-0.5*x);
%   res = npfqn_gtmtst_fluid({@(t) 110+0*t, @(t) 0*t}, {100,80}, {1,1}, {fc,fc}, ...
%                            [0 1; 0 0], 60);
%
% Reference: Y. Liu, W. Whitt (2014). Algorithms for time-varying networks of
% many-server fluid queues. INFORMS J. on Computing 26(1), 59-73.
%
% See also QSYS_GTMTST_FLUID, QSYS_GGISGI_FLUID.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

options = struct('dt', [], 'b0', [], 'w0', [], 'tol', 1e-6, 'maxiter', 100);
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

m = numel(lambdaFuns);
if numel(sFuns) ~= m || numel(muFuns) ~= m || numel(patienceCcdfs) ~= m
    line_error(mfilename, ['Every queue needs an arrival rate, a staffing, a service rate ' ...
        'and a patience law.']);
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
B0 = options.b0;
if isempty(B0)
    B0 = zeros(1, m);
end
w0 = options.w0;
if isempty(w0)
    w0 = zeros(1, m);
end

ext = zeros(m, n);
for i = 1:m
    ext(i,:) = npfqn_gtmtst_fluid_eval(lambdaFuns{i}, t);
end
Pgrid = zeros(n, m, m);
if isa(P, 'function_handle')
    for k = 1:n
        Pgrid(k,:,:) = P(t(k));
    end
else
    if any(size(P) ~= [m m])
        line_error(mfilename, 'The routing matrix must be m x m.');
    end
    for k = 1:n
        Pgrid(k,:,:) = P;
    end
end
if any(Pgrid(:) < -1e-12) || any(sum(Pgrid, 3) > 1 + 1e-9, 'all')
    line_error(mfilename, 'The routing matrix must be substochastic.');
end

lam = ext;
queues = cell(1, m);
residual = Inf;
iter = 0;
for iter = 1:options.maxiter
    sigma = zeros(m, n);
    for i = 1:m
        lamRow = lam(i,:);
        fun = @(u) interp1(t, lamRow, min(max(u, t(1)), t(end)), 'linear');
        queues{i} = qsys_gtmtst_fluid(fun, sFuns{i}, muFuns{i}, patienceCcdfs{i}, T, ...
            'dt', dt, 'B0', B0(i), 'w0', w0(i));
        sigma(i,:) = queues{i}.sigma;
    end
    newlam = ext;
    for j = 1:m
        for i = 1:m
            newlam(j,:) = newlam(j,:) + sigma(i,:) .* reshape(Pgrid(:,i,j), 1, n);
        end
    end
    residual = max(abs(newlam(:) - lam(:)));
    lam = newlam;
    if residual < options.tol
        break
    end
end

result.times = t;
result.queues = queues;
result.arrivalRates = lam;
result.iterations = iter;
result.residual = residual;
end

function y = npfqn_gtmtst_fluid_eval(f, t)
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
