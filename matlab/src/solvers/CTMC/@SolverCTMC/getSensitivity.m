function [S, SS, dpi, pi] = getSensitivity(self, param, reward, method)
% [S, SS, DPI, PI] = GETSENSITIVITY(PARAM, REWARD, METHOD)
%
% Parametric sensitivity of a steady-state reward to a scalar model
% parameter, following Trivedi and Bobbio (2017), Sec. 9.7.
%
% PARAM describes the parameter theta and how to set it on the model:
%   param.name    identifier used in reports
%   param.value   nominal value theta
%   param.set     handle (model, value) -> void, applying theta to the model
%   param.step    optional finite-difference step, default value*1e-6
%
% METHOD selects how the derivative is taken:
%
%   'fd' (default)  The generator derivative dQ/dtheta is obtained by central
%       differences on the rate with the state space held fixed. This is exact
%       to O(step^2) and requires no symbolic differentiation of the rate
%       assembly; the state space is unaffected because it depends on the
%       topology and the cutoff, not on rate values. The steady-state
%       sensitivity then follows from one linear solve, see ctmc_sens.
%
%   'symbolic'  The stationary distribution is solved as a rational function
%       of the event rate symbols x1..xE and differentiated exactly with
%       respect to each of them by the computer algebra backend (SAGE.m),
%       then combined by the chain rule
%           d(pi)/d(theta) = sum_e d(pi)/d(x_e) * d(x_e)/d(theta).
%       Only the rate map x_e(theta) is still differenced, and that map is
%       affine in theta in the common cases (a rate set to theta, or scaled by
%       it), where the central difference reproduces it exactly. The whole
%       O(step^2) error of 'fd' comes from differencing through the solve,
%       which this avoids entirely. Requires a symbolic backend, and refuses
%       rather than approximates when perturbing theta reshapes an event's
%       filtration instead of merely scaling it.
%
% REWARD is either a function handle mapping the state space to a reward rate
% vector, or a numeric reward rate vector over the states. If omitted, dpi is
% returned and S is empty.
%
% Note: this returns d(E[r])/dtheta with dr/dtheta = 0, i.e. it assumes the
% reward rates do not themselves depend on theta. Rewards that depend on
% theta need the second term of Eq. (9.83) and are not handled here.
%
% @param param Struct describing the parameter, see above
% @param reward Reward rate vector or handle over the state space (optional)
% @param method 'fd' (default) or 'symbolic' (optional)
% @return S Unscaled sensitivity d(E[r])/dtheta, Eq. (9.79)
% @return SS Scaled sensitivity (theta/E[r]) d(E[r])/dtheta, Eq. (9.80)
% @return dpi Sensitivity of the steady-state distribution (1 x n)
% @return pi Steady-state distribution (1 x n)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~isstruct(param) || ~isfield(param, 'set') || ~isfield(param, 'value')
    line_error(mfilename, 'param must be a struct with fields value and set');
end
if nargin < 3
    reward = [];
end
if nargin < 4 || isempty(method)
    method = 'fd';
end
if ~any(strcmpi(method, {'fd', 'symbolic'}))
    line_error(mfilename, sprintf('unknown method ''%s''; expected ''fd'' or ''symbolic''', method));
end

theta = param.value;
if isfield(param, 'step') && ~isempty(param.step)
    h = param.step;
else
    h = max(abs(theta), 1) * 1e-6;
end

% Nominal generator and state space
[Q, ~, ~] = self.getGenerator();
Q = full(Q);
space = self.getStateSpace();
n = size(Q, 1);

if strcmpi(method, 'symbolic')
    [dpi, pi] = symbolicSensitivity(self, param, theta, h, n);
else
    % Central differences on theta with the state space fixed
    Qp = perturbedGenerator(self, param, theta + h);
    Qm = perturbedGenerator(self, param, theta - h);
    if size(Qp, 1) ~= n || size(Qm, 1) ~= n
        line_error(mfilename, ['Perturbing the parameter changed the state space size, so the ', ...
            'generators cannot be differenced. This happens when the parameter switches a ', ...
            'transition on or off (e.g. a zero rate or an immediate transition).']);
    end
    dQ = (Qp - Qm) / (2 * h);

    pi = ctmc_solve(Q);
    dpi = ctmc_sens(Q, dQ, pi);
end

S = [];
SS = [];
if isempty(reward)
    return;
end

if isa(reward, 'function_handle')
    r = reward(space);
else
    r = reward;
end
r = r(:);
if length(r) ~= n
    line_error(mfilename, 'reward must have one entry per state');
end

% Eq. (9.83) with dr/dtheta = 0
S = dpi * r;
Er = pi * r;
if abs(Er) > GlobalConstants.Zero
    SS = (theta / Er) * S;
else
    SS = NaN;
end
end

function [dpi, pi] = symbolicSensitivity(self, param, theta, h, n)
% [DPI, PI] = SYMBOLICSENSITIVITY(SELF, PARAM, THETA, H, N)
%
% Exact d(pi)/d(x_e) from the computer algebra backend, combined with a
% differenced rate map d(x_e)/d(theta) by the chain rule.
%
% The split matters: the stationary distribution is a rational function of the
% rates of high degree, and differencing through it is where the O(h^2) error
% of the 'fd' method comes from. The rate map, by contrast, is affine in theta
% whenever theta is a rate or scales one, and a central difference is exact on
% an affine map. What is left is exact in those cases and no worse otherwise.

% Symbolic generator and the nominal rate of each event. The symbols scale
% filtrations that were normalized by their own minimum positive rate, so the
% nominal value of x_e is that minimum rate.
infGen = self.getSymbolicGenerator();
[~, F] = self.getGenerator();
nEvents = numel(F);
[rate0, shape0] = eventRates(F);

% see _kb/06-solver-catalog.md (CTMC section, symbolic sensitivity) for rationale
hRate = max(abs(theta), 1) * 1e-3;
[Qp, Fp] = perturbedGenerator(self, param, theta + hRate);
[Qm, Fm] = perturbedGenerator(self, param, theta - hRate);
if size(Qp, 1) ~= n || size(Qm, 1) ~= n
    line_error(mfilename, ['Perturbing the parameter changed the state space size, so the ', ...
        'generators cannot be differenced. This happens when the parameter switches a ', ...
        'transition on or off (e.g. a zero rate or an immediate transition).']);
end
if numel(Fp) ~= nEvents || numel(Fm) ~= nEvents
    line_error(mfilename, 'Perturbing the parameter changed the number of events.');
end
[ratep, shapep] = eventRates(Fp);
[ratem, shapem] = eventRates(Fm);
for e = 1:nEvents
    if isempty(shape0{e})
        continue
    end
    if isempty(shapep{e}) || isempty(shapem{e}) || ...
            ~isequal(size(shapep{e}), size(shape0{e})) || ...
            max(max(abs(shapep{e} - shape0{e}))) > 1e-8 || ...
            max(max(abs(shapem{e} - shape0{e}))) > 1e-8
        line_error(mfilename, ['Perturbing the parameter reshapes the filtration of event ', ...
            num2str(e), ' rather than scaling it, so the generator is not linear in a single ', ...
            'rate per event and the symbolic chain rule does not apply. Use the ''fd'' method ', ...
            'for this parameter.']);
    end
end
% see _kb/06-solver-catalog.md (CTMC section, symbolic sensitivity) for rationale
curvature = abs(ratep + ratem - 2 * rate0);
scale = max(1, max(abs(rate0)));
if max(curvature) > 1e-9 * scale
    [~, Fp] = perturbedGenerator(self, param, theta + h);
    [~, Fm] = perturbedGenerator(self, param, theta - h);
    ratep = eventRates(Fp);
    ratem = eventRates(Fm);
    drate = (ratep - ratem) / (2 * h);
else
    drate = (ratep - ratem) / (2 * hRate);
end

% Symbolic stationary distribution, then one exact derivative per symbol.
url = SAGE.resolve(SolverCTMC.symbolicBackend(self));
if isempty(url)
    url = SAGE.require();
end
symbols = cell(1, nEvents);
for e = 1:nEvents
    if ~isempty(shape0{e})
        symbols{e} = ['x', num2str(e)];
    end
end
active = find(~cellfun(@isempty, symbols));
piExpr = SAGE.solveCTMC(infGen, symbols(active), url);
piExpr = SAGE.toExpressionList(piExpr);

assignment = struct();
for e = active
    assignment.(symbols{e}) = rate0(e);
end
pi = SAGE.eval(piExpr, assignment, url);
pi = reshape(pi, 1, []);

dpi = zeros(1, n);
for e = active
    if drate(e) == 0
        % This event does not depend on theta, so its term is zero and the
        % derivative is not worth a round trip.
        continue
    end
    dExpr = SAGE.diff(piExpr, symbols{e}, 1, url);
    dvals = SAGE.eval(dExpr, assignment, url);
    dpi = dpi + drate(e) * reshape(dvals, 1, []);
end
end

function [rates, shapes] = eventRates(F)
% [RATES, SHAPES] = EVENTRATES(F)
% Minimum positive rate of each event filtration, and the filtration
% normalized by it. An event with no positive rate contributes neither.
rates = zeros(1, numel(F));
shapes = cell(1, numel(F));
for e = 1:numel(F)
    Fe = full(F{e});
    pos = Fe(Fe > 0);
    if isempty(pos)
        continue
    end
    rates(e) = min(pos);
    shapes{e} = Fe / rates(e);
end
end

function [Q, F] = perturbedGenerator(self, param, value)
% [Q, F] = PERTURBEDGENERATOR(SELF, PARAM, VALUE)
% Rebuild the generator with theta set to VALUE, on a copy of the model so
% the caller's model is left untouched.
%
% The hard refresh is required, not defensive: setService and setArrival
% deliberately leave the cached struct in place, so a copy that inherited a
% built struct would report the old rate and the difference quotient would
% silently come out as zero.

modelCopy = self.model.copy();
param.set(modelCopy, value);
modelCopy.refreshStruct(true);
solverCopy = SolverCTMC(modelCopy, self.getOptions());
[Q, F] = solverCopy.getGenerator();
Q = full(Q);
end
