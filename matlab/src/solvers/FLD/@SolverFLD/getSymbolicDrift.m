function [rhs, vars, sys] = getSymbolicDrift(self, options)
% [RHS, VARS, SYS] = GETSYMBOLICDRIFT(OPTIONS)
%
% Right-hand side of the mean-field ODE system as expression strings, one per
% state variable, together with the variable names they are written in.
%
% This is the input the computer algebra backend needs to produce a Jacobian
% or an equilibrium (see getJacobian), and it is the same system
% solver_fluid_symodes describes and exportODEs typesets, written out
% variable by variable instead of in matrix form.
%
% ONLY SMOOTH DRIFTS ARE EXPORTED. The default, matrix, closing and statedep
% methods scale rates by min(n_i, S_i), which is not differentiable at
% n_i = S_i, so their Jacobian does not exist there; emitting a one-sided
% derivative would be a silent lie exactly at the regime switch that matters.
% Use the p-norm smoothing (options.config.pstar, method matrix or pnorm) or
% the softmin method, whose drifts are smooth everywhere, and this function
% refuses the others by name.
%
% @param options Solver options (optional, defaults to the solver's own)
% @return rhs Cell array of expression strings, one per state variable
% @return vars Cell array of variable names, x1 ... xn
% @return sys The structural description from solver_fluid_symodes
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2 || isempty(options)
    options = self.getOptions();
end
sn = self.getStruct();
sys = solver_fluid_symodes(sn, options);

n = sys.nstates;
vars = cell(1, n);
for s = 1:n
    vars{s} = sprintf('x%d', s);
end

% see _kb/06-solver-catalog.md for rationale
eps0 = num2char(GlobalConstants.FineTol);

switch sys.form
    case 'W'
        if ~strcmp(sys.smoothing, 'pnorm')
            line_error(mfilename, ['The drift of this method scales rates by min(n_i, S_i), ', ...
                'which is not differentiable at n_i = S_i, so it has no Jacobian there. ', ...
                'Set options.config.pstar to use the p-norm smoothing, or use the ''softmin'' method.']);
        end
        rhs = wform_rhs(sys, vars, eps0);
    case 'J'
        rhs = jform_rhs(sys, vars, eps0);
    otherwise
        line_error(mfilename, sprintf('unsupported ODE form ''%s''', sys.form));
end
end

function rhs = wform_rhs(sys, vars, eps0)
% dx/dt = W' * theta(x) + Alambda, with the p-norm smoothed
%   theta_s = x_s / (1 + (n_i/S_i)^p_i)^(1/p_i),  theta_s = 0 at a Source,
% mirroring pnorm_ode in solver_fluid_matrix.
n = sys.nstates;
theta = cell(1, n);
for s = 1:n
    if sys.isSource(s)
        theta{s} = '0';
        continue
    end
    i = sys.stateStation(s);
    ni = stationSum(sys, i, vars, eps0);
    S = sys.S(i);
    p = sys.pstar(i);
    if S <= 0 || p <= 0
        theta{s} = vars{s};
    else
        theta{s} = sprintf('%s/(1 + (%s/%s)^%s)^(1/%s)', vars{s}, ni, ...
            num2char(S), num2char(p), num2char(p));
    end
end

rhs = cell(1, n);
for s = 1:n
    terms = {};
    for t = 1:n
        w = sys.W(t, s);
        if w == 0 || strcmp(theta{t}, '0')
            continue
        end
        terms{end+1} = sprintf('(%s)*(%s)', num2char(w), theta{t}); %#ok<AGROW>
    end
    if sys.Alambda(s) ~= 0
        terms{end+1} = num2char(sys.Alambda(s)); %#ok<AGROW>
    end
    if isempty(terms)
        rhs{s} = '0';
    else
        rhs{s} = strjoin(terms, ' + ');
    end
end
end

function rhs = jform_rhs(sys, vars, eps0)
% dx/dt = J * r(x), with r_e = coeff(e) * factor_e(x). Only the smooth factor
% types are exportable: 'min' (PS/FCFS under closing and statedep), 'fcfsw'
% (statedep FCFS) and 'dpspw' (piecewise DPS) all carry a min or a branch.
smooth = {'lin', 'ext1', 'dps', 'fcfsws'};
rate = cell(1, sys.nevents);
for e = 1:sys.nevents
    ftype = sys.factorType{e};
    if ~any(strcmp(ftype, smooth))
        line_error(mfilename, sprintf(['Event %d scales its rate by the non-smooth factor ', ...
            '''%s'', which has no derivative where the regime switches, so the system has no ', ...
            'Jacobian. Use the ''softmin'' method, or the p-norm smoothing of the ''matrix'' ', ...
            'method.'], e, ftype));
    end
    v = vars{sys.eventVar(e)};
    fdata = sys.factorData{e};
    switch ftype
        case 'lin'
            factor = v;
        case 'ext1'
            % 1 - sum of the class's phases 2..end at the source
            if isempty(fdata.others)
                factor = '1';
            else
                parts = cell(1, numel(fdata.others));
                for k = 1:numel(fdata.others)
                    parts{k} = vars{fdata.others(k)};
                end
                factor = sprintf('(1 - (%s))', strjoin(parts, ' + '));
            end
        case 'dps'
            % ode_rates_closing seeds the denominator with mean(w) and adds
            % no FineTol, so neither does this.
            ntilde = weightedStationSum(sys, fdata.station, sys.dpsw(fdata.station, :), vars, '0');
            factor = sprintf('%s/(%s + %s)', v, num2char(fdata.c0), ntilde);
        case 'fcfsws'
            i = fdata.station;
            % ode_softmin: ni is the raw station total, wni carries FineTol.
            ni = stationSum(sys, i, vars, '0');
            nhat = phaseWeightedStationSum(sys, i, vars, eps0);
            factor = sprintf('%s*(%s)/(%s)', v, softminExpr(ni, num2char(sys.S(i)), sys.alpha), nhat);
    end
    rate{e} = sprintf('(%s)*(%s)', num2char(sys.coeff(e)), factor);
end

rhs = cell(1, sys.nstates);
for s = 1:sys.nstates
    terms = {};
    for e = 1:sys.nevents
        j = sys.J(s, e);
        if j == 0
            continue
        end
        terms{end+1} = sprintf('(%s)*(%s)', num2char(j), rate{e}); %#ok<AGROW>
    end
    if isempty(terms)
        rhs{s} = '0';
    else
        rhs{s} = strjoin(terms, ' + ');
    end
end
end

function s = softminExpr(x, y, alpha)
% Smooth minimum in its weighted-average form,
%   (x e^{-a x} + y e^{-a y}) / (e^{-a x} + e^{-a y}),
% which is what softmin.m computes; softmin.m rewrites it as
% lo + gap*w/(1+w) only to keep the exponent argument non-positive, an
% overflow guard that is meaningless symbolically and would introduce the
% min/max branch this export exists to avoid.
a = num2char(alpha);
s = sprintf('((%s)*exp(-(%s)*(%s)) + (%s)*exp(-(%s)*(%s)))/(exp(-(%s)*(%s)) + exp(-(%s)*(%s)))', ...
    x, a, x, y, a, y, a, x, a, y);
end

function s = stationSum(sys, i, vars, offset)
% Total fluid mass at station i, plus the offset its consumer uses.
idx = find(sys.stateStation == i);
parts = cell(1, numel(idx));
for k = 1:numel(idx)
    parts{k} = vars{idx(k)};
end
if strcmp(offset, '0')
    s = sprintf('(%s)', strjoin(parts, ' + '));
else
    s = sprintf('(%s + %s)', offset, strjoin(parts, ' + '));
end
end

function s = weightedStationSum(sys, i, w, vars, offset)
% sum_r w_ir * n_ir over the classes of station i (DPS denominator).
parts = {};
for k = 1:sys.nstates
    if sys.stateStation(k) == i
        wt = w(sys.stateClass(k));
        if wt ~= 0
            parts{end+1} = sprintf('(%s)*%s', num2char(wt), vars{k}); %#ok<AGROW>
        end
    end
end
s = joinSum(parts, offset);
end

function s = phaseWeightedStationSum(sys, i, vars, offset)
% sum_u w_u x_u over the states of station i, with w_u the mean phase time
% weights (nhat in ode_softmin).
parts = {};
for k = 1:sys.nstates
    if sys.stateStation(k) == i
        wt = sys.fcfsPhaseW(k);
        if wt ~= 0
            parts{end+1} = sprintf('(%s)*%s', num2char(wt), vars{k}); %#ok<AGROW>
        end
    end
end
s = joinSum(parts, offset);
end

function s = joinSum(parts, offset)
% Sum of PARTS, with OFFSET added only when it is not the literal zero.
if isempty(parts)
    s = sprintf('(%s)', offset);
elseif strcmp(offset, '0')
    s = sprintf('(%s)', strjoin(parts, ' + '));
else
    s = sprintf('(%s + %s)', offset, strjoin(parts, ' + '));
end
end

function s = num2char(v)
% Decimal text the symbolic backend reads as an exact rational.
if v == round(v) && abs(v) < 1e15
    s = sprintf('%d', round(v));
else
    s = sprintf('%.17g', v);
end
end
