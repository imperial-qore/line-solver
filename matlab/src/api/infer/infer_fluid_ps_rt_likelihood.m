function [ode_h, q_indices, augPhases, LIKE] = infer_fluid_ps_rt_likelihood(sn, taggedClass, y0_levels, Rsampled)
% INFER_FLUID_PS_RT_LIKELIHOOD Fluid-based response time likelihood.
%
% Builds an augmented model with tagged class K+1 by expanding the
% NetworkStruct arrays (following solver_fluid_passage_time's pattern),
% then uses solver_fluid_odes to obtain the ODE handle.
%
% Usage modes:
%   [ode_h, q_indices, augPhases] = infer_fluid_ps_rt_likelihood(sn, taggedClass)
%       Build augmented model and return ODE handle + indexing info.
%
%   [~, ~, ~, LIKE] = infer_fluid_ps_rt_likelihood(sn, taggedClass, y0_levels, Rsampled)
%       Build augmented model, solve ODE, and return likelihood.
%
% Inputs:
%   sn           - NetworkStruct (from model.getStruct())
%   taggedClass  - class index of the tagged job
%   y0_levels    - M x K matrix of fluid levels per station per class (optional)
%   Rsampled     - observed response time (optional)
%
% The augmented model adds class K+1 (tagged) with the same service rates
% as taggedClass. At the reference (queue) station, the tagged class absorbs
% by switching to taggedClass and routing to the delay station. The likelihood
% is extracted from the ODE derivative at the terminal state.
%
% Reference: Casale, G. et al., "Fluid Analysis of Queueing in
% Processor Sharing Systems"
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
% This code is released under the 3-Clause BSD License.

M = sn.nstations;
K = sn.nclasses;
N = sn.nclosedjobs;
Kc = K + 1;

Lambda = sn.mu;
Pi = sn.phi;
PH = sn.proc;
rt = sn.rt;
S = sn.nservers;

% Find reference station (PS queue) and set inf servers to population
refIdx = 0;
for i = 1:M
    if sn.sched(i) ~= SchedStrategy.INF
        refIdx = i;
    end
    if isinf(S(i))
        S(i) = N;
    end
end

%% Build augmented arrays at struct level (following solver_fluid_passage_time)
new_mu = cell(M, 1);
new_pi = cell(M, 1);
new_proc = PH;
for j = 1:M
    new_mu{j} = cell(1, Kc);
    new_pi{j} = cell(1, Kc);
    for k = 1:K
        new_mu{j}{k} = Lambda{j}{k};
        new_pi{j}{k} = Pi{j}{k};
    end
    % Tagged class: copy from taggedClass
    new_mu{j}{Kc} = Lambda{j}{taggedClass};
    new_pi{j}{Kc} = Pi{j}{taggedClass};
    new_proc{j}{Kc} = PH{j}{taggedClass};
end

% Handle NaN/disabled entries
for j = 1:M
    for c = 1:Kc
        if isscalar(new_mu{j}{c}) && isnan(new_mu{j}{c})
            new_mu{j}{c} = [];
            new_pi{j}{c} = [];
        end
    end
end

%% Expand routing table from M*K to M*Kc
new_rt = zeros(M * Kc, M * Kc);

% Copy original routing among base classes
for l = 1:K
    for m = 1:K
        new_rt(l:Kc:end, m:Kc:end) = rt(l:K:end, m:K:end);
    end
end

% Tagged class routes like taggedClass at all stations
new_rt(Kc:Kc:end, Kc:Kc:end) = rt(taggedClass:K:end, taggedClass:K:end);

% Absorption at refIdx: tagged class switches back to original classes
chainIdx = find(sn.chains(:, taggedClass) == 1, 1);
idxClassesInChain = find(sn.chains(chainIdx, :) == 1);
for l = idxClassesInChain
    for j = 1:M
        new_rt((refIdx-1)*Kc + Kc, (j-1)*Kc + l) = rt((refIdx-1)*K + taggedClass, (j-1)*K + l);
    end
end
% Zero out tagged->tagged transitions at refIdx
for j = 1:M
    new_rt((refIdx-1)*Kc + Kc, (j-1)*Kc + Kc) = 0;
end

%% Build ODE using solver_fluid_odes
ode_opts = struct('method', 'default', 'config', struct('hide_immediate', false));
[ode_h, q_indices] = solver_fluid_odes(sn, N, new_mu', new_pi', new_proc, new_rt, S, sn.sched, sn.schedparam, ode_opts);

% Compute augmented phases
augPhases = zeros(M, Kc);
for i = 1:M
    for c = 1:Kc
        if ~isempty(new_mu{i}{c})
            augPhases(i, c) = length(new_mu{i}{c});
        end
    end
end

%% If y0_levels and Rsampled provided, solve ODE and compute likelihood
LIKE = NaN;
if nargin >= 4
    newFluid = 1;

    % Build initial condition from fluid levels
    totalPhases = sum(augPhases(:));
    y0 = zeros(1, totalPhases);

    % Set fluid levels for original classes
    for i = 1:M
        for k = 1:K
            y0(q_indices(i, k)) = y0_levels(i, k);
        end
    end

    % Move fluid from taggedClass to Tagged at refNode
    y0(q_indices(refIdx, taggedClass)) = y0(q_indices(refIdx, taggedClass)) - newFluid;
    y0(q_indices(refIdx, Kc)) = newFluid;

    % Solve ODE from 0 to Rsampled
    refTagIdx = q_indices(refIdx, Kc);
    opt = odeset('AbsTol', 1e-8, 'RelTol', 1e-5, 'NonNegative', 1:length(y0), ...
        'Events', @(t,y) ode_events(t, y, refTagIdx));
    [t, yt] = ode15s(ode_h, [0 Rsampled], y0, opt);

    % Extract likelihood from terminal state derivative
    if Rsampled <= t(end)
        lastState = yt(end,:);
        lastRates = ode_h(t(end), lastState');
        LIKE = -lastRates(refTagIdx) / newFluid;
    else
        LIKE = 0;
    end
end

end

function [value, isterminal, direction] = ode_events(~, y, idx)
    value = y(idx);
    isterminal = 1;
    direction = 0;
end
