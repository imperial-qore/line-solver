function demandEst = infer_fmlps(model, node, rt, class, ql, W)
% INFER_FMLPS FMLPS demand estimation using sn struct-level operations.
%
% Estimates service demands at a PS queue using the Fluid Maximum
% Likelihood for Processor Sharing method. Uses sn_set_service for
% fast parameter updates, avoiding model.reset()/getStruct() in the
% optimization loop.
%
% Inputs:
%   model  - LINE Network model with delay rates set and queue rates to estimate
%   node   - PS queue node (Station object)
%   rt     - response time samples (column vector, n x 1)
%   class  - class of each sample (column vector, n x 1)
%   ql     - queue lengths at arrival (n x R matrix, per-class)
%   W      - total population (number of threads/jobs)
%
% Returns:
%   demandEst - 1 x R vector of estimated mean service demands
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
% This code is released under the 3-Clause BSD License.

sn = model.getStruct();
R = sn.nclasses;
V = node.getNumberOfServers();
stIdx = node.stationIndex;

xLB = min(rt) * ones(1, R) / W;
xUB = max(rt) * ones(1, R);

% Initial point estimate
meanQL = mean(sum(ql, 2));
Vtilde = min(meanQL, V);
x0 = zeros(1, R);
for j = 1:R
    if ~isempty(rt(class == j))
        x0(j) = Vtilde * mean(rt(class == j)) / meanQL;
    else
        x0(j) = xLB(j);
    end
end

% Cache station indices (invariant across iterations)
delayIdx = 0;
refIdx = 0;
for ii = 1:sn.nstations
    if sn.sched(ii) == SchedStrategy.INF
        delayIdx = ii;
    else
        refIdx = ii;
    end
end

% Cache delay rates (invariant across iterations)
delayRates = zeros(1, R);
for kk = 1:R
    delayRates(kk) = sn.mu{delayIdx}{kk}(1);
end

%% Optimization options
options = optimset();
options.Display = 'iter';
options.Algorithm = 'interior-point';
options.LargeScale = 'on';
options.MaxIter = 1e10;
options.MaxFunEvals = 1e10;
options.MaxSQPIter = 5000;
options.TolCon = 1e-6;
options.TolX = 1e-10;

[demandEst, ~] = fmincon(@objfun, x0, [], [], [], [], xLB, xUB, [], options);

    function f = objfun(x)
        TOL = 1e-6;

        % Update service rates (cell-of-cells format for solver compatibility)
        for r = 1:R
            sn = sn_set_service_coc(sn, stIdx, r, 1/x(r));
        end

        % Build augmented ODE for each unique tagged class (once per iteration)
        uniqueTC = unique(class);
        ftemp = zeros(size(rt));

        for u = 1:length(uniqueTC)
            tc = uniqueTC(u);
            mask = (class == tc);
            rt_tc = rt(mask);
            ql_tc = ql(mask, :);

            % Build augmented model and ODE handle from sn struct directly
            [ode_h, q_idx, augPhases] = infer_fluid_ps_rt_likelihood(sn, tc);

            % Solve each sample reusing the same ODE handle
            ftemp_tc = zeros(size(rt_tc));
            for rr = 1:length(rt_tc)
                like = solve_fmlps_sample(ode_h, q_idx, augPhases, ...
                    delayRates, delayIdx, refIdx, R, ...
                    ql_tc(rr,:), W, tc, rt_tc(rr));
                ftemp_tc(rr) = log(TOL + like);
            end
            ftemp(mask) = ftemp_tc;
        end
        f = -sum(ftemp);
    end

end


function LIKE = solve_fmlps_sample(ode_h, q_indices, augPhases, ...
    delayRates, delayIdx, refIdx, K, aQueue, W, taggedClass, Rsampled)
% Solve the fluid ODE for a single sample and extract the response time likelihood.
%
% Uses the pre-built ODE handle from the augmented model. Only computes
% the sample-specific initial condition and solves the ODE.

newK = K + 1;
newFluid = 1;

% Compute initial fluid levels
y0_levels = zeros(size(augPhases, 1), K);

% Delay station: distribute remaining fluid proportionally to delay rates
delayJobs = (W - sum(aQueue)) * delayRates / sum(delayRates);
y0_levels(delayIdx, :) = delayJobs;

% Queue station: observed queue lengths per class
y0_levels(refIdx, :) = aQueue;

% Build state vector from fluid levels
totalPhases = sum(augPhases(:));
y0 = zeros(1, totalPhases);
M = size(augPhases, 1);
for i = 1:M
    for k = 1:K
        if augPhases(i, k) > 0
            y0(q_indices(i, k)) = y0_levels(i, k);
        end
    end
end

% Move fluid from taggedClass to Tagged at refNode
y0(q_indices(refIdx, taggedClass)) = y0(q_indices(refIdx, taggedClass)) - newFluid;
y0(q_indices(refIdx, newK)) = newFluid;

% Solve ODE from 0 to Rsampled
refTagIdx = q_indices(refIdx, newK);
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


function [value, isterminal, direction] = ode_events(~, y, idx)
    value = y(idx);
    isterminal = 1;
    direction = 0;
end
