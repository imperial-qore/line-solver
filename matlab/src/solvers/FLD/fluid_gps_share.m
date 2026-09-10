function [s, ds] = fluid_gps_share(xk, wk, vk)
% [S, DS] = FLUID_GPS_SHARE(XK, WK, VK)
%
% Expected capacity share of a GPS station under a normal marginal.
%
% GPS divides the server by WEIGHT among the classes that are BACKLOGGED, and
% then equally among that class's own jobs (State.afterEventStation, GPS
% branch: `cir = min(nir,1)`, share `= w_r/(w*cir')`). The share therefore
% depends on the backlog INDICATOR vector, not on the populations, and that is
% what makes GPS unreachable for a first-order closure: with continuous
% x_k > 0 every class is always backlogged, the indicator is identically one,
% and the share collapses to the constant w_r/sum_j w_j regardless of load.
% That constant is the heavy-traffic limit and is wrong everywhere else, so for
% GPS the second moment is not a correction, it is the entire mechanism.
%
% The closure is an EXACT enumeration rather than an expansion. The share is
% piecewise CONSTANT over the 2^K backlog patterns, so
%
%   E[S_r] = sum_{A ni r} P(backlog set = A) * w_r / sum_{j in A} w_j
%
% carries no truncation error once the pattern probabilities are given. Those
% come from the marginals, P(N_k >= 1) = Phi((x_k - 1/2)/sigma_k) with the
% continuity correction for an integer population, multiplied as if the
% backlogs were independent. That independence is the one approximation
% here and it is not innocuous: in a closed network the station coordinates are
% NEGATIVELY correlated through population conservation, so the exact treatment
% would need multivariate-normal orthant probabilities (closed form to K = 3,
% numerical beyond).
%
% The empty pattern contributes zero share, so sum_r E[S_r] = 1 - P(all
% classes empty) rather than 1. That is deliberate and is how the idle server
% is represented: GPS is single-server, so the backlog indicator plays the role
% that min(n,c) plays at a PS station, and no separate capacity term is applied
% by the caller.
%
% Unlike the DPS ratio closure the expansion is NOT perturbative in sigma: as
% sigma_k -> 0 the probability tends to a step at x_k = 1 and its derivative
% phi(.)/sigma_k diverges, so the Jacobian stiffens at low variance. VK = 0
% falls back to the hard indicator with zero derivative, which is the correct
% mean-field starting point for the outer iteration.
%
% Parameters:
%   xk - (K x 1) per-class populations at the station
%   wk - (K x 1) per-class GPS weights, normalised internally
%   vk - (K x 1) per-class population variances; 0 selects the hard indicator
%
% Returns:
%   s  - (K x 1) expected capacity shares, summing to 1 - P(station empty)
%   ds - (K x K) Jacobian ds_r/dx_m
%
% See also ODE_RATES_CLOSING_FACTORS, FLUID_DRIFT_JACOBIAN, FLUID_SHARE_CLOSURE.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

xk = xk(:);
wk = wk(:);
vk = vk(:);
K = numel(xk);
wantJac = nargout > 1;
s = zeros(K,1);
if wantJac
    ds = zeros(K,K);
end

% the enumeration is 2^K, so refuse rather than crawl; K here is the number of
% classes AT ONE STATION, which is small in every practical model
if K > 12
    line_error(mfilename,sprintf(['GPS closes its capacity share by enumerating the 2^K backlog patterns ' ...
        'of a station, and this station carries %d classes. Above 12 the enumeration is no longer ' ...
        'tractable; use options.method=''closing'' with a DPS station instead.'], K));
end

sw = sum(wk);
if sw <= 0
    return
end
wk = wk / sw;

p = zeros(K,1);
dp = zeros(K,1);
% P(class k backlogged) = P(N_k >= 1) for an INTEGER population, so the normal
% approximation needs the continuity correction P(N_k > 1/2); thresholding at 1
% instead systematically understates the backlog probability, and at sigma = 0
% it makes the share vanish for every class with x_k < 1, which stalls the
% server completely and is an ABSORBING state for the ODE (measured: pop [1 1]
% returned QLen 1.0000/1.0000, every job stuck at the server, against an exact
% 0.6600/0.5400). The sigma = 0 fallback is the fluid limit x_k > 0, matching
% the mean-field statement that positive fluid mass is backlogged.
for k = 1:K
    if vk(k) > 0
        sd = sqrt(vk(k));
        z = (xk(k) - 0.5)/sd;
        p(k) = 0.5*erfc(-z/sqrt(2));
        dp(k) = exp(-0.5*z^2)/(sqrt(2*pi)*sd);
    else
        p(k) = double(xk(k) > 0);
        dp(k) = 0;
    end
end

dsdp = zeros(K,K);
for mask = 1:(2^K - 1)
    A = bitget(mask, 1:K) == 1;
    W = sum(wk(A));
    if W <= 0
        continue % every backlogged class in this pattern carries zero weight
    end
    shareA = wk(A)/W;
    q = p;
    q(~A) = 1 - p(~A); % per-class factor of P(A)
    s(A) = s(A) + prod(q)*shareA;
    if wantJac
        for m = 1:K
            qm = q;
            qm(m) = 1; % product over j ~= m
            dsdp(A,m) = dsdp(A,m) + (2*A(m)-1)*prod(qm)*shareA;
        end
    end
end

if wantJac
    ds = dsdp .* dp(:)'; % chain rule through p_m = Phi((x_m-1)/sigma_m)
end
end
