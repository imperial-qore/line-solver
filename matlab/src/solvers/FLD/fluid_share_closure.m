function [s, ds, cn, dcn] = fluid_share_closure(x, wv, C)
% [S, DS, CN, DCN] = FLUID_SHARE_CLOSURE(X, WV, C)
%
% Second-order closure of the capacity share of a sharing discipline. A
% DPS station gives coordinate j the fraction
%
%   S_j = w_j X_j / sum_m w_m X_m
%
% of its capacity, and PS is the same expression with unit weights. The
% first-order fluid closure evaluates that ratio at the mean, which is not
% E[S_j]: the map is a ratio, so Jensen's inequality biases it towards the
% coordinates carrying the LARGER weight. Writing u_j = w_j X_j and
% v = sum_m u_m, the delta method gives
%
%   E[u_j/v] = mu_j/v - Cov(u_j,v)/v^2 + mu_j*Var(v)/v^3 + O(sigma^3)
%
% evaluated at the means. The correction is exactly capacity-conserving:
% sum_j Cov(u_j,v) = Cov(v,v) = Var(v), so the two correction terms cancel in
% the sum and sum_j S_j = 1 identically, as it must for a work-conserving
% discipline. That identity is the invariant to check on any change here.
%
% The expansion is local and fails when v is small against its own standard
% deviation, where the exact expectation is a Cauchy-like integral with no
% finite mean. There a raw share can come out negative; it is clipped at zero
% and the survivors renormalised, which preserves the conservation identity.
%
% THE EXPANSION CARRIES ITS OWN CONVERGENCE RATIO, and both corrections below
% are admitted only while that ratio is below one. Every second-order term
% here is a term of the series for E[1/v], whose successive terms are in the
% ratio
%
%   ratio = Var(v)/v^2
%
% so the truncation is meaningful below 1 and the terms GROW above it.
% Nothing in the algebra notices: at a near-empty station the corrections are
% simply evaluated far outside the region where they mean anything, and come
% back larger than the quantity they correct. On layer 7 of test_LQN_13
% (arbitraryMultiplicity) the station holds 0.5388 jobs with Var(N) = 2, so
% the ratio is 6.9, and Cov(S_j,N) came back at -1.86 where 0.71 is
% all any joint distribution can produce: S_j lies in [0,1], so its variance
% is at most 1/4 and Cauchy-Schwarz caps the covariance at sqrt(Var(N))/2.
% The drift that follows is not integrable -- one mean solve took 287306 drift
% evaluations and 26 s against 774 and 0.03 s for the first-order one, and the
% moment fixed point runs up to ITER_MAX of those.
%
% LOCAL_EXPANSION_WEIGHT therefore scales both corrections by a factor that is
% exactly one while the ratio is at most 1, falling smoothly to zero by 4 (a standard
% deviation twice the mean, where a non-negative v has no mass left near it),
% and is C^1 at both ends so the drift stays differentiable and the integrator
% keeps its step. Outside the region the closure degrades to the first-order
% share u/v, which is what SOLVER_FLUID_CLOSING answers with. Note the factor
% multiplies the SHARE correction and CN together, so the cancellation that
% makes the joint closure exact at an unsaturated station survives it
% identically: there S_j*N = X_j whatever the factor is.
%
% Parameters:
%   x  - (n x 1) coordinate means of one station block
%   wv - (n x 1) per-coordinate weight, constant within a class
%   C  - (n x n) covariance of the same coordinates; empty or zero selects the
%        first-order (plug-in) share
%
% THE SHARE IS ONLY HALF OF THE RATE. What a station clears is S_j*psi(N), and
% the two factors are correlated through N, so E[S_j*psi] is not E[S_j]*E[psi].
% CN returns the missing Cov(S_j, N), from which ODE_RATES_CLOSING_FACTORS
% builds the joint closure E[S_j*psi] = E[S_j]*E[psi] + psi'(n)*Cov(S_j,N).
% That term is what makes the product exact where min() is the identity: at an
% unsaturated station psi = N and psi' = 1, and S_j*N = X_j identically, so the
% two second-order corrections must cancel and leave the plain mean. Separately
% closed they do not, and a five-server station holding 0.77 jobs reported a
% 12.4% queueing delay that cannot exist. CN sums to zero over j, because the
% shares sum to one at every point and Cov(1,N) = 0, so the joint closure is
% work-conserving for exactly the reason the plain share closure is.
%
% Parameters:
%   x  - (n x 1) coordinate means of one station block
%   wv - (n x 1) per-coordinate weight, constant within a class
%   C  - (n x n) covariance of the same coordinates; empty or zero selects the
%        first-order (plug-in) share
%
% Returns:
%   s   - (n x 1) expected shares, summing to one
%   ds  - (n x n) Jacobian ds_j/dx_m, with C held fixed
%   cn  - (n x 1) Cov(S_j, N) at the means, summing to zero
%   dcn - (n x n) Jacobian dcn_j/dx_m, with C held fixed
%
% See also ODE_RATES_CLOSING_FACTORS, FLUID_DRIFT_JACOBIAN, FLUID_MIN_CLOSURE.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

x = x(:);
wv = wv(:);
n = numel(x);
wantJac = nargout > 1;
wantCov = nargout > 2;
u = wv .* x;
v = sum(u);
if wantCov
    cn = zeros(n,1);
    dcn = zeros(n,n);
end
if v <= 0
    s = zeros(n,1);
    if wantJac, ds = zeros(n,n); end
    return
end

s = u / v;
if wantJac
    ds = diag(wv)/v - (u/v^2) * wv';
end

if nargin < 3 || isempty(C) || ~any(C(:))
    return
end

cuv = wv .* (C * wv);   % Cov(u_j, v)
cvv = wv' * C * wv;     % Var(v)

% how far into the series this point sits, and how much of the second-order
% correction that leaves admissible; tau depends on x through v alone, since C
% is held fixed here exactly as the Jacobians below are
[tau, dtau_deps] = local_expansion_weight(cvv / v^2);
if tau <= 0 && ~wantJac
    return  % first-order share, and s/ds are already exactly that
end
dtau = (dtau_deps * (-2*cvv/v^3)) * wv';   % (1 x n) d tau / d x

if wantCov
    % Cov(S_j, N) at the means, from grad(S_j)'*C*1: S_j = w_j X_j / V, so
    % dS_j/dX_a = w_j*delta_aj/V - u_j*w_a/V^2 and the two pieces contract
    % against C*1 = Cov(X,N) and w'*C*1 = Cov(V,N).
    an = C * ones(n,1);     % Cov(X_j, N)
    cvn = wv' * an;         % Cov(V, N)
    cn0 = wv .* an / v - u * (cvn / v^2);
    cn = tau * cn0;
    if nargout > 3
        dcn = tau * ( -(wv .* an) * (wv'/v^2) - (cvn/v^2) * diag(wv) ...
              + (2*cvn/v^3) * u * wv' ) + cn0 * dtau;
    end
end

scorr = -cuv/v^2 + (u*cvv)/v^3;
s = s + tau * scorr;
if wantJac
    ds = ds + tau * ( (2/v^3) * cuv * wv' + (cvv/v^3) * diag(wv) ...
         - (3*cvv/v^4) * u * wv' ) + scorr * dtau;
end

% A COORDINATE CARRYING NO MASS MUST NOT DECIDE THE CLIP. With u_j = 0 the
% plug-in share is zero and the correction leaves s_j = -Cov(u_j,v)/v^2, a
% quantity of the order of rounding whose SIGN is not meaningful. Letting it
% select the branch below zeroes that coordinate's whole Jacobian row, and a
% zero row is an exact zero eigenvalue: on oqn-14-sparse-fcfs, where each class
% visits one station and so leaves the other coordinates empty, the reduced
% spectrum went from -4.4195 to 0 on a 1e-16 change in the closure variance,
% FLUID_LYAPUNOV then declared the fixed point non-hyperbolic, and SolverFLD
% silently fell back from 'minnormal' to 'matrix' -- 0.5 against 0.7503 on the
% queue length. Clip only a share that is negative BEYOND the numerical zero.
if all(s >= -GlobalConstants.Zero)
    s(s < 0) = 0;
    return
end

% the expansion has left the region where it is valid for at least one
% coordinate; clip and renormalise the survivors so the shares still sum to one
act = s > 0;
if ~any(act)
    s = u / v;
    if wantJac
        ds = diag(wv)/v - (u/v^2) * wv';
    end
    return
end
T = sum(s(act));
snew = zeros(n,1);
snew(act) = s(act)/T;
if wantJac
    dT = sum(ds(act,:), 1);
    dsnew = zeros(n,n);
    dsnew(act,:) = ds(act,:)/T - (s(act)/T^2) * dT;
    ds = dsnew;
end
s = snew;
end

function [tau, dtau] = local_expansion_weight(ratio)
% [TAU, DTAU] = LOCAL_EXPANSION_WEIGHT(RATIO) how much of the second-order
% correction the series ratio RATIO = Var(v)/v^2 admits, and d TAU / d RATIO.
%
% One on [0, 1], zero from 4 up, and the C^1 smoothstep between them. Both
% ends matter. The lower one has to be EXACTLY one on the whole convergent
% region, so every model already inside it is bit-identical; the upper one has
% to be reached with a vanishing derivative, because the drift is integrated
% and a kink in it is what collapses the step size.
%
% The two thresholds are the series, not a tuning: at ratio = 1 successive
% terms stop shrinking, at ratio = 4 the standard deviation of v is twice its mean,
% so a non-negative v has essentially no mass near the point being expanded
% about.
lo = 1; hi = 4;
if ratio <= lo
    tau = 1; dtau = 0;
elseif ratio >= hi
    tau = 0; dtau = 0;
else
    t = (ratio - lo)/(hi - lo);
    tau = 1 - t*t*(3 - 2*t);
    dtau = -6*t*(1 - t)/(hi - lo);
end
end
