function [lower, upper] = sim_quest_heuristic_ci(bqe, centre, Ap, Np, nstar, alpha, useAutocorr)
% OA_QUEST_HEURISTIC_CI Fallback interval used when a QUEST stage test fails.
%
% [LOWER, UPPER] = OA_QUEST_HEURISTIC_CI(BQE, CENTRE, AP, NP, NSTAR, ALPHA,
% USEAUTOCORR) builds the heuristic interval OA_FQUEST and OA_FIRQUEST deliver
% when the sample is too small for the four stage tests to pass at the smallest
% admissible batch count. BQE holds the K batched quantile estimators, pooled
% over replications for OA_FIRQUEST, CENTRE is the full-sample empirical
% quantile, AP and NP the two single-component variance-parameter estimators,
% and NSTAR the number of observations they were computed from.
%
% Three intervals are formed and the smallest interval containing all of them is
% returned, which is the article's prescription:
%
%   Two symmetric intervals of half-width
%     h = max(t_{1-alpha/2,K} sqrt(AP/NSTAR), t_{1-alpha/2,K-1} sqrt(NP/NSTAR)),
%   one about CENTRE and one about the average of BQE. Taking the wider of the
%   two variance components is deliberately conservative, since neither can be
%   trusted once a stage test has failed.
%
%   Willink's asymmetric interval, which corrects the batched quantile
%   estimators for skewness through the cube-root transform
%   G(zeta) = ([1+6 gamma(zeta-gamma)]^(1/3)-1)/(2 gamma) with
%   gamma = skewness/(6 sqrt(K)), evaluated at both t-quantiles so the two arms
%   differ.
%
% USEAUTOCORR additionally scales the asymmetric arms by
% max(sqrt((1+phi1)/(1-phi1)), 1), with phi1 the lag-1 autocorrelation of BQE.
% Pass true for OA_FQUEST, where the batch quantiles come from one sample path
% and can stay correlated, and false for OA_FIRQUEST, where they come from
% independent replications and the article drops the correction.
%
% Reference: R. Willink, "A Confidence Interval and Test for the Mean of an
% Asymmetric Distribution", Commun. Statist. Theory Methods 34, 2005;
% A. Lolos et al., Proc. Winter Simulation Conference, 2023, step 10, and
% Proc. Winter Simulation Conference, 2025, equations 8 to 10.
%
% See also OA_FQUEST, OA_FIRQUEST
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

bqe = bqe(:);
K = numel(bqe);
if K < 3
    line_error(mfilename, ...
        'the heuristic interval needs at least 3 batch quantiles, got %d', K);
end

half = max(sim_tinv(1 - alpha / 2, K) * sqrt(Ap / nstar), ...
    sim_tinv(1 - alpha / 2, K - 1) * sqrt(Np / nstar));

bqeBar = mean(bqe);
S2 = sum((bqe - bqeBar).^2) / (K - 1);
S2tilde = sum((bqe - centre).^2) / (K - 1);

lower = min(centre - half, bqeBar - half);
upper = max(centre + half, bqeBar + half);

if S2 <= 0
    return
end

skew = (K / ((K - 1) * (K - 2))) * sum(((bqe - bqeBar) / sqrt(S2)).^3);
gamma = skew / (6 * sqrt(K));

varphi = 1;
if useAutocorr
    phi1 = sum((bqe(1:end - 1) - bqeBar) .* (bqe(2:end) - bqeBar)) / ((K - 1) * S2);
    if abs(phi1) < 1
        varphi = max(sqrt((1 + phi1) / (1 - phi1)), 1);
    end
end

tq = sim_tinv(1 - alpha / 2, K - 1);
scale = varphi * sqrt(S2tilde / K);
G1 = sim_willink_g(tq, gamma) * scale;
G2 = sim_willink_g(-tq, gamma) * scale;

lower = min([lower, centre - G1, centre - G2]);
upper = max([upper, centre - G1, centre - G2]);
end

function z = sim_willink_g(zeta, gamma)
% Willink's skewness-adjustment transform, the identity for tiny skewness.
if abs(gamma) <= 0.001
    z = zeta;
    return
end
arg = 1 + 6 * gamma * (zeta - gamma);
% the cube root is taken on the reals, the argument may turn negative for a
% strongly skewed and small batch sample
if arg >= 0
    root = arg^(1 / 3);
else
    root = -((-arg)^(1 / 3));
end
z = (root - 1) / (2 * gamma);
end
