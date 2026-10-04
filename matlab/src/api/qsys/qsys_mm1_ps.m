function [W, W2, alpha] = qsys_mm1_ps(lambda, mu)
% [W,W2,ALPHA] = QSYS_MM1_PS(LAMBDA, MU)
%
% Exact sojourn-time moments of the multiclass M/M/1-PS queue.
%
% Class j arrives in a Poisson stream of rate LAMBDA(j) and requires an
% exponential amount of service with rate MU(j). The processor is shared
% equally by all jobs in service, so the class of a job affects its sojourn
% time both through its own service rate and through the mix of rates of the
% jobs it shares the processor with. With ALPHA = 1 - sum_j LAMBDA(j)/MU(j)
% the unutilized fraction of the processor, the moments of the sojourn time
% W_r of a tagged class-r job are
%
%   E[W_r]   = 1/(ALPHA*MU(r))
%   E[W_r^2] = 2/(ALPHA*MU(r))^2 * [1 - sum_j lambda_j (mu_j-mu_r)/(mu_j(mu_j+mu_r))]
%                                / [1 - sum_j lambda_j/(mu_j+mu_r)]
%
% which is equation (7) of Mitra and Morrison (1983). Both are exact, not
% asymptotic: the open system is the N -> infinity limit of the closed
% terminal-driven system whose moments that paper expands in 1/N, and the
% leading term of the expansion is exact in the limit. For a single class the
% second moment reduces to the classical 4/(mu^2 (1-rho)^2 (2-rho)) of Coffman,
% Muntz and Trotter (1970).
%
% Input:
%   lambda - per-class Poisson arrival rates (1,R), non-negative
%   mu     - per-class exponential service rates (1,R), positive
%
% Output:
%   W     - per-class mean sojourn times (1,R)
%   W2    - per-class second moments of the sojourn time (1,R)
%   alpha - unutilized fraction of the processor, 1 - sum_j lambda_j/mu_j
%
% Reference: D. Mitra, J. A. Morrison, "Asymptotic Expansions of Moments of
% the Waiting Time in Closed and Open Processor-Sharing Systems with Multiple
% Job Classes", Adv. Appl. Prob. 15(4), 1983, equation (7).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

lambda = lambda(:)';
mu = mu(:)';
R = numel(lambda);
if numel(mu) ~= R
    line_error(mfilename, 'lambda and mu must have the same number of classes');
end
if any(~isfinite(lambda)) || any(lambda < 0)
    line_error(mfilename, 'lambda must be finite and non-negative');
end
if any(~isfinite(mu)) || any(mu <= 0)
    line_error(mfilename, 'mu must be finite and positive');
end

alpha = 1 - sum(lambda ./ mu);
if alpha <= 0
    line_error(mfilename, sprintf('System is unstable: utilization %.6f >= 1', 1 - alpha));
end

W = zeros(1, R);
W2 = zeros(1, R);
for r = 1:R
    mur = mu(r);
    num = 1 - sum(lambda .* (mu - mur) ./ (mu .* (mu + mur)));
    den = 1 - sum(lambda ./ (mu + mur));
    W(r) = 1 / (alpha * mur);
    W2(r) = 2 / (alpha * mur)^2 * num / den;
end
end
