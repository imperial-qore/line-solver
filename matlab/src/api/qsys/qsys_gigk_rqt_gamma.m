function [Gamma_s,theta] = qsys_gigk_rqt_gamma(rho,mu,Gamma_a,sigma_s,k,alpha_a,regime)
% [GAMMA_S,THETA]=QSYS_GIGK_RQT_GAMMA(RHO,MU,GAMMA_A,SIGMA_S,K,ALPHA_A,REGIME)
%
% Service variability parameter of the Robust Queueing Theory (RQT) framework,
% obtained from the first two moments by the adaptation of Section 7.1,
%
%   Gamma_s = (2 (theta0 + theta1 sigma_s^2/k + theta2 Gamma_a^2 rho^2 k))^((a-1)/a)
%             - Gamma_a k^((a-1)/a),
%
% where (theta0,theta1,theta2) are regressed so that the worst-case system time
% of Theorem 3 approximates the MEAN system time of the corresponding stochastic
% queue. The arrival side needs no adaptation: Gamma_a = sigma_a for an external
% renewal stream, and the network traffic equations (NPFQN_TRAFFIC_RQT) carry it
% to the internal streams. Since the last term cancels Gamma_a at alpha=2, the
% adaptation acts on the sum Gamma_a + Gamma_s/k^(1/alpha) that Theorem 3 reads.
%
% THE FACTOR 2 IS NOT IN THE PRINTED FORMULA and is restored here. Section 7.1
% states that the functional form is motivated by Kingman's bound, which the
% alpha=2 bound of Theorem 3 reproduces when (Gamma_a+Gamma_s)^2 = 2(sigma_a^2 +
% sigma_s^2); the published (theta0,theta1,theta2) are all near unity, i.e. they
% are corrections to that bound rather than a substitute for its factor 2.
% Dropping the factor puts M/M/1 about 40% BELOW its exact mean system time at
% rho=0.9 and M/M/3 about 45% below, which contradicts the errors of at most
% 9.5% that Tables 2-3 report; restoring it gives +4.7% and +8.8%, inside that
% envelope. See _kb/06-solver-catalog.md.
%
% CAUTION: even so the form is not dimensionally homogeneous, since theta0 is an
% additive constant on a scale of variances, so it is only valid in the time unit
% the regression was run in. This routine therefore evaluates it in units of the
% mean service time, 1/mu = 1, and converts the result back. Do not call it with
% a time scale of your own choosing.
%
% Inputs:
%   RHO     - Traffic intensity lambda/(k*mu)
%   MU      - Service rate of each server, which sets the time unit
%   GAMMA_A - Variability parameter of the arrival uncertainty set
%   SIGMA_S - Standard deviation of the service time
%   K       - Number of servers (default 1)
%   ALPHA_A - Effective arrival tail coefficient in (1,2] (default 2)
%   REGIME  - Adaptation regime, Table 1: 'independent' (default, service
%             distribution unknown), 'normal' or 'pareto'
%
% Returns:
%   GAMMA_S - Variability parameter of the service uncertainty set
%   THETA   - The (theta0,theta1,theta2) triple used
%
% Reference: C. Bandi, D. Bertsimas, N. Youssef (2015). Robust Queueing Theory.
% Operations Research 63(3), 676-700, Section 7.1 and Table 1.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5 || isempty(k), k = 1; end
if nargin < 6 || isempty(alpha_a), alpha_a = 2; end
if nargin < 7 || isempty(regime), regime = 'independent'; end

switch lower(char(regime))
    case {'pareto'}
        theta = [-0.05, 1.09, 1.11];
    case {'normal'}
        theta = [-0.02, 1.03, 1.04];
    case {'independent','default'}
        theta = [-0.06, 1.07, 1.07];
    otherwise
        line_error(mfilename, 'Unknown RQT adaptation regime: %s', char(regime));
end

% evaluate in units of the mean service time, then convert back
ga = Gamma_a * mu;
ss = sigma_s * mu;
e = (alpha_a-1)/alpha_a;
b = 2*(theta(1) + theta(2)*ss^2/k + theta(3)*ga^2*rho^2*k);
gs = max(b,0)^e - ga*k^e;
Gamma_s = gs / mu;
end
