function result = sim_asymvar_mm1(lambda, mu)
% SIM_ASYMVAR_MM1 Asymptotic variance of the M/M/1 number-in-system process.
%
% THE QUANTITY THAT MATTERS FOR PLANNING is not the variance of the process but
% its ASYMPTOTIC VARIANCE, sigma^2 = lim t Var(time-average over [0,t]), which is
% twice the integral of the autocovariance. It is what says how long a run must
% be, because the time average of a positively correlated process converges at
% rate sigma^2/t, not at Var(X)/t.
%
% For M/M/1 with utilization rho,
%   E[N] = rho/(1-rho),  Var(N) = rho/(1-rho)^2,
%   sigma^2 = 2 rho (1+rho)/(mu (1-rho)^4).
%
% The FOURTH power is the whole story: the variance of the process grows like
% (1-rho)^-2, but the run length needed to average it away grows like (1-rho)^-4
% over the squared mean, i.e. like (1-rho)^-2. A queue at rho = 0.9 needs about
% 100 times the run of one at rho = 0.
%
% Returns a struct with fields mean, variance, asymptoticVariance and
% relaxationTime (sigma^2/Var, the correlation time scale).
%
% Reference: W. Whitt (1989). Planning queueing simulations. Management Science
% 35(11), 1341-1366; J. Abate, W. Whitt (1988). The correlation functions of RBM
% and M/M/1. Stochastic Models 4(2), 315-359.
%
% See also SIM_RUNLENGTH, SIM_ASYMVAR_CTMC.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if lambda <= 0 || mu <= 0
    line_error(mfilename, 'The arrival and service rates must be positive.');
end
rho = lambda/mu;
if rho >= 1
    line_error(mfilename, 'The queue must be stable, rho < 1.');
end
result.mean = rho/(1-rho);
result.variance = rho/(1-rho)^2;
result.asymptoticVariance = 2*rho*(1+rho)/(mu*(1-rho)^4);
if result.variance > 0
    result.relaxationTime = result.asymptoticVariance/result.variance;
else
    result.relaxationTime = 0;
end
end
