function [W,rhohat,Sworst] = qsys_gig1_rqt(lambda,mu,Gamma_a,Gamma_s,alpha_a,alpha_s)
% [W,RHOHAT,SWORST]=QSYS_GIG1_RQT(LAMBDA,MU,GAMMA_A,GAMMA_S,ALPHA_A,ALPHA_S)
%
% Robust Queueing Theory (RQT) worst-case system time of a G/G/1 FCFS queue,
% the single-server case of QSYS_GIGK_RQT. The closed-form bound is Theorem 2,
%
%   W <= (alpha-1)/alpha^(alpha/(alpha-1)) * lambda^(1/(alpha-1))
%        (Gamma_a+Gamma_s)^(alpha/(alpha-1)) / (1-rho)^(1/(alpha-1)) + 1/lambda,
%
% and SWORST is the exact worst case over the uncertainty sets, eq. (12).
%
% Inputs:
%   LAMBDA  - Arrival rate
%   MU      - Service rate
%   GAMMA_A - Variability parameter of the arrival uncertainty set
%   GAMMA_S - Variability parameter of the service uncertainty set
%   ALPHA_A - Arrival tail coefficient in (1,2] (default 2)
%   ALPHA_S - Service tail coefficient in (1,2] (default 2)
%
% Returns:
%   W       - Closed-form bound on the system time (Theorem 2)
%   RHOHAT  - Modified utilization (so that M/M/1 formulas still hold)
%   SWORST  - Exact worst-case system time over the uncertainty sets, eq. (12)
%
% Reference: C. Bandi, D. Bertsimas, N. Youssef (2015). Robust Queueing Theory.
% Operations Research 63(3), 676-700.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5 || isempty(alpha_a), alpha_a = 2; end
if nargin < 6 || isempty(alpha_s), alpha_s = 2; end

if nargout < 3
    [W,rhohat] = qsys_gigk_rqt(lambda,mu,Gamma_a,Gamma_s,1,alpha_a,alpha_s);
else
    [W,rhohat,Sworst] = qsys_gigk_rqt(lambda,mu,Gamma_a,Gamma_s,1,alpha_a,alpha_s);
end
end
