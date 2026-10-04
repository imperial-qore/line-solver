function [W,rhohat]=qsys_gigk_approx_whitt(lambda,mu,ca,cs,k)
% [W,RHOHAT]=QSYS_GIGK_APPROX_WHITT(LAMBDA,MU,CA,CS,K)
%
% GI/G/k approximation of Whitt (1993), eqs. (2.16)-(2.25):
%   Wq = phi(rho,ca^2,cs^2,k) * ((ca^2+cs^2)/2) * Wq(M/M/k)
% where phi interpolates the Cosmetatos M/D/k (phi1) and D/M/k (phi3)
% correction factors. Exact for M/M/k; reduces to the Cosmetatos M/D/k
% approximation for cs=0. Implements eq. (2.25) as printed, which was
% validated here against the paper's Tables 5-7 (New column).
%
% Inputs:
%   LAMBDA - Arrival rate
%   MU     - Service rate per server
%   CA     - Coefficient of variation of inter-arrival time
%   CS     - Coefficient of variation of service time
%   K      - Number of servers
%
% Returns:
%   W      - Average time in system (response time)
%   RHOHAT - Modified utilization (so that M/M/1 formulas still hold)
%
% Reference: Whitt, W. (1993). Approximations for the GI/G/m queue.
% Production and Operations Management 2(2), 114-161.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

ca2 = ca^2;
cs2 = cs^2;
rho = lambda / (k * mu);

% Exact M/M/k baseline (Erlang-C based)
W_mmk = qsys_mmk(lambda, mu, k);
Wq_mmk = W_mmk - 1/mu;

% Cosmetatos correction, as modified by Whitt (1993), eq. (2.17)
gamma = min(0.24, (1-rho)*(k-1)*(sqrt(4+5*k)-2)/(16*k*rho));
phi1 = 1 + gamma;                       % M/D/k factor, eq. (2.16)
phi2 = 1 - 4*gamma;                     % eq. (2.18)
phi3 = phi2 * exp(-2*(1-rho)/(3*rho));  % D/M/k factor, eq. (2.20)
phi4 = min(1, (phi1+phi3)/2);           % eq. (2.21)

c2 = (ca2+cs2)/2;
if c2 >= 1
    psi = 1;                            % eq. (2.22)
else
    psi = phi4^(2*(1-c2));
end

if abs(ca2-cs2) < 1e-12
    phi = psi;                          % eq. (2.25) reduces to psi
elseif ca2 > cs2
    phi = (4*(ca2-cs2)/(4*ca2-3*cs2))*phi1 + (cs2/(4*ca2-3*cs2))*psi;
else
    phi = ((cs2-ca2)/(2*(ca2+cs2)))*phi3 + ((cs2+3*ca2)/(2*(ca2+cs2)))*psi;
end

Wq = phi * c2 * Wq_mmk;                 % eq. (2.24)
W = Wq + 1/mu;
rhohat = W*lambda/(1+W*lambda); % so that M/M/1 formulas still hold

end
