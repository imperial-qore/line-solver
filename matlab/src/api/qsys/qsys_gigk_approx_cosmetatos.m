function [W,rhohat]=qsys_gigk_approx_cosmetatos(lambda,mu,ca,cs,k)
% [W,RHOHAT]=QSYS_GIGK_APPROX_COSMETATOS(LAMBDA,MU,CA,CS,K)
%
% GI/G/k approximation by interpolation of the M/M/k, M/D/k and D/M/k
% queues (Cosmetatos 1982; Page 1982):
%   Wq = [ca^2*cs^2 + ca^2*(1-cs^2)*phi1/2 + (1-ca^2)*cs^2*phi3/2]*Wq(M/M/k)
% where phi1 and phi3 are the Cosmetatos (1975) correction factors for
% M/D/k and D/M/k, with the safeguards of Whitt (1993). The D/D/k corner
% has Wq=0. The interpolation requires ca^2<=1 and cs^2<=1; outside this
% region the Lee-Longton scaling Wq = ((ca^2+cs^2)/2)*Wq(M/M/k) is used.
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
% References: Cosmetatos, G.P. (1975). Approximate explicit formulae for
% the average queueing time in the processes (M/D/r) and (D/M/r).
% INFOR 13, 328-331. Page, E. (1982). Tables of waiting times for M/M/n,
% M/D/n and D/M/n and their use to give approximate waiting times in more
% general queues. J. Opl. Res. Soc. 33, 453-473.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

ca2 = ca^2;
cs2 = cs^2;
rho = lambda / (k * mu);

% Exact M/M/k baseline waiting time (Erlang-C based)
W_mmk = qsys_mmk(lambda, mu, k);
Wq_mmk = W_mmk - 1/mu;

if ca2 <= 1 && cs2 <= 1
    % Cosmetatos correction, as modified by Whitt (1993), eq. (2.17)
    gamma = min(0.24, (1-rho)*(k-1)*(sqrt(4+5*k)-2)/(16*k*rho));
    phi1 = 1 + gamma;                       % M/D/k factor
    phi3 = (1 - 4*gamma) * exp(-2*(1-rho)/(3*rho)); % D/M/k factor
    Wq = (ca2*cs2 + ca2*(1-cs2)*phi1/2 + (1-ca2)*cs2*phi3/2) * Wq_mmk;
else
    % Interpolation weights are invalid outside the unit box
    Wq = ((ca2+cs2)/2) * Wq_mmk;
end
W = Wq + 1/mu;

rhohat = W*lambda/(1+W*lambda); % so that M/M/1 formulas still hold

end
