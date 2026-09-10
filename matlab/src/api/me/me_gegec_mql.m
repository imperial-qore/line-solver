function L = me_gegec_mql(lambda, Ca, mu, Cs, c)
%ME_GEGEC_MQL Mean queue length of a stable infinite-capacity GE/GE/c/FCFS queue
%
% Exact ME solution of Kouvatsos (1994), equation (3.9). Shared by the
% infinite-capacity branch of ME_OQN_BLK; ME_OQN carries a numerically
% identical local copy because it is compiled to a MEX file by MEXIFY.
%
% INPUTS:
%   lambda - Arrival rate
%   Ca     - Squared coefficient of variation of the interarrival times
%   mu     - Service rate of one server
%   Cs     - Squared coefficient of variation of the service times
%   c      - Number of servers (finite, c >= 1)
%
% OUTPUT:
%   L      - Mean number of jobs in the queue
%
% Reference:
%   D.D. Kouvatsos, "Entropy Maximisation and Queueing Network Models",
%   Annals of Operations Research, 48:63-126, 1994, equation (3.9).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

alpha2 = 2 / (Cs + 1);
alpha1 = 1 - alpha2;
beta2 = 2 / (Ca + 1);
beta1 = 1 - beta2;
lambda2 = beta2 * lambda;
mu2 = alpha2 * mu;
g = zeros(c, 1);
for j = 1:(c - 1)
    g(j) = (lambda2 + (j - 1) * mu2 * beta1) * alpha2 / (j * mu2 * (1 - alpha1 * beta1));
end
g(c) = (lambda2 + (c - 1) * mu2 * beta1) * alpha2 / (lambda2 * alpha1 + c * mu2);
x = (lambda2 + c * mu2 * beta1) / (lambda2 * alpha1 + c * mu2);
Gn = cumprod(g);
Z = 1 + sum(Gn(1:c-1)) + Gn(c) / (1 - x);
S1 = 0;
for n = 1:(c - 1)
    S1 = S1 + n * Gn(n);
end
S2 = Gn(c) * (c / (1 - x) + x / (1 - x)^2);
L = (S1 + S2) / Z;
end
