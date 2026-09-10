function [W,rhohat]=qsys_gg1(lambda,mu,ca2,cs2)
% [W,RHOHAT]=QSYS_GG1(LAMBDA,MU,CA2,CS2) analyzes a G/G/1 queue.
%
% Uses exact methods for special cases (M/M/1, M/G/1, G/M/1) and
% Allen-Cunneen approximation for the general case. In the G/M/1 case,
% the interarrival-time distribution is fitted from (LAMBDA,CA2) by a
% two-moment renewal process (H2 with balanced means for CA2>1, mixed
% Erlang for CA2<1) and sigma is the root of sigma = A*(mu*(1-sigma)),
% with A* the interarrival-time LST.
%
% Inputs:
%   LAMBDA - Arrival rate
%   MU     - Service rate
%   CA2    - Squared coefficient of variation of inter-arrival time
%   CS2    - Squared coefficient of variation of service time
%
% Returns:
%   W      - Average time in system (response time)
%   RHOHAT - Modified utilization (so that M/M/1 formulas still hold)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

tol = 1e-8;

if abs(ca2 - 1) < tol && abs(cs2 - 1) < tol
    % M/M/1 case
    [W,rhohat] = qsys_mm1(lambda, mu);
elseif abs(ca2 - 1) < tol
    % M/G/1 case (ca2 = 1)
    [W,rhohat] = qsys_mg1(lambda, mu, sqrt(cs2));
elseif abs(cs2 - 1) < tol
    % G/M/1 case (cs2 = 1)
    sigma = qsys_gm1_sigma(lambda, mu, ca2);
    W = qsys_gm1(sigma, mu);
    rhohat = W * lambda / (1 + W * lambda);
else
    % General G/G/1 case - use Allen-Cunneen approximation
    [W,rhohat] = qsys_gig1_approx_allencunneen(lambda, mu, sqrt(ca2), sqrt(cs2));
end

end

function sigma = qsys_gm1_sigma(lambda, mu, ca2)
% Root in (0,1) of sigma = A*(mu*(1-sigma)) for a two-moment fit of the
% interarrival-time LST A*. (Handle-free so that MATLAB Coder can mexify.)
jj = 0; p = 0; nu = 0; p1 = 0; l1 = 0; l2 = 0;
if ca2 >= 1
    % hyperexponential H2 with balanced means
    p1 = (1 + sqrt((ca2-1)/(ca2+1)))/2;
    l1 = 2*p1*lambda;
    l2 = 2*(1-p1)*lambda;
elseif ca2 >= 1e-6
    % mixed Erlang(j-1,j) with common rate (Tijms, 1994)
    jj = ceil(1/ca2);
    p = (jj*ca2 - sqrt(jj*(1+ca2) - jj^2*ca2))/(1+ca2);
    nu = (jj - p)*lambda;
end
% fixed-point iteration; T(x)=A*(mu*(1-x)) is increasing with the queue
% root as its smallest fixed point, so iterates converge monotonically
sigma = lambda/mu;
for it=1:100000
    s = mu*(1-sigma);
    if ca2 < 1e-6
        % deterministic interarrival times
        signew = exp(-s/lambda);
    elseif ca2 < 1
        signew = p*(nu/(s+nu))^(jj-1) + (1-p)*(nu/(s+nu))^jj;
    else
        signew = p1*l1/(s+l1) + (1-p1)*l2/(s+l2);
    end
    if abs(signew-sigma) < 1e-13
        sigma = signew;
        return
    end
    sigma = signew;
end
end
