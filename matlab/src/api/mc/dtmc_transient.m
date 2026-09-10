function pi_t = dtmc_transient(P, pi0, steps)
% PI_T = DTMC_TRANSIENT(P, PI0, STEPS)
%
% Transient distribution of a DTMC with transition matrix P, i.e. the rows
% PI(k) = PI0*P^k for k = 0,...,STEPS. Twin of the Python api.mc.dtmc_transient.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

n = size(P,1);
if nargin < 2 || isempty(pi0)
    pi0 = ones(1,n)/n;
end
pi0 = reshape(pi0, 1, n);
if nargin < 3 || isempty(steps)
    steps = 1;
end
if steps < 0 || steps ~= round(steps)
    line_error(mfilename,'The number of steps must be a non-negative integer.');
end

pi_t = zeros(steps+1, n);
pik = pi0;
pi_t(1,:) = pik;
for k = 1:steps
    pik = pik * P;
    pi_t(k+1,:) = pik;
end
end
