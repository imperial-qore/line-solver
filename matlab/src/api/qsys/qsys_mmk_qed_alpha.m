function alpha = qsys_mmk_qed_alpha(beta)
% ALPHA = QSYS_MMK_QED_ALPHA(BETA)
%
% The Halfin-Whitt delay-probability function
%
%   alpha(beta) = [ 1 + beta Phi(beta)/phi(beta) ]^(-1),   beta > 0,
%
% with phi and Phi the standard normal density and cdf. It is the limit of the
% Erlang C delay probability of the M/M/s queue as s -> Inf with
% beta = (1-rho)sqrt(s) held fixed, and it decreases strictly from 1 at beta = 0
% to 0 as beta -> Inf, which is what makes it invertible for staffing.
%
% BETA may be an array. Non-positive entries return 1: with no server slack
% every arrival is delayed.
%
% Evaluated as phi/(phi + beta*Phi) rather than as the reciprocal of
% 1 + beta*Phi/phi. The two are the same function, but the quotient Phi/phi
% overflows once phi underflows (beta beyond about 38), whereas this form
% degrades to 0/(0+beta) = 0, which is the correct limit.
%
% Reference: S. Halfin, W. Whitt (1981). Heavy-traffic limits for queues with
% many exponential servers. Operations Research 29(3), 567-588.
%
% See also QSYS_MMK_QED, QSYS_MMK_QED_STAFFING.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

alpha = ones(size(beta));
pos = beta > 0;
b = beta(pos);
phi = exp(-b.^2/2) / sqrt(2*pi);
Phi = erfc(-b/sqrt(2)) / 2;
alpha(pos) = phi ./ (phi + b .* Phi);
end
