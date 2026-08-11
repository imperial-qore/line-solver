function d2 = da_traffic_superpos(lambda, a2)
% D2 = DA_TRAFFIC_SUPERPOS(LAMBDA, A2)
%
% Asymptotic-method superposition of independent flows with rates LAMBDA
% and squared coefficients of variation A2: returns the rate-weighted SCV
% mixture of the merged flow (Whitt's QNA stationary-interval formula).
% Entries with non-finite rates are ignored.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

a2 = a2(isfinite(lambda));
lambda = lambda(isfinite(lambda));
d2 = a2(:)' * lambda(:) / sum(lambda);
end
