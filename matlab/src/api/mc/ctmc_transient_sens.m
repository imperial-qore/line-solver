function [dpi, pi, t] = ctmc_transient_sens(Q, dQ, pi0, t0, t1)
% [DPI, PI, T] = CTMC_TRANSIENT_SENS(Q, DQ, PI0, T0, T1)
%
% Sensitivity of the transient distribution of a CTMC to a scalar parameter
% theta, given the generator Q and its derivative DQ = dQ/dtheta.
%
% Differentiating the forward equations d pi(t)/dt = pi(t) Q with respect to
% theta, and assuming the initial vector does not depend on theta, gives
%
%   d/dt (dpi(t)/dtheta) = (dpi(t)/dtheta) Q + pi(t) (dQ/dtheta),
%   dpi(0)/dtheta = 0,
%
% i.e. Trivedi and Bobbio (2017), Eq. (9.82). The state and its sensitivity
% are integrated as one augmented system of size 2n, since the sensitivity
% equation is driven by pi(t) and the two cannot be advanced separately.
%
% @param Q Generator matrix (n x n)
% @param dQ Derivative of the generator with respect to theta (n x n)
% @param pi0 Initial distribution (1 x n); uniform if empty
% @param t0 Initial time; 0 if omitted
% @param t1 Final time
% @return dpi Sensitivity of the distribution at each time point (length(T) x n)
% @return pi Distribution at each time point (length(T) x n)
% @return t Column vector of time points
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

n = size(Q, 1);
if nargin == 3
    t1 = pi0;
    t0 = 0;
    pi0 = [];
end
if nargin == 4
    t1 = t0;
    t0 = 0;
end
if isempty(pi0)
    pi0 = ones(1, n) / n;
end
if any(size(dQ) ~= size(Q))
    line_error(mfilename, 'dQ must have the same size as Q');
end

% Augmented state v = [pi, dpi], with dpi(0) = 0 since pi(0) does not depend
% on theta
v0 = [pi0(:); zeros(n, 1)];

[t, v] = ode23(@augmentedode, [t0, t1], v0);

pi = v(:, 1:n);
dpi = v(:, (n+1):(2*n));

    function dvdt = augmentedode(~, v)
        % DVDT = AUGMENTEDODE(T, V)

        p = v(1:n)';
        s = v((n+1):(2*n))';
        dpdt = p * Q;
        dsdt = s * Q + p * dQ;
        dvdt = [dpdt(:); dsdt(:)];
    end

end
