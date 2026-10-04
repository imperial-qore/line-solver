function [tout, xtraj, sj] = cache_sojourn_ode(drift, tspan, x0, wgen, odeopt)
% [TOUT, XTRAJ, SJ] = CACHE_SOJOURN_ODE(DRIFT, TSPAN, X0, WGEN, ODEOPT)
%
% Integrates the cache mean-field drift dx/dt = DRIFT(t,x) over TSPAN from X0
% and, when WGEN is given, the SOJOURN-WEIGHTED AVERAGE of the trajectory
%
%     XBAR = int_{t0}^{t1} x(t) g(t) dt / int_{t0}^{t1} g(t) dt,
%     g(t) = phi(t) * c,   dphi/dt = phi * A,   phi(t0) = WGEN.phi0,
%
% i.e. the average of x against the density g of a phase-type clock. The
% integral is carried by the integrator itself as extra state components
% (phi, int g x, int g), so its value is accurate to the ODE tolerance and does
% NOT depend on where the solver chose to emit output points. A left or right
% Riemann sum over the output grid did, and that is what made the adaptive
% (MATLAB) and fixed-grid (JAR/Python) editions of the ENV cache mean field
% disagree in the third digit.
%
% WGEN: struct with A (nph x nph), phi0 (1 x nph), c (nph x 1). For a stage
%       whose holding time is the MAP {D0,D1} and whose drift runs in a time
%       unit that is LAM times real time: A = D0/LAM, c = -D0*1/LAM and
%       phi0 = map_pie * expm(A*t0), which makes int g = F(t1/LAM) - F(t0/LAM)
%       with F the holding-time CDF.
%
% Returns the grid TOUT and occupancy XTRAJ (dim x nt) of the x part, and SJ
% with fields XBAR (dim x 1, empty when int g = 0) and WTOT (= int g).
% SJ is [] when WGEN is empty, and the integration is then the plain drift.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

x0 = x0(:);
if nargin < 4 || isempty(wgen)
    [tout, xtraj] = ode15s(drift, tspan, x0, odeopt);
    xtraj = xtraj';
    sj = [];
    return
end
n = numel(x0);
A = wgen.A;
c = wgen.c(:);
nph = size(A, 1);
phi0 = wgen.phi0(:);
At = A';
iphi = n + (1:nph);
iI = n + nph + (1:n);
iW = n + nph + n + 1;
    function dz = rhs(t, z)
        x = z(1:n);
        phi = z(iphi);
        g = phi' * c;
        dz = [drift(t, x); At * phi; g * x; g];
    end
z0 = [x0; phi0; zeros(n, 1); 0];
[tout, Z] = ode15s(@rhs, tspan, z0, odeopt);
xtraj = Z(:, 1:n)';
zend = Z(end, :)';
sj = struct('xbar', [], 'wtot', zend(iW));
if zend(iW) > 0
    sj.xbar = zend(iI) / zend(iW);
end
end
