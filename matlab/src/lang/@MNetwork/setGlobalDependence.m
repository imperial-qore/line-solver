function setGlobalDependence(self, phi, peakRate, wireCutoff)
% SETGLOBALDEPENDENCE(PHI, PEAKRATE, WIRECUTOFF)
%
% Declares a globally state-dependent service-rate scaling phi(n), where n is
% the FULL (nstations x nclasses) matrix of per-class populations, i.e. the
% whole network state rather than the population local to one station. This is
% the primitive behind Whittle networks: a rate that reads the entire state and,
% when it satisfies the balance property
%
%   phi_s(n) phi_t(n-e_s) = phi_t(n) phi_s(n-e_t),
%
% yields the reversible product form pi(n) ~ Phi(n) prod_s rho_s^{n_s} and
% insensitivity to the service-time distribution beyond its mean. It also
% expresses bandwidth-sharing allocations, where one route holds several links
% simultaneously and no per-station scaling can reproduce the coupling.
%
% PHI returns either a scalar (broadcast to every station and class), an
% (nstations x 1) column (per station, broadcast over classes) or an
% (nstations x nclasses) matrix. The effective service rate of class r at
% station i is its base rate times phi(i,r), composing multiplicatively with any
% setLoadDependence / setClassDependence / setJointDependence already declared.
%
% PEAKRATE is REQUIRED and normalizes utilization as Util = T*S/peak, matching
% the T*S/c convention of ordinary multiserver stations. Pass a scalar
% (identical peak everywhere) or an (nstations x nclasses) matrix.
%
% WIRECUTOFF (optional, default 10) is the per-slot OPEN-class truncation used
% when PHI is materialized onto the JSON wire by LINEMODEL_SAVE; closed classes
% are tabulated up to their own population instead. It plays no part in solving
% and exists because PHI, being a handle, cannot cross a language boundary: the
% writer needs to know how far the lattice extends, and a handle defined only up
% to some population (as a balanced-fairness recursion is) must say so. Set it to
% the same cutoff the model is solved at.
%
% Unlike setClassDependence (product-form beta_{i,r}), a global dependence is
% NON-product-form in general and only SolverCTMC declares support for it.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~isa(phi,'function_handle')
    line_error(mfilename, 'Global dependence must be specified through a function handle.');
end
if nargin < 3 || isempty(peakRate)
    line_error(mfilename, 'Global dependence requires an explicit peak rate: setGlobalDependence(phi, peakRate). Pass a scalar or an (nstations x nclasses) matrix.');
end
if ~isnumeric(peakRate) || any(peakRate(:) <= 0)
    line_error(mfilename, 'peakRate must be positive.');
end
if nargin < 4 || isempty(wireCutoff)
    wireCutoff = 10;
end
if ~isnumeric(wireCutoff) || ~isscalar(wireCutoff) || wireCutoff < 1 || wireCutoff ~= round(wireCutoff)
    line_error(mfilename, 'wireCutoff must be a positive integer scalar.');
end

M = getNumberOfStations(self);
K = getNumberOfClasses(self);

% Probe the handle now so a wrong output shape is refused at declaration time
% rather than midway through state-space generation.
probe = {zeros(M,K), ones(M,K)};
for p = 1:numel(probe)
    try
        v = phi(probe{p});
    catch ME
        line_error(mfilename, sprintf('The global dependence handle failed to evaluate on a %dx%d population matrix: %s', M, K, ME.message));
    end
    if ~isnumeric(v) || any(~isfinite(v(:))) || any(v(:) < 0)
        line_error(mfilename, 'The global dependence handle must return finite nonnegative numeric scalings.');
    end
    if ~(isscalar(v) || isequal(size(v),[M,1]) || isequal(size(v),[M,K]))
        line_error(mfilename, sprintf('The global dependence handle must return a scalar, an (%d x 1) column or an (%d x %d) matrix; it returned a %s array.', M, M, K, mat2str(size(v))));
    end
end

if isscalar(peakRate)
    peak = peakRate * ones(M,K);
elseif isequal(size(peakRate),[M,1])
    peak = repmat(peakRate(:),1,K);
elseif isequal(size(peakRate),[M,K])
    peak = peakRate;
else
    line_error(mfilename, sprintf('peakRate must be a scalar, an (%d x 1) column or an (%d x %d) matrix.', M, M, K));
end

self.gdScaling = phi;
self.gdScalingPeak = peak;
self.gdScalingCutoff = double(wireCutoff);
self.resetStruct();
end
