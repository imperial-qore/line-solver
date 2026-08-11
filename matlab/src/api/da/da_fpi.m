function [x, it, converged] = da_fpi(iterfun, x0, options)
% [X,IT,CONVERGED] = DA_FPI(ITERFUN, X0, OPTIONS)
%
% Generic damped successive-substitution driver for decomposition-
% aggregation (DA) fixed-point iterations. Each call to ITERFUN performs
% one DA sweep: solve the isolated submodels given the current coupling
% iterate X, exchange flows or rates, and return the updated iterate.
%
% ITERFUN: function handle [XNEW, XREF] = ITERFUN(X, IT) evaluating one DA
%          sweep from iterate X at iteration count IT. XREF is the baseline
%          for the convergence test; return XREF = X for a standard
%          successive-substitution test, or a mid-sweep checkpoint when the
%          method compares against a renormalized iterate.
% X0:      initial iterate (any numeric array).
% OPTIONS: solver options struct; uses iter_max and iter_tol, plus the
%          optional fields config.da_damping in (0,1] (default 1, i.e.
%          undamped), config.da_norm (function handle mapping the iterate
%          difference to a scalar, default @(d) max(abs(d(:))); a two-
%          argument handle is called as da_norm(xnew, xref) instead, e.g.
%          for relative-difference tests), config.da_miniter (default 1):
%          convergence is not tested before this sweep count, and
%          config.da_nanstop (default false): when true, a NaN convergence
%          measure terminates the iteration, replicating legacy while-loop
%          drivers whose "continue while delta > tol" test exits on NaN;
%          when false a NaN measure keeps iterating, as in legacy
%          "break if delta < tol" drivers.
%
% Returns the final iterate X, the number IT of sweeps executed, and a
% CONVERGED flag (false if the iteration stopped at iter_max).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

omega = 1;
normfun = @(d) max(abs(d(:)));
nanstop = false;
miniter = 1;
if isfield(options, 'config')
    if isfield(options.config, 'da_damping') && ~isempty(options.config.da_damping)
        omega = options.config.da_damping;
    end
    if isfield(options.config, 'da_norm') && ~isempty(options.config.da_norm)
        normfun = options.config.da_norm;
    end
    if isfield(options.config, 'da_nanstop') && ~isempty(options.config.da_nanstop)
        nanstop = options.config.da_nanstop;
    end
    if isfield(options.config, 'da_miniter') && ~isempty(options.config.da_miniter)
        miniter = options.config.da_miniter;
    end
end
twoargnorm = nargin(normfun) == 2;

x = x0;
it = 0;
converged = false;
for it = 1:options.iter_max
    [xnew, xref] = iterfun(x, it);
    if omega ~= 1
        xnew = (1 - omega) * xref + omega * xnew;
    end
    if twoargnorm
        delta = normfun(xnew, xref);
    else
        delta = normfun(xnew - xref);
    end
    x = xnew;
    if it >= miniter
        if delta < options.iter_tol
            converged = true;
            break
        elseif nanstop && isnan(delta)
            break
        end
    end
end
end
