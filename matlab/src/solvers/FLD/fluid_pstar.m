function [usePnorm, pstar] = fluid_pstar(method, options, M)
% [USEPNORM, PSTAR] = FLUID_PSTAR(METHOD, OPTIONS, M)
%
% @brief The p-norm smoothing exponent in force for a matrix-family fluid drift.
%
% Ruuskanen et al., PEVA 151 (2021), eq. (26)-(27) replace the hard min(n,c) of
% the matrix drift by ghat = (1 + (n/c)^p)^(-1/p), so p selects WHICH DRIFT is
% integrated rather than tuning one.
%
% TWO WAYS IN, AND THE METHOD NAME IS ONE OF THEM. options.pstar or
% options.config.pstar switch the smoothing on under any matrix-family method,
% which is what a caller naming an exponent means. Failing that the method
% 'pnorm' supplies its own default of 20, the exponent whose behaviour matches
% softmin at alpha = 20 and the one SolverFluid.java and the C++ FluidOptions
% already default to. Without it the name 'pnorm' integrated the UNSMOOTHED
% matrix drift, so the label named a method that never ran.
%
% One rule, two callers: SOLVER_FLUID_MATRIX integrates the drift this selects
% and SOLVER_FLUID_SYMODES exports it, so the exported system and the
% integrated one cannot drift apart.
%
% @param method the requested method name, qualified ('fluid.pnorm') or not
% @param options solver options, read for pstar and config.pstar
% @param M number of stations, the length of the returned vector
% @return usePnorm true when the p-norm drift is the one to build
% @return pstar M-by-1 exponent per station, empty when USEPNORM is false

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

pstarVal = [];
if isfield(options,'pstar') && ~isempty(options.pstar)
    pstarVal = options.pstar;
elseif isfield(options,'config') && isfield(options.config,'pstar') ...
        && ~isempty(options.config.pstar)
    pstarVal = options.config.pstar;
elseif ischar(method) || isstring(method)
    if strcmp(strrep(char(method), 'fluid.', ''), 'pnorm')
        pstarVal = 20;
    end
end
usePnorm = ~isempty(pstarVal);
if usePnorm
    if isscalar(pstarVal)
        pstarVal = pstarVal * ones(M, 1);
    end
    pstar = pstarVal(:);
else
    pstar = [];
end
end
