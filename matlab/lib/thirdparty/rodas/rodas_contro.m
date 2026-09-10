function v = rodas_contro(dense, i, x)
% RODAS_CONTRO  CONTRO: the third-order interpolant RODAS carries per step
%
%   V = RODAS_CONTRO(DENSE, I, X)
%
%   PORTED THIRD-PARTY CODE (rodas.f, Hairer & Wanner) -- see RODAS_CORE.
%
%   Component I of the solution at X, valid over the step just accepted. A
%   Rosenbrock method chooses its step from the local error, so its accepted
%   points are wherever the stiffness put them and never the ones a caller
%   asked for; this is what RODAS carries so that an arbitrary output grid
%   costs no extra step and no interpolation of the caller's own.
%
%   DENSE is the struct handed to the SOLOUT callback. On the FIRST call --
%   made before any step, with NR = 1 -- CONT holds no coefficients yet and the
%   state at that point is Y itself, so callers must not interpolate there.
%
%   See also RODAS_CORE, RODAS.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

n = dense.n;
cont = dense.cont;
s = (x - dense.xold)/dense.h;
v = cont(i)*(1 - s) + s*(cont(i+n) + (1 - s)*(cont(i+2*n) + s*cont(i+3*n)));
end
