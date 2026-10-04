function bool = ctmc_isfeasible(Q, tol)
% BOOL = CTMC_ISFEASIBLE(Q, TOL)
%
% True when Q is a valid infinitesimal generator: square, non-negative
% off-diagonal entries, non-positive diagonal, and zero row sums, each up to
% TOL (default 1e-10). Twin of the Python api.mc.ctmc_isfeasible; note that
% dtmc_isfeasible instead returns a precision level rather than a flag.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    tol = 1e-10;
end

bool = false;
if size(Q,1) ~= size(Q,2) || isempty(Q)
    return
end
Q = full(Q);
offDiag = Q - diag(diag(Q));
if any(offDiag(:) < -tol)
    return
end
if any(diag(Q) > tol)
    return
end
if any(abs(sum(Q,2)) > tol)
    return
end
bool = true;
end
