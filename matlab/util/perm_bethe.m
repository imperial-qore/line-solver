function val = perm_bethe(A, options)
% VAL = PERM_BETHE(A)
% VAL = PERM_BETHE(A, OPTIONS)
%
% Bethe (sum product) approximation of the permanent of a nonnegative matrix.
% Two families of messages, the right going messages R and the left going
% messages L, are iterated to a fixed point on the elementwise square root of
% A, and the Bethe free energy of that fixed point is exponentiated.
%
% Twin of jline.lib.perm.BethePermanent and of the python
% line_solver.api.perm.BethePermanent. All three exclude the DIAGONAL entry
% from each message denominator rather than the target entry; see
% _kb/03-api-layer.md for why that is kept as is in every codebase.
%
% Input:
%   A       - nonnegative square matrix
%   OPTIONS - optional struct with fields
%             epsilon      squared message change below which the iteration
%                          stops (default 1e-3)
%             maxIteration maximum number of message passing sweeps
%                          (default 200000)
%
% Output:
%   VAL - approximation of the permanent
%
% The Bethe permanent is a LOWER BOUND of the permanent for a nonnegative
% matrix, inside the Gurvits sqrt(2)^n factor. The gap grows with n and is not
% an error: do not read it as a defect, and do not expect it to cancel when
% several estimates are normalized against each other, since it is
% state-dependent.
%
% Reference:
%   P. O. Vontobel, "The Bethe Permanent of a Nonnegative Matrix", IEEE
%   Transactions on Information Theory, 59(3), 2013.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    options = struct();
end
if ~isfield(options, 'epsilon'), options.epsilon = 1e-3; end
if ~isfield(options, 'maxIteration'), options.maxIteration = 200000; end

if any(A(:) < 0)
    line_error(mfilename, 'Matrix must be non-negative.');
end
if size(A,1) ~= size(A,2)
    line_error(mfilename, 'Matrix must be square.');
end

% MINVAL guards the LOGARITHMS of message products below against underflow. It
% is deliberately NOT applied to A: flooring the input is what fabricates a
% permanent of n!*eps where the truth is zero, so a non-positive entry is
% refused outright instead. See PERM_REQUIRE_SUPPORT.
MINVAL = 2.220446049250314e-16;

n = size(A, 1);
if n == 0
    val = 1;
    return
end

perm_require_support(A, mfilename);

S = sqrt(A);

rpast = ones(n, n);
lpast = ones(n, n);
[r, l] = sub_update(S, lpast, n);

iteration = 0;
while sum(sum((rpast-r).^2 + (lpast-l).^2)) > options.epsilon && iteration < options.maxIteration
    iteration = iteration + 1;
    rpast = r;
    lpast = l;
    [r, l] = sub_update(S, lpast, n);
end

term1 = max(sum(S.*l, 2), MINVAL);
term2 = max(sum(S.*r, 1), MINVAL);
term3 = max(r.*l + 1, MINVAL);

logval = sum(log(term1)) + sum(log(term2)) - sum(sum(log(term3)));
val = exp(logval);
if ~isfinite(val)
    val = 0;
end
end

function [r1, l1] = sub_update(S, l, n)
% One sweep of the sum product recursion. The denominators drop the diagonal
% term, matching the JAR and python twins.
rowTerms = S .* l;
denomR = sum(rowTerms, 2) - diag(rowTerms);
r1 = S ./ repmat(denomR, 1, n);

colTerms = S .* r1;
denomL = sum(colTerms, 1) - diag(colTerms)';
l1 = S ./ repmat(denomL, n, 1);
end
