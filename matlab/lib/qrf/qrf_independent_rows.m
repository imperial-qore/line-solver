function keep = qrf_independent_rows(A, tol)
% QRF_INDEPENDENT_ROWS  Indices of a maximal independent subset of A's rows.
%
%   KEEP = QRF_INDEPENDENT_ROWS(A) selects by column-pivoted QR of A', so the
%   retained rows are the numerically best-conditioned independent set rather
%   than simply the first ones encountered.
%
%   The QRF equality block is heavily redundant, carrying about twice as many
%   rows as its rank, and both quadprog (MATLAB) and SLSQP (Python) degrade on
%   the rank-deficient system. Dropping dependent rows changes no feasible
%   point: they are exact linear combinations of the retained ones. This is the
%   MATLAB counterpart of qrf_noblo_common.independent_rows.

if isempty(A)
    keep = zeros(0,1);
    return
end
[~, R, P] = qr(full(A)', 0);
d = abs(diag(R));
if isempty(d)
    keep = zeros(0,1);
    return
end
if nargin < 2 || isempty(tol)
    tol = max(size(A)) * eps * d(1);
end
keep = sort(P(1:sum(d > tol)));
keep = keep(:);
end
