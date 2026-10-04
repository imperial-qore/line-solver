function res = laplace_invert_pdf(F, tset, method, n)
% RES = LAPLACE_INVERT_PDF(F, TSET, METHOD, N)
%
% Invert the transform of a DENSITY on a grid. Values at t <= 0 are 0 and the
% result is clamped at zero. Twin of the native Python
% api.lti.laplace_invert_pdf.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(method)
    method = 'euler';
end
if nargin < 4
    n = [];
end

res = zeros(size(tset));
if strcmpi(method, 'weeks')
    % One coefficient set serves the whole grid; this is the point of Weeks.
    if isempty(n), n = 200; end
    [sigma, b, q] = laplace_weeks_scaling(F, n);
    res = laplace_invert_weeks(F, tset, q, sigma, b);
else
    for i = 1:numel(tset)
        if tset(i) > 0
            res(i) = laplace_invert(F, tset(i), method, n);
        end
    end
end
res = max(res, 0);
end
