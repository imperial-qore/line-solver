function res = laplace_invert_cdf(F, tset, method, n)
% RES = LAPLACE_INVERT_CDF(F, TSET, METHOD, N)
%
% Invert the transform of a CDF on a grid. F IS THE TRANSFORM OF THE DENSITY,
% not of the CDF: this function forms F(s)/s itself. Values at t <= 0 are 0;
% the result is clipped to [0,1] and made monotone, since neither property is
% guaranteed by a truncated inversion. Twin of the native Python
% api.lti.laplace_invert_cdf.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(method)
    method = 'euler';
end
if nargin < 4
    n = [];
end

Fcdf = @(s) local_ratio(F, s);
res = zeros(size(tset));
if strcmpi(method, 'weeks')
    if isempty(n), n = 200; end
    [sigma, b, q] = laplace_weeks_scaling(Fcdf, n);
    res = laplace_invert_weeks(Fcdf, tset, q, sigma, b);
else
    for i = 1:numel(tset)
        if tset(i) > 0
            res(i) = laplace_invert(Fcdf, tset(i), method, n);
        end
    end
end
res = min(max(res, 0), 1);
res = cummax(res);
end

function v = local_ratio(F, s)
if abs(s) < 1e-15
    v = 1.0;
else
    v = F(s) / s;
end
end
