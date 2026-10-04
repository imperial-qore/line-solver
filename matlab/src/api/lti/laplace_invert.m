function res = laplace_invert(F, t, method, n)
% RES = LAPLACE_INVERT(F, T, METHOD, N)
%
% Invert a Laplace transform at a single time point T. METHOD is one of
% 'euler' (default), 'talbot', 'gaver-stehfest', 'cme' or 'weeks'. N is the
% number of terms and defaults per method: 41 euler, 32 talbot, 12
% gaver-stehfest, 25 cme, 200 weeks (the p0 of the coefficient search).
%
% 'weeks' recomputes its coefficients on every call, which wastes the one
% property that method has. Call LAPLACE_WEEKS_SCALING once and pass the
% coefficients to LAPLACE_INVERT_WEEKS when inverting on a grid.
%
% Twin of the native Python api.lti.laplace_invert.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(method)
    method = 'euler';
end
if nargin < 4
    n = [];
end

switch lower(method)
    case 'euler'
        if isempty(n), n = 41; end
        res = laplace_invert_euler(F, t, n);
    case 'talbot'
        if isempty(n), n = 32; end
        res = laplace_invert_talbot(F, t, n);
    case {'gaver-stehfest','gaver_stehfest','gaver'}
        if isempty(n), n = 12; end
        res = laplace_invert_gaver_stehfest(F, t, n);
    case 'cme'
        if isempty(n), n = 25; end
        res = laplace_invert_cme(F, t, n);
    case 'weeks'
        if isempty(n), n = 200; end
        [sigma, b, q] = laplace_weeks_scaling(F, n);
        res = laplace_invert_weeks(F, t, q, sigma, b);
    otherwise
        line_error(mfilename, sprintf('Unknown Laplace inversion method: %s. Supported: euler, talbot, gaver-stehfest, cme, weeks.', method));
end
end
