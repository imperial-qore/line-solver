function res = laplace_invert_cme(F, t, maxFnEvals)
% RES = LAPLACE_INVERT_CME(F, T, MAXFNEVALS)
%
% Invert a Laplace transform at T by the Concentrated Matrix Exponential
% method of Horvath, Horvath and Telek, using the pre-computed coefficient
% tables shipped with LINE. Delegates to the vendored MATLAB_ILT, which owns
% the tables; this wrapper exists so that api/lti presents the same surface in
% all four codebases.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(maxFnEvals)
    maxFnEvals = 25;
end
if t <= 0
    line_error(mfilename, 'The Laplace inversion time point must be positive.');
end

res = matlab_ilt(F, t, maxFnEvals, 'cme');
res = res(1);
end
