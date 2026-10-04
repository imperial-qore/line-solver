function z = sim_norminv(p)
% OA_NORMINV Standard normal quantile function.
%
% Z = OA_NORMINV(P) returns Phi^{-1}(P) elementwise for P in (0,1).
% Implemented on ERFINV, which is in base MATLAB, so the output-analysis
% routines in this folder do not pull in the Statistics and Machine Learning
% Toolbox. P = 0 gives -Inf and P = 1 gives +Inf.
%
% Examples:
%   sim_norminv(0.975)   % 1.9600
%   sim_norminv(0.5)     % 0
%
% See also OA_NORMCDF, OA_TINV
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if any(p(:) < 0) || any(p(:) > 1)
    line_error(mfilename, 'Probability must lie in [0,1].');
end
z = sqrt(2) * erfinv(2 * p - 1);
end
