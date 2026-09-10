function p = sim_normcdf(z)
% OA_NORMCDF Standard normal cumulative distribution function.
%
% P = OA_NORMCDF(Z) returns Phi(Z) elementwise. Implemented on ERFC, which is
% in base MATLAB, so the output-analysis routines in this folder do not pull in
% the Statistics and Machine Learning Toolbox.
%
% Examples:
%   sim_normcdf(0)       % 0.5
%   sim_normcdf(1.96)    % 0.9750
%
% See also OA_NORMINV, OA_TINV
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

p = 0.5 * erfc(-z ./ sqrt(2));
end
