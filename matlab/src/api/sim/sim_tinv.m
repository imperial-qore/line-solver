function t = sim_tinv(p, nu)
% OA_TINV Quantile function of Student's t distribution.
%
% T = OA_TINV(P, NU) returns the P-quantile of the t distribution with NU
% degrees of freedom. P is a scalar or array, NU a positive scalar.
%
% Implemented on the incomplete beta inverse, which is in base MATLAB, so the
% output-analysis routines in this folder do not pull in the Statistics and
% Machine Learning Toolbox. The identity used is
%   P(|T| > t) = betainc(nu/(nu+t^2), nu/2, 1/2),
% inverted for the two-sided tail 2(1-P) and mapped back with
%   t = sqrt(nu (1-z)/z),   z = betaincinv(2(1-P), nu/2, 1/2).
% Agreement with the toolbox TINV is to within 1e-10 over the range used here.
%
% Examples:
%   sim_tinv(0.975, 10)   % 2.2281
%   sim_tinv(0.975, 29)   % 2.0452
%
% See also OA_NORMCDF, OA_NORMINV
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~isscalar(nu) || ~isreal(nu) || nu <= 0
    line_error(mfilename, 'nu must be a positive real scalar');
end
if any(p(:) < 0) || any(p(:) > 1)
    line_error(mfilename, 'Probability must lie in [0,1]');
end

t = zeros(size(p));
for i = 1:numel(p)
    pi_ = p(i);
    if pi_ == 0.5
        t(i) = 0;
    elseif pi_ <= 0
        t(i) = -Inf;
    elseif pi_ >= 1
        t(i) = Inf;
    else
        % reflect the lower half onto the upper half, the law is symmetric
        flip = pi_ < 0.5;
        if flip
            pi_ = 1 - pi_;
        end
        z = betaincinv(2 * (1 - pi_), nu / 2, 0.5);
        ti = sqrt(nu * (1 - z) / z);
        if flip
            ti = -ti;
        end
        t(i) = ti;
    end
end
end
