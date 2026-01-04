function lst = aoi_lst_exp(mu)
%AOI_LST_EXP Laplace-Stieltjes transform for exponential distribution
%
% lst = aoi_lst_exp(mu)
%
% Returns a function handle for the LST of an exponential distribution
% with rate mu.
%
% Parameters:
%   mu (double): Rate parameter (mean = 1/mu)
%
% Returns:
%   lst (function_handle): LST function @(s) mu/(mu+s)
%
% The LST of Exp(mu) is: H*(s) = mu / (mu + s)
%
% See also: aoi_lst_erlang, aoi_lst_det, aoi_lst_ph

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if mu <= 0
    line_error(mfilename, 'Rate mu must be positive');
end

lst = @(s) mu ./ (mu + s);

end
