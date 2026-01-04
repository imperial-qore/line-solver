function lst = aoi_lst_erlang(k, mu)
%AOI_LST_ERLANG Laplace-Stieltjes transform for Erlang distribution
%
% lst = aoi_lst_erlang(k, mu)
%
% Returns a function handle for the LST of an Erlang-k distribution
% with rate parameter mu.
%
% Parameters:
%   k (integer): Shape parameter (number of phases)
%   mu (double): Rate parameter per phase (mean = k/mu)
%
% Returns:
%   lst (function_handle): LST function @(s) (mu/(mu+s))^k
%
% The LST of Erlang(k, mu) is: H*(s) = (mu / (mu + s))^k
%
% See also: aoi_lst_exp, aoi_lst_det, aoi_lst_ph

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if k < 1 || k ~= round(k)
    line_error(mfilename, 'Shape k must be a positive integer');
end
if mu <= 0
    line_error(mfilename, 'Rate mu must be positive');
end

lst = @(s) (mu ./ (mu + s)).^k;

end
