function lst = aoi_lst_det(d)
%AOI_LST_DET Laplace-Stieltjes transform for deterministic (constant) distribution
%
% lst = aoi_lst_det(d)
%
% Returns a function handle for the LST of a deterministic distribution
% with constant value d.
%
% Parameters:
%   d (double): Constant value (must be positive)
%
% Returns:
%   lst (function_handle): LST function @(s) exp(-s*d)
%
% The LST of a deterministic random variable X = d is: H*(s) = exp(-s*d)
%
% See also: aoi_lst_exp, aoi_lst_erlang, aoi_lst_ph

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if d <= 0
    line_error(mfilename, 'Constant d must be positive');
end

lst = @(s) exp(-s .* d);

end
