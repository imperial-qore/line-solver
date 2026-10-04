function lst = aoi_lst_ph(alpha, T)
%AOI_LST_PH Laplace-Stieltjes transform for phase-type distribution
%
% lst = aoi_lst_ph(alpha, T)
%
% Returns a function handle for the LST of a phase-type (PH) distribution
% with initial probability vector alpha and sub-generator matrix T.
%
% Parameters:
%   alpha (row vector): Initial probability vector (1 x n)
%   T (matrix): Sub-generator matrix (n x n)
%
% Returns:
%   lst (function_handle): LST function @(s) alpha * inv(s*I - T) * (-T*e)
%
% The LST of PH(alpha, T) is: H*(s) = alpha * (s*I - T)^{-1} * t
% where t = -T * ones(n,1) is the exit rate vector.
%
% For numerical stability, we use: H*(s) = alpha * ((s*I - T) \ t)
%
% See also: aoi_lst_exp, aoi_lst_erlang, aoi_lst_det

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Validate inputs
alpha = alpha(:)';  % Ensure row vector
n = length(alpha);

if size(T, 1) ~= n || size(T, 2) ~= n
    line_error(mfilename, 'T must be %d x %d to match alpha', n, n);
end

% Exit rate vector
t = -T * ones(n, 1);

% Identity matrix
I = eye(n);

% Return function handle
lst = @(s) ph_lst_eval(s, alpha, T, t, I);

end

function val = ph_lst_eval(s, alpha, T, t, I)
%PH_LST_EVAL Evaluate PH LST at point(s) s

    if isscalar(s)
        val = alpha * ((s * I - T) \ t);
    else
        % Handle array input
        val = zeros(size(s));
        for i = 1:numel(s)
            val(i) = alpha * ((s(i) * I - T) \ t);
        end
    end
end
