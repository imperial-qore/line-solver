function sz = moment_tensorsize(A)
% sz = moment_tensorsize(A)
%
% Size of a joint moment array with the trailing singleton dimensions removed.
%
% MATLAB cannot represent a trailing singleton dimension, so a column vector of
% n+1 elements and an (n+1)x1 joint moment array of a degenerate second class
% are the same object. The convention taken throughout the joint conversions is
% the first one: trailing singletons are dropped, so that a column vector is
% treated as the d = 1 case.
%
% Input:
%   A: joint moment array
%
% Output:
%   sz: row vector of the extents, with the trailing singleton dimensions
%       removed
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003.

sz = size(A);
while numel(sz) > 1 && sz(end) == 1
    sz(end) = [];
end
end
