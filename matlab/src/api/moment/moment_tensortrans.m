function B = moment_tensortrans(A, T, mode)
% B = moment_tensortrans(A, T, mode)
%
% Applies a conversion matrix along one dimension of a joint moment array.
%
% This is the mode product of the array with the matrix: every fibre of A along
% the given dimension is replaced by T times that fibre. Applying it once per
% dimension realises the Kronecker product of the univariate conversions, which
% is the structure of every separable edge of the house of moments.
%
% Input:
%   A: joint moment array of size (n_1+1)x...x(n_d+1)
%   T: (n_mode+1)x(n_mode+1) conversion matrix
%   mode: dimension to transform (1 <= mode <= d)
%
% Output:
%   B: array of the same size as A
%
% Example:
%   B = moment_tensortrans(A, moment_stirling1(size(A,1)-1), 1);
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003.

sz = moment_tensorsize(A);
d = numel(sz);
if ~isscalar(mode) || mode < 1 || mode > d || mode ~= round(mode)
    line_error(mfilename,'The mode must be a dimension index in 1,...,ndims(A).');
end
if size(T,2) ~= sz(mode)
    line_error(mfilename,'The matrix T must have as many columns as the extent of the transformed dimension.');
end
% MATLAB has no 1-D array, so a d = 1 input is reshaped to (n+1)x1 and the
% permutation is padded to the two dimensions permute insists on.
dfull = max(d, 2);
szf = [sz, ones(1, dfull - d)];
perm = [mode, setdiff(1:dfull, mode)];
Ap = permute(reshape(A, szf), perm);
szp = size(Ap);
szp = [szp, ones(1, dfull - numel(szp))];
Bp = reshape(T * reshape(Ap, szf(mode), []), [size(T,1), szp(2:end)]);
B = reshape(ipermute(Bp, perm), size(A));
end
