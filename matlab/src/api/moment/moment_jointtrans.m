function B = moment_jointtrans(A, edge)
% B = moment_jointtrans(A, edge)
%
% Applies the conversion matrix of one edge of the house of moments along every
% dimension of a joint moment array. This is the separable (Kronecker product)
% form shared by all the joint conversions except the cumulant and the central
% ones.
%
% Input:
%   A: joint moment array of size (n_1+1)x...x(n_d+1)
%   edge: edge label accepted by moment_housematrix
%
% Output:
%   B: array of the same size as A holding the converted joint moments
%
% Example:
%   f = moment_jointtrans(m, 'factorial_from_raw');
%
% Reference:
% A. Heindl and A. van de Liefvoort. Moment conversions for discrete
% distributions. PMCCS, 2003.

sz = moment_tensorsize(A);
B = A;
for mode = 1:numel(sz)
    B = moment_tensortrans(B, moment_housematrix(edge, sz(mode)-1), mode);
end
end
