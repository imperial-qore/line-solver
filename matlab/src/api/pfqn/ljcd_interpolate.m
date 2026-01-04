function Xval = ljcd_interpolate(nvec, cutoffs, table, K)
% LJCD_INTERPOLATE Multi-linear interpolation for LJCD scaling tables
%
% Xval = LJCD_INTERPOLATE(nvec, cutoffs, table, K)
%
% Performs multi-linear interpolation of a throughput value from an LJCD
% scaling table for non-integer population vectors.
%
% Parameters:
%   nvec    - [n1, n2, ..., nK] continuous population vector (clamped to cutoffs)
%   cutoffs - [N1, N2, ..., NK] per-class cutoffs
%   table   - Linearized throughput table (1-D array indexed by ljd_linearize)
%   K       - Number of classes
%
% Returns:
%   Xval    - Interpolated throughput value
%
% For K classes, interpolates between 2^K corner points of the hypercube
% containing the population vector using multi-linear interpolation.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Get floor and ceiling for each dimension
nFloor = floor(nvec);
nCeil = ceil(nvec);

% Clamp to valid range
nFloor = max(0, min(nFloor, cutoffs));
nCeil = max(0, min(nCeil, cutoffs));

% Compute fractional parts (weights)
frac = nvec - nFloor;

% Handle edge case: all integer values
if all(abs(frac) < 1e-10)
    idx = ljd_linearize(nFloor, cutoffs);
    if idx <= length(table)
        Xval = table(idx);
    else
        Xval = 0;
    end
    return;
end

% Multi-linear interpolation over 2^K corners
Xval = 0;
numCorners = 2^K;

for corner = 0:(numCorners - 1)
    % Build corner point: bit i determines floor (0) or ceil (1) for class i
    cornerPoint = nFloor;
    weight = 1.0;

    for i = 1:K
        if bitget(corner, i)
            % Use ceiling for this dimension
            cornerPoint(i) = nCeil(i);
            weight = weight * frac(i);
        else
            % Use floor for this dimension
            weight = weight * (1 - frac(i));
        end
    end

    % Look up table value at corner point
    idx = ljd_linearize(cornerPoint, cutoffs);
    if idx <= length(table)
        Xval = Xval + weight * table(idx);
    end
end

end
