function C = ldesTrajCell(raw, M, R)
% C = LDESTRAJCELL(RAW, M, R)
% Normalize a jsondecode'd transient trajectory (logical nesting
% [station][class][row][2]) into an M x R cell array C where C{i,r} is the
% (nRows x 2) [value, time] matrix for class r at the i-th stateful node, or []
% if absent/empty.
%
% jsondecode collapses this nesting inconsistently depending on raggedness:
%   - fully rectangular  -> numeric 4-D [M x R x nRows x 2]
%   - ragged outer       -> cell{M}, each element numeric 3-D [R x nRows x 2]
%                           (or [] for an empty station, or a nested cell)
%   - R == 1 collapses    -> a station element may be a plain [nRows x 2] matrix
% This helper resolves all of these to a uniform M x R cell of 2-D matrices.

C = cell(M, R);
for i = 1:M
    for r = 1:R
        C{i, r} = [];
    end
end
if isempty(raw)
    return;
end

stations = splitOuter(raw, M);
for i = 1:M
    classes = splitOuter(stations{i}, R);
    for r = 1:R
        m = classes{r};
        if ~isempty(m) && ismatrix(m) && size(m, 2) == 2
            C{i, r} = double(m);
        end
    end
end
end

function parts = splitOuter(x, n)
% SPLITOUTER Split the leading logical dimension of a jsondecode'd value into a
% 1 x n cell. Handles cell inputs, numeric N-D arrays (leading dim = list
% index), a single collapsed 2-D matrix (n==1), and empties.
parts = cell(1, n);
if isempty(x)
    return;
end
if iscell(x)
    for i = 1:min(n, numel(x))
        parts{i} = x{i};
    end
    return;
end
if isnumeric(x)
    nd = ndims(x);
    sz = size(x);
    if nd >= 3
        % Leading dimension indexes the list; drop it for each element.
        for i = 1:min(n, sz(1))
            sub = x(i, :, :, :);
            parts{i} = reshape(sub, sz(2:end));
        end
    else % 2-D
        if n == 1
            parts{1} = x;                 % single element is the matrix itself
        elseif size(x, 1) == n
            for i = 1:n
                parts{i} = x(i, :);       % each row is one (vector) element
            end
        else
            parts{1} = x;                 % best effort
        end
    end
end
end
