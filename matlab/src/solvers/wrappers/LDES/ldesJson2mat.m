function m = ldesJson2mat(x, nr, nc)
% M = LDESJSON2MAT(X, NR, NC)
% Convert a jsondecode'd matrix into a numeric NR x NC matrix. Handles numeric
% 2-D arrays, cell arrays (ragged / containing JSON null, which jsondecode maps
% to []), and empties (returns []). JSON null entries become NaN. When NR,NC are
% supplied the result is oriented to NR x NC. Shared by the LDES wrapper methods
% that parse the fully JSON-mediated result of the LDES engine.
if nargin < 2, nr = []; end
if nargin < 3, nc = []; end
if isempty(x)
    m = [];
    return;
end
if iscell(x)
    m = [];
    for i = 1:numel(x)
        ri = x{i};
        if iscell(ri)
            vals = nan(1, numel(ri));
            for j = 1:numel(ri)
                v = ri{j};
                if ~isempty(v) && isnumeric(v)
                    vals(j) = double(v);
                end
            end
        elseif isnumeric(ri)
            vals = reshape(double(ri), 1, []);
        else
            vals = NaN;
        end
        if i == 1
            m = vals;
        elseif size(vals, 2) == size(m, 2)
            m = [m; vals]; %#ok<AGROW>
        else
            m = [];   % ragged beyond repair
            return;
        end
    end
else
    m = double(x);
end
% jsondecode collapses a 1x1 array to a scalar; ensure requested orientation.
if isvector(m) && ~isempty(nr) && ~isempty(nc) && numel(m) == nr * nc
    m = reshape(m, nc, nr).';
end
% Orient to nr x nc if a plain transpose matches.
if ~isempty(nr) && ~isempty(nc) && ~isequal(size(m), [nr nc]) && isequal(size(m), [nc nr])
    m = m.';
end
end
