function cells = tbi_partition(sn, options)
% CELLS = TBI_PARTITION(SN, OPTIONS)
%
% Partition the station set into cells for trajectory-based iteration.
% Honors options.config.tbi_cells, a cell array of disjoint station index
% vectors covering 1:nstations. Otherwise stations are agglomerated
% greedily on the symmetrized station-level routing weights, targeting
% options.config.tbi_cellsize stations per cell (default 5), so that
% strongly coupled stations share a cell and connecting flows stay weak.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;
K = sn.nclasses;

if isfield(options.config,'tbi_cells') && ~isempty(options.config.tbi_cells)
    cells = options.config.tbi_cells;
    covered = sort(cell2mat(cellfun(@(x) x(:)', cells(:)', 'UniformOutput', false)));
    if numel(covered) ~= M || any(covered(:)' ~= 1:M)
        line_error(mfilename,sprintf('options.config.tbi_cells must be a partition of the station set 1:%d.', M));
    end
    return
end

if isfield(options.config,'tbi_cellsize') && ~isempty(options.config.tbi_cellsize)
    cellsize = options.config.tbi_cellsize;
else
    cellsize = 5;
end

% station-level coupling weights, aggregated over classes and symmetrized
rt = sn.rt;
W = zeros(M);
for i = 1:M
    for j = 1:M
        W(i,j) = sum(sum(rt((i-1)*K+(1:K), (j-1)*K+(1:K))));
    end
end
A = W + W';
A(1:M+1:end) = 0;

ncells_target = max(1, ceil(M/cellsize));
cells = num2cell(1:M);
C = A; % inter-cell coupling weights
while numel(cells) > ncells_target
    n = numel(cells);
    szs = cellfun(@numel, cells);
    % most coupled cell pair, respecting the cell size cap
    maxc = -1; besta = 0; bestb = 0;
    for a = 1:n
        for b = a+1:n
            if szs(a)+szs(b) <= 2*cellsize && C(a,b) > maxc
                maxc = C(a,b);
                besta = a;
                bestb = b;
            end
        end
    end
    if besta == 0 % every merge exceeds the size cap: merge the two smallest
        [~,ord] = sort(szs);
        besta = min(ord(1),ord(2));
        bestb = max(ord(1),ord(2));
    end
    cells{besta} = [cells{besta}, cells{bestb}];
    C(besta,:) = C(besta,:) + C(bestb,:);
    C(:,besta) = C(:,besta) + C(:,bestb);
    C(besta,besta) = 0;
    cells(bestb) = [];
    C(bestb,:) = [];
    C(:,bestb) = [];
end
end
