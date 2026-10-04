function idx = ag_blk(n, m)
% Row/column range of level N (0-based) in a component with M phases per level.
idx = (n * m + 1):(n * m + m);
end
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
