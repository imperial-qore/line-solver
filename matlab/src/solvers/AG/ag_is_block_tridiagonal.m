function result = ag_is_block_tridiagonal(Q, lvl)
% IS_BLOCK_TRIDIAGONAL True when every transition of Q stays within the
% neighbouring level, LVL being the level index of each state.
n = size(Q, 1);
result = true;
for i = 1:n
    for j = 1:n
        if abs(lvl(i) - lvl(j)) > 1 && abs(Q(i, j)) > 1e-14
            result = false;
            return;
        end
    end
end
end
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
