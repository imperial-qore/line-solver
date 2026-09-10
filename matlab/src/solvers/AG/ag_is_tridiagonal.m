function result = ag_is_tridiagonal(Q)
% IS_TRIDIAGONAL Check if a matrix is tridiagonal
n = size(Q, 1);
result = true;
for i = 1:n
    for j = 1:n
        if abs(i - j) > 1 && abs(Q(i, j)) > 1e-14
            result = false;
            return;
        end
    end
end
end
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
