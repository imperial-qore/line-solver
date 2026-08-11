function Q=ctmc_makeinfgen(Q)
% Q=CTMC_MAKEINFGEN(Q)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
A=Q-diag(diag(Q)); Q=A-diag(sum(A,2));
end
