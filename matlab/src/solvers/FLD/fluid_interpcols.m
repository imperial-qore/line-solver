function b = fluid_interpcols(tg, B, tt)
% B = FLUID_INTERPCOLS(TG, B, TT)
% Clamped piecewise-linear interpolation of the columns of B at scalar time tt.
%
% B is (nrows x ngrid) with column j sampled at time tg(j); tg is a strictly
% increasing grid. The result b is (nrows x 1). Times outside [tg(1),tg(end)]
% are clamped to the boundary columns (zero-order hold outside the grid). This
% is the shared time-varying-input evaluator used both by the TBI transient
% (cross-cell inflow drift) and by the time-varying rate-multiplier closure in
% solver_fluid_odes.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if tt <= tg(1)
    b = B(:,1);
    return
end
if tt >= tg(end)
    b = B(:,end);
    return
end
j = find(tg <= tt, 1, 'last');
wj = (tt - tg(j)) / (tg(j+1) - tg(j));
b = (1-wj)*B(:,j) + wj*B(:,j+1);
end
