function [m1, scv, m2] = lqn_ph_moments(alpha, T)
% [M1, SCV, M2] = LQN_PH_MOMENTS(ALPHA, T)
%
% First two moments of the phase-type law (ALPHA,T) without building a
% Distribution object, which is what the layered fixed point needs at every
% iteration for every composed entry law. A defective ALPHA carries an atom at
% zero and contributes nothing to either moment.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

alpha = reshape(alpha, 1, []);
e = ones(size(T,1), 1);
x1 = -(T \ e);
m1 = alpha * x1;
x2 = -(T \ x1);
m2 = 2 * (alpha * x2);
if ~isfinite(m1) || m1 <= GlobalConstants.FineTol
    m1 = GlobalConstants.FineTol;
    m2 = 2 * m1^2;
    scv = 1.0;
    return
end
scv = m2 / m1^2 - 1;
if ~isfinite(scv) || scv <= GlobalConstants.FineTol
    scv = GlobalConstants.FineTol;
end
end
