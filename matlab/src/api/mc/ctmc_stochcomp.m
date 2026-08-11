function [S,Q11,Q12,Q21,Q22,T]=ctmc_stochcomp(Q,I)
% [S,Q11,Q12,Q21,Q22,T] = CTMC_STOCHCOMP(Q,I)
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if nargin==1
    I=1:ceil(length(Q)/2);
end
isSelected = false(1, length(Q));
for idx = 1:length(I)
    isSelected(I(idx)) = true;
end
Ic = zeros(1, length(Q) - length(I));
icPos = 0;
for idx = 1:length(Q)
    if ~isSelected(idx)
        icPos = icPos + 1;
        Ic(icPos) = idx;
    end
end
Q11 = Q(I,I);
Q12 = Q(I,Ic);
Q21 = Q(Ic,I);
Q22 = Q(Ic,Ic);
% see _kb/03-api-layer.md (mc/ additions) for the ILUT/fill-in rationale
GMRES_MIN_STATES = 6000;
T = [];
if size(Q22,1) > GMRES_MIN_STATES
    [T,gflag] = ctmc_gmres_multi(-Q22, Q21);
    if gflag ~= 0
        T = [];
    end
end
if isempty(T)
    T = (-Q22) \ Q21;
    if ~all(isfinite(T(:)))
        % Backslash returns NaN on a rank-deficient complement, which then
        % contaminates every downstream metric. Fall back to the minimum-norm
        % least-squares solution, as the Java and Python kernels do.
        T = lsqminnorm(full(-Q22), full(Q21));
    end
end
T = Q12*T;
S = Q11+T;
end
