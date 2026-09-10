function SS = spaceClosedMulti(M, N, caps)
% SS = SPACECLOSEDMULTI(M, N)
% SS = SPACECLOSEDMULTI(M, N, CAPS)
%
% CAPS, when given, is an MxR per (node, class) bound. It is handed straight
% through to spaceClosedSingle, which prunes the per-class enumeration at the
% branch: a (node, class) pair the class never visits has capacity 0 and
% contributes only the zero slot.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

R = length(N);
if nargin < 3
    caps = [];
end
if M == 0
    SS = [];
    return
end
SS = [];
for r = 1:R
    sr = State.spaceClosedSingle(M, N(r), colcap(caps,r));
    if isempty(sr)
        % THE WIDTH MUST SURVIVE AN INFEASIBLE BLOCK. This class's population
        % cannot be placed within its caps, so the whole joint distribution is
        % infeasible -- but the caller CONCATENATES these blocks vertically, and
        % handing it a 0x0 makes that a dimension mismatch rather than a no-op.
        % Return no rows at the FULL width instead.
        SS = zeros(0, M*R);
        return
    end
    if r == 1
        SS = sr;
    else
        SS = State.cartesian(SS, sr);
    end
end
if isempty(SS)
    SS = zeros(0, M*R);
end
end

function c = colcap(caps, r)
if isempty(caps) || size(caps,2) < r
    c = [];
else
    c = caps(:,r)';
end
end
