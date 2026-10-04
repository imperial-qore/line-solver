function P = spaceCachePendings(s, maxPending)
% P = SPACECACHEPENDINGS(S, MAXPENDING)
%
% Enumerate the ways in which up to MAXPENDING secondary (delayed-hit) requests
% can be distributed over S concurrent in-flight fetches. Returns one row per
% assignment, each row holding S non-negative counts summing to at most
% MAXPENDING. For S = 0 the single empty assignment is returned.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if s == 0
    P = zeros(1,0);
    return
end
if maxPending <= 0
    P = zeros(1,s);
    return
end

P = zeros(1,s); % all fetches with no merged request
for total = 1:maxPending
    P = [P; State.spaceCacheCompositions(total, s)];
end

end
