function C = spaceCacheCompositions(total, parts)
% C = SPACECACHECOMPOSITIONS(TOTAL, PARTS)
%
% Enumerate the weak compositions of TOTAL into exactly PARTS non-negative
% integers, one composition per row.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if parts == 1
    C = total;
    return
end
C = [];
for first = 0:total
    tail = State.spaceCacheCompositions(total - first, parts - 1);
    C = [C; [repmat(first, size(tail,1), 1), tail]];
end

end
