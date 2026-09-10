function tags = tagPad(tags, ntot, R)
% TAGS = TAGPAD(TAGS, NTOT, R)
%
% Bring a START/PREEMPT tag matrix up to NTOT rows, one per successor, so that
% it can be indexed with the same row indices as outspace. Arcs that carry no
% tag are zero rows, which is why the tag sites only ever write the rows they
% touch (see State.tagArc).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isempty(tags)
    tags = zeros(ntot, R);
    return
end
if size(tags,2) ~= R
    tags = [tags, zeros(size(tags,1), R-size(tags,2))];
end
if size(tags,1) < ntot
    tags = [tags; zeros(ntot-size(tags,1), R)];
elseif size(tags,1) > ntot
    tags = tags(1:ntot,:);
end
end
