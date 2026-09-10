function tags = tagArc(tags, ntot, nblk, R, cls)
% TAGS = TAGARC(TAGS, NTOT, NBLK, R, CLS)
%
% Record a START or PREEMPT tag on the block of NBLK successor rows that ends
% at row NTOT of outspace, i.e. on the rows just appended by the caller. CLS
% is the class the tag refers to: a scalar shared by the whole block, or one
% entry per row of the block, with 0 meaning "no tag on this row".
%
% Rows appended earlier and left untagged are padded with zeros here rather
% than at their own append site, so only the arcs that actually carry a tag
% need to say anything. The caller pads the tail once, after the last append,
% with State.tagPad.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isempty(tags)
    tags = zeros(0,R);
end
if nblk <= 0
    return
end
nbefore = ntot - nblk;
if size(tags,1) < nbefore
    tags = [tags; zeros(nbefore-size(tags,1), R)];
end
blk = zeros(nblk, R);
if isscalar(cls)
    if cls > 0
        blk(:,cls) = 1;
    end
else
    cls = cls(:);
    if numel(cls) ~= nblk
        line_error(mfilename, sprintf(['A tag class list of %d entries cannot annotate a block of %d successor rows: ' ...
            'pass one class per row, or a single class for the whole block.'], numel(cls), nblk));
    end
    for j = 1:nblk
        if cls(j) > 0
            blk(j,cls(j)) = blk(j,cls(j)) + 1;
        end
    end
end
tags = [tags; blk];
end
