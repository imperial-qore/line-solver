function space = pollingSpace(sn, ind, space)
% SPACE = POLLINGSPACE(SN, IND, SPACE)
%
% Append the polling controller columns to the rows of SPACE, which must hold
% the [buffer, server, routing-variable] layout of a polling station. Rows are
% expanded into one row per controller configuration the discipline can
% occupy, and rows for which no configuration exists are dropped.
%
% Enumerating the controller per row rather than as a blind cartesian product
% is what keeps the state space tight and the chain irreducible: pos is pinned
% to the class in service, a switchover excludes a busy service facility, and a
% park excludes a non-empty station. A cartesian product would instead admit
% states such as "serving buffer 1 while a class-2 job occupies the server",
% which no transition can reach or leave in a way consistent with the marginals.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

pinfo = State.pollingInfo(sn, ind);
if isempty(pinfo) || isempty(space)
    return
end

R = sn.nclasses;
ist = sn.nodeToStation(ind);
K = zeros(1,R);
for r=1:R
    if isempty(sn.proc{ist}{r})
        K(r) = 1;
    else
        K(r) = length(sn.proc{ist}{r}{1});
    end
end
nbufcols = 1:R;
nsrvcols = R + (1:sum(K));

out = [];
for row=1:size(space,1)
    nbuf = space(row, nbufcols);
    srv = space(row, nsrvcols);
    % class occupying the single service facility, 0 when it is empty
    srvclass = 0;
    for r=1:R
        if sum(srv((sum(K(1:r-1))+1):sum(K(1:r)))) > 0
            srvclass = r;
            break
        end
    end
    blocks = State.pollingProject(pinfo, State.pollingBlocks(pinfo, srvclass, nbuf, R));
    for b=1:size(blocks,1)
        out(end+1,:) = [space(row,:), blocks(b,:)]; %#ok<AGROW>
    end
end
if isempty(out)
    space = zeros(0, size(space,2) + pinfo.width);
else
    space = out;
end
end
