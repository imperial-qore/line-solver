function [hashid, sn] = getHashOrAdd(sn, ind, inspace)
% [HASHID, QN] = GETHASHORADD(QN, IND, INSPACE)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isempty(inspace)
    hashid = -1;
    return
end

% ind: node index
%ist = sn.nodeToStation(ind);
isf = sn.nodeToStateful(ind);

if isempty(sn.space{isf})
    line_error(mfilename,'Station state space is not initialized. Use setStateSpace method.\n');
end

% resize
if size(inspace,2) < size(sn.space{isf},2)
    inspace = [zeros(size(inspace,1),size(sn.space{isf},2)-size(inspace,2)), inspace];
elseif size(inspace,2) > size(sn.space{isf},2)
    sn.space{isf} = [zeros(size(sn.space{isf},1),size(inspace,2)-size(sn.space{isf},2)),sn.space{isf}];
end

hashid = matchrows(sn.space{isf},inspace);
%hashid=zeros(size(inspace,1),1);
for j=1:size(inspace,1)
    %    hashid(j,1) = matchrow(sn.space{isf},inspace(j,:));
    if hashid(j,1) <0
        sn.space{isf}(end+1,:) = inspace(j,:);
        hashid(j,1) = size(sn.space{isf},1);
    end
end
end
