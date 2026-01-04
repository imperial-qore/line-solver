function hashid = getHash(sn, ind, inspace)
% HASHID = GETHASH(QN, IND, INSPACE)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isempty(inspace)
    hashid = -1;
    return
end

% ind: node index
%ist = sn.nodeToStation(ind);
isf = sn.nodeToStateful(ind);

inspace = [zeros(size(inspace,1),size(sn.space{isf},2)-size(inspace,2)), inspace];
if isempty(sn.space{isf})
    line_error(mfilename,'Station state space is not initialized. Use setStateSpace method.\n');
end

hashid=zeros(size(inspace,1),1);
for j=1:size(inspace,1)
    if size(sn.space{isf},2)<size(inspace(j,:),2) % cutoff exceeded
        hashid(j,1) = -1 ;
    else
        hashid(j,1) = matchrow(sn.space{isf},inspace(j,:));
    end
end
end
