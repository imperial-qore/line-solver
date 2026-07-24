function SS = spaceClosedMultiCS(M, N, chains)
% SS = SPACECLOSEDMULTICS(M, N, CHAINS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% State space for closed multiclass CQN with class-switching

C = size(chains,1);
chainInitPos = cell(1,C);
for c=1:C
    inchain=find(chains(c,:));
    chainInitPos{c} = multichoose(length(inchain),sum(N(inchain)));
end
SS = [];
chainInitPosLen = cellfun(@(c) size(c,1),chainInitPos)-1;
v = pprod(chainInitPosLen);
while v>=0
    subN = [];
    for c=1:C
        subN((end+1):(end+size(chainInitPos{c},2))) = chainInitPos{c}(v(c)+1,:);
    end
    SS = [SS; State.spaceClosedMulti(M, subN)];
    v = pprod(v,chainInitPosLen);
end
end
