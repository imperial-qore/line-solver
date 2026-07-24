function [qnchains,chains] = getChains(self, rt)
% [QNCHAINS,CHAINS] = GETCHAINS(RT)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.getStruct;
chains = sn.chains;
nchains = size(chains,1);
qnchains = cell(1, size(chains,1));
for c=1:nchains
    qnchains{c} = Chain(['Chain',num2str(c)]);
    for r=find(chains(c,:))
        qnchains{c}.addClass(self.classes{r}, sn.visits{c}(:,r), r);
    end
end
end
