function [qnchains,chains] = getChains(self, rt)
% [QNCHAINS,CHAINS] = GETCHAINS(RT)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

if nargin<2 %~exist('rt','var')
    rt = getRoutingMatrix(self);
end

sn = self.getStruct;
M = sn.nstations;
K = sn.nclasses;
chains = sn.chains;
refstat = sn.refstat;
nchains = size(chains,1);
inchain = cell(1,nchains);
for c=1:nchains
    inchain{c} = find(chains(c,:));
    if sum(refstat(inchain{c}) == refstat(inchain{c}(1))) ~= length(inchain{c})
        refstat(inchain{c}) = refstat(inchain{c}(1));
        %        line_error(mfilename,sprintf('Classes in chain %d have different reference stations. Chain %d classes: %s', c, c, int2str(inchain{c})));
    end
end

visits = cell(size(chains,1),1); % visits{c}(i,j) is the number of visits that a chain-c job pays at node i in class j
for c=1:size(chains,1)
    inchain{c} = find(chains(c,:));
    cols = [];
    for i=1:M
        for k=inchain{c}(:)'
            cols(end+1) = (i-1)*K+k;
        end
    end
    Pchain = rt(cols,cols); % routing probability of the chain
    visited = sum(Pchain,2,"omitnan") > 0;
    
    %                Pchain(visited,visited)
    %                if ~dtmc_isfeasible(Pchain(visited,visited))
    %                    line_error(mfilename,sprintf('The routing matrix in chain %d is not stochastic. Chain %d classes: %s',c, c, int2str(inchain{c})));
    %                end    
    alpha_visited = dtmc_solve(Pchain(visited,visited));
    alpha = zeros(1,M*K); alpha(visited) = alpha_visited;    
    if max(alpha)>=1-GlobalConstants.FineTol
        line_error(mfilename,'One chain has an absorbing state.');
    end
    
    visits{c} = zeros(M,K);
    for i=1:M
        for k=1:length(inchain{c})
            visits{c}(i,inchain{c}(k)) = alpha((i-1)*length(inchain{c})+k);
        end
    end
    visits{c} = visits{c} / sum(visits{c}(refstat(inchain{c}(1)),inchain{c}));
end
qnchains = cell(1, size(chains,1));
for c=1:size(chains,1)
    qnchains{c} = Chain(['Chain',num2str(c)]);
    for r=find(chains(c,:))
        qnchains{c}.addClass(self.classes{r}, visits{c}(:,r), r);
    end
end

end
