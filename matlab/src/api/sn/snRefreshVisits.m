function [visits, nodevisits, sn] = snRefreshVisits(sn, chains, rt, rtnodes)
% [VISITS, NODEVISITS, SN] = SNREFRESHVISITS(sn, CHAINS)
%
% Solver traffic equations for the model average number visits to nodes
% and stations.
%
% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

I = sn.nnodes;
M = sn.nstateful;
K = sn.nclasses;
refstat = sn.refstat;
nchains = size(chains,1);

%% obtain chain characteristics
inchain = sn.inchain;
for c=1:nchains
    if sum(refstat(inchain{c}) == refstat(inchain{c}(1))) ~= length(inchain{c})
        refstat(inchain{c}) = refstat(inchain{c}(1));
        %        line_error(mfilename,sprintf('Classes in chain %d have different reference stations. Chain %d classes: %s', c, c, int2str(inchain{c})));
    end
end

%% generate station visits
visits = cell(nchains,1); % visits{c}(i,j) is the number of visits that a chain-c job pays at node i in class j
for c=1:nchains
    cols = zeros(1,M*length(inchain{c}));
    for i=1:M
        nIC = length(inchain{c});
        for ik=1:nIC
            cols(1,(i-1)*nIC+ik) = (i-1)*K+inchain{c}(ik);
        end
    end    
    
    Pchain = rt(cols,cols); % routing probability of the chain
    visited = sum(Pchain,2) > 0;
    %                if ~dtmc_isfeasible(Pchain(visited,visited))
    %                    line_error(mfilename,sprintf('The routing matrix in chain %d is not stochastic. Chain %d classes: %s',c, c, int2str(inchain{c})));
    %                end
    
    alpha_visited = dtmc_solve(Pchain(visited,visited));
    alpha = zeros(1,M*K); alpha(visited) = alpha_visited;
    if max(alpha)>=1-GlobalConstants.FineTol
        %disabled because a self-looping customer is an absorbing chain
        %line_error(mfilename,'One chain has an absorbing state.');
    end
    visits{c} = zeros(M,K);
    for i=1:M
        for k=1:length(inchain{c})
            visits{c}(i,inchain{c}(k)) = alpha((i-1)*length(inchain{c})+k);
        end
    end
    visits{c} = visits{c} / sum(visits{c}(sn.stationToStateful(refstat(inchain{c}(1))),inchain{c}));
    visits{c} = abs(visits{c});
end

%% generate node visits
nodevisits = cell(1,nchains);
for c=1:nchains
    nodes_cols = zeros(1,I*length(inchain{c}));
    for i=1:I
        nIC = length(inchain{c});
        for ik=1:nIC
            nodes_cols(1,(i-1)*nIC+ik) = (i-1)*K+inchain{c}(ik);
        end
    end
    nodes_Pchain = rtnodes(nodes_cols, nodes_cols); % routing probability of the chain
    nodes_visited = sum(nodes_Pchain,2) > 0;

    nodes_alpha_visited = dtmc_solve(nodes_Pchain(nodes_visited,nodes_visited));
    nodes_alpha = zeros(1,I); nodes_alpha(nodes_visited) = nodes_alpha_visited;
    nodevisits{c} = zeros(I,K);
    for i=1:I
        for k=1:length(inchain{c})
            nodevisits{c}(i,inchain{c}(k)) = nodes_alpha((i-1)*length(inchain{c})+k);
        end
    end
    nodevisits{c} = nodevisits{c} / sum(nodevisits{c}(sn.statefulToNode(refstat(inchain{c}(1))),inchain{c}));
    nodevisits{c}(nodevisits{c}<0) = 0; % remove small numerical perturbations
end

for c=1:nchains
    nodevisits{c}(isnan(nodevisits{c})) = 0;
end
%% save results in sn
sn.visits = visits;
sn.nodevisits = nodevisits;
sn.inchain = inchain;
end