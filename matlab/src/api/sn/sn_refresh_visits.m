%{ @file sn_refresh_visits.m
 %  @brief Solves traffic equations to compute visit ratios
 %
 %  @author LINE Development Team
 %  Copyright (c) 2012-2026, Imperial College London
 %  All rights reserved.
%}

%{
 % @brief Solves traffic equations to compute visit ratios
 %
 % @details
 % This function solves the traffic equations to compute the average number
 % of visits to nodes and stations for each chain in the network.
 %
 % @par Syntax:
 % @code
 % [visits, nodevisits, sn] = sn_refresh_visits(sn, chains, rt, rtnodes)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>chains<td>Chain definitions
 % <tr><td>rt<td>Station routing matrix
 % <tr><td>rtnodes<td>Node routing matrix
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>visits<td>Cell array of visit ratios at stations per chain
 % <tr><td>nodevisits<td>Cell array of visit ratios at nodes per chain
 % <tr><td>sn<td>Updated network structure with visit fields populated
 % </table>
%}
function [visits, nodevisits, sn] = sn_refresh_visits(sn, chains, rt, rtnodes)

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
    for ist=1:M
        nIC = length(inchain{c});
        for ik=1:nIC
            cols(1,(ist-1)*nIC+ik) = (ist-1)*K+inchain{c}(ik);
        end
    end    
    
    Pchain = rt(cols,cols); % routing probability of the chain

    % Handle NaN values in routing matrix (e.g., from Cache class switching)
    % For visits calculation, replace NaN with equal probabilities
    for row = 1:size(Pchain,1)
        nan_cols = isnan(Pchain(row,:));
        if any(nan_cols)
            % Get the non-NaN sum for this row
            non_nan_sum = sum(Pchain(row, ~nan_cols));
            % Distribute remaining probability equally among NaN entries
            remaining_prob = max(0, 1 - non_nan_sum);
            n_nan = sum(nan_cols);
            if n_nan > 0 && remaining_prob > 0
                Pchain(row, nan_cols) = remaining_prob / n_nan;
            else
                Pchain(row, nan_cols) = 0;
            end
        end
    end

    visited = sum(Pchain,2) > 0;

    % Normalize routing matrix for Fork-containing models
    % Fork nodes have row sums > 1 (sending to all branches with prob 1 each)
    % which causes dtmc_solve_reducible to fail. Normalize to make stochastic.
    % Record original row sums to correct visit ratios after DTMC solve.
    row_sums = ones(size(Pchain,1), 1);
    if any(sn.nodetype == NodeType.Fork)
        for row = 1:size(Pchain,1)
            rs = sum(Pchain(row,:));
            row_sums(row) = rs;
            if rs > GlobalConstants.FineTol
                Pchain(row,:) = Pchain(row,:) / rs;
            end
        end
    end

    % Use dtmc_solve as primary, fallback to dtmc_solve_reducible for chains with transient states
    Pchain_visited = Pchain(visited,visited);
    alpha_visited = dtmc_solve(Pchain_visited);
    % Fallback to dtmc_solve_reducible if dtmc_solve fails (e.g., reducible chain)
    if all(alpha_visited == 0) || any(isnan(alpha_visited))
        [alpha_visited, ~, ~, ~, ~] = dtmc_solve_reducible(Pchain_visited, [], struct('tol', GlobalConstants.FineTol));
    end
    alpha = zeros(1,M*K); alpha(visited) = alpha_visited;
    if max(alpha)>=1-GlobalConstants.FineTol
        %disabled because a self-looping customer is an absorbing chain
        %line_error(mfilename,'One chain has an absorbing state.');
    end

    % Apply Fork fanout correction: multiply fork-scope station visits by fanout
    % This corrects for the normalization of rows with sum > 1 (caused by Fork
    % routing being collapsed into the stateful routing matrix).
    % Block propagation through Join nodes to prevent the correction from
    % leaking back through cycles (which would cancel out during normalization).
    if any(sn.nodetype == NodeType.Fork)
        n = size(Pchain, 1);
        nIC = length(inchain{c});
        % Build adjacency but block outgoing edges from Join nodes
        adj = Pchain > GlobalConstants.FineTol;
        for isf = 1:M
            nd = sn.statefulToNode(isf);
            if sn.nodetype(nd) == NodeType.Join
                for k = 1:nIC
                    adj((isf-1)*nIC + k, :) = false;
                end
            end
        end
        % Compute transitive closure with join-blocked adjacency
        reachable = adj;
        for iter = 1:ceil(log2(n))
            reachable = reachable | (reachable * reachable > 0);
        end

        % For each fork row, multiply only fork-scope reachable states by fanout
        for row = 1:n
            fanout = row_sums(row);
            if fanout > 1 + GlobalConstants.FineTol
                for col = 1:n
                    if reachable(row, col)
                        alpha(col) = alpha(col) * fanout;
                    end
                end
            end
        end
    end

    visits{c} = zeros(M,K);
    for ist=1:M
        for k=1:length(inchain{c})
            visits{c}(ist,inchain{c}(k)) = alpha((ist-1)*length(inchain{c})+k);
        end
    end
    normSum = sum(visits{c}(sn.stationToStateful(refstat(inchain{c}(1))),inchain{c}));
    if normSum > GlobalConstants.FineTol
        visits{c} = visits{c} / normSum;
    end
    visits{c} = abs(visits{c});
end

%% generate node visits
nodevisits = cell(1,nchains);
for c=1:nchains
    nodes_cols = zeros(1,I*length(inchain{c}));
    for ind=1:I
        nIC = length(inchain{c});
        for ik=1:nIC
            nodes_cols(1,(ind-1)*nIC+ik) = (ind-1)*K+inchain{c}(ik);
        end
    end
    nodes_Pchain = rtnodes(nodes_cols, nodes_cols); % routing probability of the chain

    % Handle NaN values in routing matrix (e.g., from Cache class switching)
    % For visits calculation, replace NaN with equal probabilities
    for row = 1:size(nodes_Pchain,1)
        nan_cols = isnan(nodes_Pchain(row,:));
        if any(nan_cols)
            % Get the non-NaN sum for this row
            non_nan_sum = sum(nodes_Pchain(row, ~nan_cols));
            % Distribute remaining probability equally among NaN entries
            remaining_prob = max(0, 1 - non_nan_sum);
            n_nan = sum(nan_cols);
            if n_nan > 0 && remaining_prob > 0
                nodes_Pchain(row, nan_cols) = remaining_prob / n_nan;
            else
                nodes_Pchain(row, nan_cols) = 0;
            end
        end
    end

    nodes_visited = sum(nodes_Pchain,2) > 0;

    % Normalize routing matrix for Fork-containing models
    % Record original row sums to correct visit ratios after DTMC solve.
    nodes_row_sums = ones(size(nodes_Pchain,1), 1);
    if any(sn.nodetype == NodeType.Fork)
        for row = 1:size(nodes_Pchain,1)
            rs = sum(nodes_Pchain(row,:));
            nodes_row_sums(row) = rs;
            if rs > GlobalConstants.FineTol
                nodes_Pchain(row,:) = nodes_Pchain(row,:) / rs;
            end
        end
    end

    % Use dtmc_solve as primary, fallback to dtmc_solve_reducible for chains with transient states
    nodes_Pchain_visited = nodes_Pchain(nodes_visited,nodes_visited);
    nodes_alpha_visited = dtmc_solve(nodes_Pchain_visited);
    % Fallback to dtmc_solve_reducible if dtmc_solve_reducible fails (e.g., reducible chain)
    if all(nodes_alpha_visited == 0) || any(isnan(nodes_alpha_visited))
        [nodes_alpha_visited, ~, ~, ~, ~] = dtmc_solve_reducible(nodes_Pchain_visited, [], struct('tol', GlobalConstants.FineTol));
    end
    nodes_alpha = zeros(1,I*K); nodes_alpha(nodes_visited) = nodes_alpha_visited;

    % Apply Fork fanout correction for node visits
    % Block propagation through Join nodes to prevent the correction from
    % leaking back through cycles.
    if any(sn.nodetype == NodeType.Fork)
        n_nodes = size(nodes_Pchain, 1);
        nIC = length(inchain{c});
        % Build adjacency but block outgoing edges from Join nodes
        nodes_adj = nodes_Pchain > GlobalConstants.FineTol;
        for nd = 1:I
            if sn.nodetype(nd) == NodeType.Join
                for k = 1:nIC
                    nodes_adj((nd-1)*nIC + k, :) = false;
                end
            end
        end
        % Compute transitive closure with join-blocked adjacency
        nodes_reachable = nodes_adj;
        for iter = 1:ceil(log2(n_nodes))
            nodes_reachable = nodes_reachable | (nodes_reachable * nodes_reachable > 0);
        end

        % For each fork row, multiply only fork-scope reachable states by fanout
        for row = 1:n_nodes
            fanout = nodes_row_sums(row);
            if fanout > 1 + GlobalConstants.FineTol
                for col = 1:n_nodes
                    if nodes_reachable(row, col)
                        nodes_alpha(col) = nodes_alpha(col) * fanout;
                    end
                end
            end
        end
    end

    nodevisits{c} = zeros(I,K);
    for ind=1:I
        for k=1:length(inchain{c})
            nodevisits{c}(ind,inchain{c}(k)) = nodes_alpha((ind-1)*length(inchain{c})+k);
        end
    end
    nodeNormSum = sum(nodevisits{c}(sn.statefulToNode(refstat(inchain{c}(1))),inchain{c}));
    if nodeNormSum > GlobalConstants.FineTol
        nodevisits{c} = nodevisits{c} / nodeNormSum;
    end
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