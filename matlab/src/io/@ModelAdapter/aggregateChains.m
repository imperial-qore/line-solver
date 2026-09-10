function [chainModel, alpha, deaggInfo] = aggregateChains(model, suffix)
% AGGREGATECHAINS Transform a multi-class model into an equivalent chain-aggregated model
%
% [chainModel, alpha, deaggInfo] = AGGREGATECHAINS(model)
% [chainModel, alpha, deaggInfo] = AGGREGATECHAINS(model, suffix)
%
% This function transforms a queueing network model with multiple classes
% into a stochastically equivalent model where each chain becomes a single
% class. Classes belonging to the same chain (i.e., classes that can switch
% into each other) are merged into one aggregate class.
%
% The aggregated model preserves:
% - Total chain population (closed chains)
% - Total arrival rate (open chains)
% - Service demands at each station
% - Routing structure at the chain level
%
% Parameters:
%   model  - Source Network model with potentially multiple classes per chain
%   suffix - Optional suffix for chain class names (default: '')
%
% Returns:
%   chainModel - New Network model with one class per chain
%   alpha      - Aggregation factors matrix (M x K) from original model
%   deaggInfo  - Struct containing all information needed to deaggregate
%                chain-level results back to class-level results
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2 || isempty(suffix)
    suffix = '';
end

% Get network structure
sn = model.getStruct();

% Extract dimensions
M = sn.nstations;
K = sn.nclasses;
C = sn.nchains;

% If each class is its own chain, just return a copy
if C == K
    chainModel = model.copy();
    alpha = eye(M, K);
    deaggInfo = struct();
    deaggInfo.alpha = alpha;
    deaggInfo.inchain = sn.inchain;
    deaggInfo.originalSn = sn;
    deaggInfo.isAggregated = false;
    return;
end

% Get aggregation parameters from existing API
[Lchain, STchain, Vchain, alpha, Nchain, SCVchain, refstatchain] = sn_get_demands_chain(sn);

% Determine which chains are open vs closed
isOpenChain = false(1, C);
lambdaChain = zeros(1, C);
sourceIdx = find(sn.nodetype == NodeType.Source);
if ~isempty(sourceIdx)
    sourceStationIdx = sn.nodeToStation(sourceIdx);
end

for c = 1:C
    inchain = sn.inchain{c};
    isOpenChain(c) = any(isinf(sn.njobs(inchain)));
    if isOpenChain(c) && ~isempty(sourceIdx)
        % Sum arrival rates for all classes in this chain
        lambdaChain(c) = sum(sn.rates(sourceStationIdx, inchain), 'omitnan');
    end
end

% Create the new aggregated network
chainModel = Network(sprintf('%s_aggregated', model.getName()));

% Map from original node index to new node
nodeMap = cell(sn.nnodes, 1);
stationMap = cell(M, 1);

% First, identify which nodes to copy (skip auto-added ClassSwitch nodes)
nodesToCopy = [];
for i = 1:length(model.nodes)
    node = model.nodes{i};
    % Skip auto-added class switch nodes
    if isa(node, 'ClassSwitch') && node.autoAdded
        continue;
    end
    nodesToCopy(end+1) = i;
end

% Create nodes in the new model
newNodeIdx = 0;
for idx = 1:length(nodesToCopy)
    i = nodesToCopy(idx);
    node = model.nodes{i};
    newNodeIdx = newNodeIdx + 1;

    switch class(node)
        case 'Source'
            nodeMap{i} = Source(chainModel, node.name);
        case 'Sink'
            nodeMap{i} = Sink(chainModel, node.name);
        case 'Queue'
            nodeMap{i} = Queue(chainModel, node.name, node.schedStrategy);
            if ~isinf(node.numberOfServers)
                nodeMap{i}.setNumberOfServers(node.numberOfServers);
            end
            if ~isempty(node.cap) && isfinite(node.cap)
                nodeMap{i}.setCapacity(node.cap);
            end
        case 'Delay'
            nodeMap{i} = Delay(chainModel, node.name);
        case 'Router'
            nodeMap{i} = Router(chainModel, node.name);
        case 'ClassSwitch'
            % User-defined ClassSwitch nodes should not exist in the
            % aggregated model since class switching is eliminated
            continue;
        otherwise
            line_warning(mfilename, sprintf('Node type %s not fully supported in chain aggregation.', class(node)));
            continue;
    end

    % Store station mapping
    if isa(nodeMap{i}, 'Station')
        stationIdx = sn.nodeToStation(i);
        if stationIdx > 0
            stationMap{stationIdx} = nodeMap{i};
        end
    end
end

% Create chain classes
chainClass = cell(1, C);
for c = 1:C
    inchain = sn.inchain{c};

    % Build chain class name from original class names
    classNames = sn.classnames(inchain);
    if length(classNames) == 1
        chainClassName = char(classNames{1});
    else
        chainClassName = sprintf('Chain%d', c);
    end
    if ~isempty(suffix)
        chainClassName = [chainClassName, suffix];
    end

    if isOpenChain(c)
        % Open chain
        chainClass{c} = OpenClass(chainModel, chainClassName);
    else
        % Closed chain - find reference station
        refStationIdx = refstatchain(c);
        refStation = stationMap{refStationIdx};
        if isempty(refStation)
            line_error(mfilename, sprintf('Reference station %d for chain %d not found in aggregated model.', refStationIdx, c));
        end
        chainClass{c} = ClosedClass(chainModel, chainClassName, Nchain(c), refStation);
    end
end

% Set arrival rates for open chains
if any(isOpenChain)
    chainSource = chainModel.getSource();
    for c = 1:C
        if isOpenChain(c) && lambdaChain(c) > 0
            chainSource.setArrival(chainClass{c}, Exp(lambdaChain(c)));
        end
    end
end

% Set service times at each station
for i = 1:M
    station = stationMap{i};
    if isempty(station) || isa(station, 'Source') || isa(station, 'Sink')
        continue;
    end

    for c = 1:C
        if STchain(i, c) > 0 && isfinite(STchain(i, c))
            scv = SCVchain(i, c);
            if isnan(scv) || scv <= 0
                scv = 1; % Default to exponential
            end

            % Choose distribution based on SCV
            if abs(scv - 1) < GlobalConstants.FineTol
                % Exponential (SCV = 1)
                dist = Exp(1/STchain(i, c));
            elseif scv < 1
                % SCV < 1: use Erlang or deterministic
                if scv < GlobalConstants.FineTol
                    dist = Det(STchain(i, c));
                else
                    % Erlang: SCV = 1/k, so k = 1/SCV
                    k = round(1/scv);
                    if k < 1
                        k = 1;
                    end
                    dist = Erlang.fitMeanAndOrder(STchain(i, c), k);
                end
            else
                % SCV > 1: use HyperExp or Cox2
                dist = HyperExp.fitMeanAndSCV(STchain(i, c), scv);
            end

            station.setService(chainClass{c}, dist);
        else
            % No service for this chain at this station
            station.setService(chainClass{c}, Disabled.getInstance());
        end
    end
end

% Build routing matrix from aggregated routing probabilities
I_new = length(chainModel.nodes);
P = chainModel.initRoutingMatrix();

% Create mapping from original station to new model node index
stationToNewNode = zeros(M, 1);
for i = 1:M
    iNode = sn.stationToNode(i);
    if ~isempty(nodeMap{iNode})
        % Find the index of nodeMap{iNode} in chainModel.nodes
        for n = 1:I_new
            if chainModel.nodes{n} == nodeMap{iNode}
                stationToNewNode(i) = n;
                break;
            end
        end
    end
end

% For each chain, compute routing probabilities between stations
for c = 1:C
    inchain = sn.inchain{c};

    % Build routing for this chain
    for i = 1:M
        iNode = sn.stationToNode(i);
        if isempty(nodeMap{iNode}) || stationToNewNode(i) == 0
            continue;
        end

        % Skip source/sink in routing calculations for closed chains
        if ~isOpenChain(c) && (sn.nodetype(iNode) == NodeType.Source || sn.nodetype(iNode) == NodeType.Sink)
            continue;
        end

        % Get stateful index for station i (sn.rt is indexed by stateful nodes)
        isf_i = sn.stationToStateful(i);

        for j = 1:M
            jNode = sn.stationToNode(j);
            if isempty(nodeMap{jNode}) || stationToNewNode(j) == 0
                continue;
            end

            % Get stateful index for station j
            isf_j = sn.stationToStateful(j);

            % Compute aggregated routing probability from station i to station j
            pij = 0;
            for k = inchain
                for s = inchain
                    % Use stateful indices to access sn.rt
                    fromIdx = (isf_i-1)*K + k;
                    toIdx = (isf_j-1)*K + s;
                    if fromIdx <= size(sn.rt, 1) && toIdx <= size(sn.rt, 2)
                        p_ks = sn.rt(fromIdx, toIdx);
                        if alpha(i, k) > 0 && p_ks > 0
                            pij = pij + alpha(i, k) * p_ks;
                        end
                    end
                end
            end

            if pij > GlobalConstants.FineTol
                iNodeNew = stationToNewNode(i);
                jNodeNew = stationToNewNode(j);
                P{c, c}(iNodeNew, jNodeNew) = pij;
            end
        end
    end

    % Normalize routing probabilities (should be close to 1 already if indexing is correct)
    for iNew = 1:I_new
        rowSum = sum(P{c, c}(iNew, :));
        if rowSum > GlobalConstants.FineTol
            if abs(rowSum - 1) > 0.01
                line_warning(mfilename, 'Large normalization correction at node %d for chain %d: rowSum=%.4f', iNew, c, rowSum);
            end
            P{c, c}(iNew, :) = P{c, c}(iNew, :) / rowSum;
        end
    end
end

% Link the model with the routing matrix
chainModel.link(P);

% Build deaggregation information
deaggInfo = struct();
deaggInfo.alpha = alpha;
deaggInfo.Vchain = Vchain;
deaggInfo.STchain = STchain;
deaggInfo.Lchain = Lchain;
deaggInfo.SCVchain = SCVchain;
deaggInfo.Nchain = Nchain;
deaggInfo.lambdaChain = lambdaChain;
deaggInfo.isOpenChain = isOpenChain;
deaggInfo.inchain = sn.inchain;
deaggInfo.refstat = sn.refstat;
deaggInfo.refstatchain = refstatchain;
deaggInfo.originalSn = sn;
deaggInfo.isAggregated = true;
deaggInfo.nclasses = K;
deaggInfo.nchains = C;

end
