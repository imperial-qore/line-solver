function [fesModel, fesStation, deaggInfo] = aggregateFES(model, stationSubset, options)
% AGGREGATEFES Replace a station subset with a Flow-Equivalent Server (FES)
%
% [fesModel, fesStation, deaggInfo] = AGGREGATEFES(model, stationSubset)
% [fesModel, fesStation, deaggInfo] = AGGREGATEFES(model, stationSubset, options)
%
% This function replaces a subset of stations in a closed product-form
% queueing network with a single Flow-Equivalent Server (FES). The FES has
% Limited Joint Dependence (LJD) service rates where the rate for class-c
% in state (n1,...,nK) equals the throughput of class-c in an isolated
% subnetwork consisting only of the subset stations.
%
% Parameters:
%   model         - Closed product-form Network model
%   stationSubset - Cell array of Station objects to aggregate
%   options       - Optional struct with fields:
%       .solver   - Solver for throughput computation ('mva' default)
%       .cutoffs  - Per-class population cutoffs (default: njobs per class)
%       .verbose  - Enable verbose output (default: false)
%
% Returns:
%   fesModel   - New Network with FES replacing the subset
%   fesStation - Reference to the FES Queue station
%   deaggInfo  - Struct containing:
%       .originalModel    - Original model reference
%       .stationSubset    - Original subset stations
%       .subsetIndices    - Original station indices
%       .throughputTable  - Computed throughputs for all states
%       .cutoffs          - Per-class cutoffs used
%       .stochCompSubset  - Stochastic complement for subset
%       .stochCompComplement - Stochastic complement for complement
%       .isolatedModel    - Isolated subnetwork model
%
% Example:
%   model = Network('Example');
%   % ... create stations and classes ...
%   [fesModel, fesStation, info] = ModelAdapter.aggregateFES(model, {queue2, queue3});
%   solver = SolverMVA(fesModel);
%   avgTable = solver.getAvgTable();
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% Input validation and defaults
if nargin < 3
    options = struct();
end
if ~isfield(options, 'solver')
    options.solver = 'mva';
end
if ~isfield(options, 'cutoffs')
    options.cutoffs = [];
end
if ~isfield(options, 'verbose')
    options.verbose = false;
end

%% Get network structure
sn = model.getStruct();
M = sn.nstations;
K = sn.nclasses;

% Get station indices for subset
modelNodes = model.getNodes();
origStationIdxs = model.getStationIndexes();
modelStations = modelNodes(origStationIdxs);
subsetIndices = zeros(1, length(stationSubset));
for i = 1:length(stationSubset)
    for j = 1:M
        if stationSubset{i} == modelStations{j}
            subsetIndices(i) = j;
            break;
        end
    end
end

% Validate inputs using sn struct and indices
[isValid, errorMsg] = fes_validate(sn, subsetIndices);
if ~isValid
    line_error(mfilename, errorMsg);
end

% Get complement indices (stations not in subset)
allIndices = 1:M;
complementIndices = setdiff(allIndices, subsetIndices);

% Set cutoffs if not provided
if isempty(options.cutoffs)
    % Default: use total jobs per class
    options.cutoffs = sn.njobs';
end
cutoffs = options.cutoffs;

if options.verbose
    fprintf('FES Aggregation: %d subset stations, %d complement stations, %d classes\n', ...
        length(subsetIndices), length(complementIndices), K);
end

%% Compute stochastic complement routing matrices
% The routing matrix sn.rt is indexed by (stateful_node-1)*K + class
% We need to partition it by stations

% Build index sets for rt matrix
subsetRtIndices = [];
for i = subsetIndices
    isf = sn.stationToStateful(i);
    subsetRtIndices = [subsetRtIndices, ((isf-1)*K + 1):(isf*K)];
end

complementRtIndices = [];
for i = complementIndices
    isf = sn.stationToStateful(i);
    complementRtIndices = [complementRtIndices, ((isf-1)*K + 1):(isf*K)];
end

% Compute stochastic complement for subset (routing within subset only)
% S = P11 + P12 * inv(I - P22) * P21
rt = sn.rt;
stochCompSubset = dtmc_stochcomp(rt, subsetRtIndices);

% Compute stochastic complement for complement
stochCompComplement = dtmc_stochcomp(rt, complementRtIndices);

if options.verbose
    fprintf('Stochastic complements computed: subset (%dx%d), complement (%dx%d)\n', ...
        size(stochCompSubset, 1), size(stochCompSubset, 2), ...
        size(stochCompComplement, 1), size(stochCompComplement, 2));
end

%% Build isolated subnetwork data and compute throughputs
[L_iso, mi_iso, visits_iso, isDelay_iso] = fes_build_isolated(sn, subsetIndices, stochCompSubset);

if options.verbose
    fprintf('Isolated subnetwork data extracted for %d stations.\n', length(subsetIndices));
end

% Compute throughputs for all population states
scalingTable = fes_compute_throughputs(L_iso, mi_iso, isDelay_iso, cutoffs, options);

if options.verbose
    fprintf('Throughput table computed for %d states.\n', prod(cutoffs + 1));
end

%% Create the FES model
fesModel = Network(sprintf('%s_FES', model.getName()));

% Copy complement stations to new model
nodeMap = cell(sn.nnodes, 1);
stationMap = cell(M, 1);

for i = complementIndices
    origStation = modelStations{i};
    nodeIdx = sn.stationToNode(i);

    switch class(origStation)
        case 'Queue'
            newStation = Queue(fesModel, origStation.name, origStation.schedStrategy);
            if ~isinf(origStation.numberOfServers)
                newStation.setNumberOfServers(origStation.numberOfServers);
            end
            if ~isempty(origStation.cap) && isfinite(origStation.cap)
                newStation.setCapacity(origStation.cap);
            end
        case 'Delay'
            newStation = Delay(fesModel, origStation.name);
        otherwise
            line_error(mfilename, sprintf('Unsupported station type %s.', class(origStation)));
    end

    nodeMap{nodeIdx} = newStation;
    stationMap{i} = newStation;
end

% Create the FES station (single Queue with PS scheduling)
fesStation = Queue(fesModel, 'FES', SchedStrategy.PS);
fesStation.setNumberOfServers(1);

% Create job classes
newClasses = cell(1, K);
modelClasses = model.classes;

% Choose reference station (FES or first complement station)
if ~isempty(complementIndices)
    refStation = stationMap{complementIndices(1)};
else
    refStation = fesStation;
end

for k = 1:K
    origClass = modelClasses{k};
    population = sn.njobs(k);
    newClass = ClosedClass(fesModel, origClass.name, population, refStation);
    newClasses{k} = newClass;
end

% Set service distributions for complement stations
for i = complementIndices
    newStation = stationMap{i};

    for k = 1:K
        origPH = sn.proc{i}{k};

        if isempty(origPH) || (iscell(origPH) && isempty(origPH{1})) || ...
           (iscell(origPH) && all(isnan(origPH{1}(:))))
            newStation.setService(newClasses{k}, Disabled.getInstance());
        else
            if iscell(origPH)
                T_matrix = origPH{1};  % Sub-generator (diagonal is -rate)
                t0_vector = origPH{2}; % Exit rate vector
                nPhases = size(T_matrix, 1);

                if nPhases == 1
                    rate = -T_matrix(1,1);
                    newStation.setService(newClasses{k}, Exp(rate));
                else
                    alpha = ones(1, nPhases) / nPhases;
                    dist = APH(alpha, T_matrix);
                    newStation.setService(newClasses{k}, dist);
                end
            else
                newStation.setService(newClasses{k}, Exp(sn.rates(i, k)));
            end
        end
    end
end

% Set LJCD on FES with per-class throughput scaling tables
% The scalingTable{k} contains throughputs for class k at each population state
%
% For multi-class FES, we use Limited Joint Class Dependence (LJCD) which allows
% each class to have its own state-dependent scaling factor.
% The FES service rate for class c in state (n1,...,nK) equals throughput_c(n).

% Set base service rate (will be scaled by LJCD)
% Use Exp(1.0) as base; the LJCD scaling provides the actual throughput rate
for k = 1:K
    fesStation.setService(newClasses{k}, Exp(1.0));
end

% Handle zeros in scaling tables (replace with small positive value)
for k = 1:K
    scalingTable{k}(scalingTable{k} < GlobalConstants.FineTol) = GlobalConstants.FineTol;
end

if options.verbose
    fprintf('Per-class scaling tables:\n');
    for k = 1:K
        fprintf('  Class %d: Max = %.6f, Min = %.6f\n', k, ...
            max(scalingTable{k}), min(scalingTable{k}(scalingTable{k} > GlobalConstants.FineTol)));
    end
end

% Set per-class joint dependence (LJCD)
% scalingTable is already a cell array {1 x K} where each cell is a linearized vector
fesStation.setJointClassDependence(scalingTable, cutoffs);

%% Build routing matrix for FES model
I_fes = length(fesModel.nodes);
P = fesModel.initRoutingMatrix();

% Get FES node index
fesNodeIdx = 0;
for n = 1:I_fes
    if fesModel.nodes{n} == fesStation
        fesNodeIdx = n;
        break;
    end
end

% Build node index mapping for complement stations
complementNodeMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
for i = complementIndices
    nodeIdx = sn.stationToNode(i);
    if ~isempty(nodeMap{nodeIdx})
        for n = 1:I_fes
            if fesModel.nodes{n} == nodeMap{nodeIdx}
                complementNodeMap(i) = n;
                break;
            end
        end
    end
end

% Set routing probabilities
% Routes within complement use stochastic complement
% Routes to/from subset go through FES
if options.verbose
    fprintf('Setting up routing for FES model:\n');
    fprintf('  FES node index: %d\n', fesNodeIdx);
    fprintf('  Complement station indices: %s\n', mat2str(complementIndices));
    fprintf('  Subset station indices: %s\n', mat2str(subsetIndices));
end

for k = 1:K
    P_k = zeros(I_fes, I_fes);

    % Compute visit ratios within subset for exit distribution
    % For closed network, solve pi * S = pi where S is stochastic complement for subset
    % For simplicity, assume uniform if can't compute (valid for symmetric cases)
    nSub = length(subsetIndices);
    visitRatios = ones(1, nSub) / nSub;  % Default: uniform

    % Routes within complement (DIRECT paths only, not through subset)
    for i = complementIndices
        if ~isKey(complementNodeMap, i)
            continue;
        end
        iNode = complementNodeMap(i);
        isf_i = sn.stationToStateful(i);

        for j = complementIndices
            if ~isKey(complementNodeMap, j)
                continue;
            end
            jNode = complementNodeMap(j);
            isf_j = sn.stationToStateful(j);

            % Use ORIGINAL routing for direct complement-to-complement paths
            rtIdx_i = (isf_i - 1) * K + k;
            rtIdx_j = (isf_j - 1) * K + k;

            if rtIdx_i <= size(rt, 1) && rtIdx_j <= size(rt, 2)
                prob = rt(rtIdx_i, rtIdx_j);
                if prob > GlobalConstants.FineTol
                    P_k(iNode, jNode) = prob;
                end
            end
        end

        % Routes from complement to subset -> FES
        for j = subsetIndices
            isf_j = sn.stationToStateful(j);
            rtIdx_i = (isf_i - 1) * K + k;
            rtIdx_j = (isf_j - 1) * K + k;

            if rtIdx_i <= size(rt, 1) && rtIdx_j <= size(rt, 2)
                prob = rt(rtIdx_i, rtIdx_j);
                if prob > GlobalConstants.FineTol
                    P_k(iNode, fesNodeIdx) = P_k(iNode, fesNodeIdx) + prob;
                end
            end
        end
    end

    % Routes from FES to complement
    % Weight by visit ratios within subset
    for j = complementIndices
        if ~isKey(complementNodeMap, j)
            continue;
        end
        jNode = complementNodeMap(j);
        isf_j = sn.stationToStateful(j);

        probSum = 0;
        for idx = 1:nSub
            i = subsetIndices(idx);
            isf_i = sn.stationToStateful(i);
            rtIdx_i = (isf_i - 1) * K + k;
            rtIdx_j = (isf_j - 1) * K + k;

            if rtIdx_i <= size(rt, 1) && rtIdx_j <= size(rt, 2)
                prob = rt(rtIdx_i, rtIdx_j);
                % Weight by visit ratio
                probSum = probSum + visitRatios(idx) * prob;
            end
        end

        if probSum > GlobalConstants.FineTol
            P_k(fesNodeIdx, jNode) = probSum;
        end
    end

    % No explicit self-loop on FES - internal routing is captured by LJD rates
    % (the state-dependent throughput already accounts for internal circulation)

    % Normalize rows
    for n = 1:I_fes
        rowSum = sum(P_k(n, :));
        if rowSum > GlobalConstants.FineTol
            P_k(n, :) = P_k(n, :) / rowSum;
        end
    end

    P{k, k} = P_k;

    if options.verbose
        fprintf('Routing matrix for class %d:\n', k);
        nodeNames = cell(1, I_fes);
        for n = 1:I_fes
            nodeNames{n} = fesModel.nodes{n}.name;
        end
        fprintf('  Nodes: %s\n', strjoin(nodeNames, ', '));
        for i = 1:I_fes
            fprintf('  %s: %s\n', nodeNames{i}, mat2str(P_k(i,:), 4));
        end
    end
end

% Link the model
fesModel.link(P);

%% Build deaggregation info
deaggInfo = struct();
deaggInfo.originalModel = model;
deaggInfo.stationSubset = stationSubset;
deaggInfo.subsetIndices = subsetIndices;
deaggInfo.complementIndices = complementIndices;
deaggInfo.throughputTable = scalingTable;
deaggInfo.cutoffs = cutoffs;
deaggInfo.stochCompSubset = stochCompSubset;
deaggInfo.stochCompComplement = stochCompComplement;
deaggInfo.isolatedDemands = L_iso;
deaggInfo.isolatedServers = mi_iso;
deaggInfo.isolatedVisits = visits_iso;
deaggInfo.isolatedIsDelay = isDelay_iso;
deaggInfo.fesNodeIdx = fesNodeIdx;

if options.verbose
    fprintf('FES model created with %d stations (1 FES + %d complement).\n', ...
        I_fes, length(complementIndices));
end

end
