%{ @file fes_compute_throughputs.m
 %  @brief Computes throughputs for all population states for FES scaling table
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes throughputs for all population states
 %
 % @details
 % Enumerates all population states (n1, n2, ..., nK) up to the specified
 % cutoffs and solves the isolated subnetwork with MVA to get per-class throughputs.
 % These throughputs are stored in a linearized scaling table for use as
 % Limited Joint Dependence (LJD) service rates in the Flow-Equivalent Server.
 %
 % @par Syntax:
 % @code
 % scalingTable = fes_compute_throughputs(L, mi, isDelay, cutoffs, options)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>L<td>Service demands matrix (M_sub x K) from fes_build_isolated
 % <tr><td>mi<td>Number of servers per station (1 x M_sub), Inf for Delay
 % <tr><td>isDelay<td>Boolean array (1 x M_sub), true if station is a Delay
 % <tr><td>cutoffs<td>Per-class population cutoffs [N1, N2, ..., NK]
 % <tr><td>options<td>Options struct with field .verbose (default: false)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>scalingTable<td>Cell array {class} of linearized throughput vectors
 % </table>
%}
function scalingTable = fes_compute_throughputs(L, mi, isDelay, cutoffs, options)

if nargin < 5
    options = struct();
end
if ~isfield(options, 'verbose')
    options.verbose = false;
end

% Get dimensions
[M_sub, K] = size(L);

% Separate Queue and Delay stations for MVA
% Queue stations go into L_queue, Delay stations contribute to Z (think time)
queueIdx = find(~isDelay);
delayIdx = find(isDelay);

M_queue = length(queueIdx);

% Extract Queue station demands and servers
if M_queue > 0
    L_queue = L(queueIdx, :);
    mi_queue = mi(queueIdx);
else
    % No Queue stations - all Delay
    L_queue = zeros(0, K);
    mi_queue = [];
end

% Compute think time Z from Delay stations
% Z(k) = sum of service demands at Delay stations for class k
if ~isempty(delayIdx)
    Z = sum(L(delayIdx, :), 1);
else
    Z = zeros(1, K);
end

% Compute table size
tableSize = prod(cutoffs + 1);

% Initialize scaling tables (one per class)
scalingTable = cell(1, K);
for k = 1:K
    scalingTable{k} = zeros(1, tableSize);
end

if options.verbose
    fprintf('Computing FES throughputs for %d states (%d Queue, %d Delay stations)...\n', ...
        tableSize, M_queue, length(delayIdx));
end

% Enumerate all states using BFS
stateQueue = {zeros(1, K)}; % Start with zero population
% Each population state maps to a unique linearized index in 1..tableSize,
% so a logical array suffices to track which states have been visited.
visited = false(1, tableSize);
nStates = 0;

while ~isempty(stateQueue)
    nvec = stateQueue{1};
    stateQueue(1) = [];

    % Linearized index uniquely identifies this state
    idx = ljd_linearize(nvec, cutoffs);
    if visited(idx)
        continue;
    end
    visited(idx) = true;
    nStates = nStates + 1;

    % Compute throughput for this state
    totalPop = sum(nvec);

    if totalPop == 0
        % Zero population: throughput is 0 for all classes
        for k = 1:K
            scalingTable{k}(idx) = 0;
        end
    else
        % PFQN_MVA's MI is NOT a server count -- it enters only as the additive
        % term of C(i,s)=L(i,s)*(mi(i)+Qarv), so passing the real multiplicity
        % INFLATES the residence time instead of adding servers. Route through
        % PFQN_MVAMS, which forwards to PFQN_MVA when every station is a single
        % server and to the load-dependent recursion with mu(i,n)=min(n,S(i))
        % when one is not. See _kb/03-api-layer.md (pfqn_mva: mi is not S).
        if M_queue > 0
            [XN, ~, ~, ~] = pfqn_mvams(zeros(1,K), L_queue, nvec, Z, ...
                ones(M_queue,1), mi_queue(:));
        else
            % Only Delay stations - throughput is N(k) / Z(k) for each class
            XN = zeros(1, K);
            for k = 1:K
                if nvec(k) > 0 && Z(k) > 0
                    XN(k) = nvec(k) / Z(k);
                end
            end
        end

        % Store throughputs. A failure here is NOT caught: an FES table silently
        % filled with zeros is a wrong aggregate, not a degraded one, and every
        % downstream beta_r(n) reads it as "the subnetwork serves nothing".
        for k = 1:K
            if nvec(k) > 0
                scalingTable{k}(idx) = XN(k);
            else
                scalingTable{k}(idx) = 0;
            end
        end
    end

    % Add neighboring states to queue
    for k = 1:K
        if nvec(k) < cutoffs(k)
            newState = nvec;
            newState(k) = newState(k) + 1;
            newIdx = ljd_linearize(newState, cutoffs);
            if ~visited(newIdx)
                stateQueue{end+1} = newState;
            end
        end
    end
end

if options.verbose
    fprintf('Computed throughputs for %d states.\n', nStates);
end

end
