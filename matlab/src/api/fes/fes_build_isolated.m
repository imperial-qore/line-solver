%{ @file fes_build_isolated.m
 %  @brief Builds isolated subnetwork data from a station subset for FES computation
 %
 %  @author LINE Development Team
%}

%{
 % @brief Builds isolated subnetwork data from a station subset
 %
 % @details
 % Extracts the data needed to analyze an isolated subnetwork containing
 % only the specified stations. This data is used with pfqn_mva to compute
 % throughputs for the Flow-Equivalent Server (FES) service rates.
 %
 % @par Syntax:
 % @code
 % [L, mi, visits, isDelay] = fes_build_isolated(sn, subsetIndices, stochCompS)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure (from model.getStruct())
 % <tr><td>subsetIndices<td>Array of station indices to include (1-based)
 % <tr><td>stochCompS<td>Stochastic complement routing matrix for subset (K*M_sub x K*M_sub)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>L<td>Service demands matrix (M_sub x K), L(i,k) = visits(i,k) / rate(i,k)
 % <tr><td>mi<td>Number of servers per station (1 x M_sub), Inf for Delay nodes
 % <tr><td>visits<td>Visit ratios matrix (M_sub x K) from stochastic complement
 % <tr><td>isDelay<td>Boolean array (1 x M_sub), true if station is a Delay node
 % </table>
%}
function [L, mi, visits, isDelay] = fes_build_isolated(sn, subsetIndices, stochCompS)

K = sn.nclasses;
M_sub = length(subsetIndices);

% Initialize outputs
L = zeros(M_sub, K);       % Service demands
mi = zeros(1, M_sub);      % Number of servers
isDelay = false(1, M_sub); % Delay node flags

% Extract station properties from sn
for i = 1:M_sub
    origIdx = subsetIndices(i);
    nodeIdx = sn.stationToNode(origIdx);
    nodeType = sn.nodetype(nodeIdx);

    % Check if this is a Delay node
    if nodeType == NodeType.Delay
        isDelay(i) = true;
        mi(i) = Inf;  % Infinite servers for Delay
    else
        isDelay(i) = false;
        % Get number of servers from nservers field
        mi(i) = sn.nservers(origIdx);
        if isnan(mi(i)) || mi(i) < 1
            mi(i) = 1;
        end
    end
end

% Compute visit ratios from stochastic complement routing matrix
% stochCompS is indexed by (station-1)*K + class
% Solve pi * S = pi for stationary distribution (by class)
visits = zeros(M_sub, K);

for k = 1:K
    % Extract routing sub-matrix for class k within subset
    % The stochCompS is organized as blocks: station i, class k -> station j, class k
    P_k = zeros(M_sub, M_sub);

    for i = 1:M_sub
        rowIdx = (i-1)*K + k;
        for j = 1:M_sub
            colIdx = (j-1)*K + k;
            if rowIdx <= size(stochCompS, 1) && colIdx <= size(stochCompS, 2)
                P_k(i, j) = stochCompS(rowIdx, colIdx);
            end
        end
    end

    % Normalize rows if needed
    for i = 1:M_sub
        rowSum = sum(P_k(i, :));
        if rowSum > GlobalConstants.FineTol && abs(rowSum - 1) > GlobalConstants.FineTol
            P_k(i, :) = P_k(i, :) / rowSum;
        elseif rowSum < GlobalConstants.FineTol
            % If no routing defined, route to self
            P_k(i, i) = 1.0;
        end
    end

    % Compute visit ratios by solving pi * P = pi
    % This is the stationary distribution of the embedded DTMC
    % Use (I - P' + e*e'/M)' * pi' = e'/M where e is ones vector
    I = eye(M_sub);
    e = ones(M_sub, 1);

    % Solve for stationary distribution
    A = I - P_k' + (e * e') / M_sub;

    if rank(A) == M_sub
        pi_k = (A' \ (e / M_sub))';
    else
        % If singular, use uniform distribution
        pi_k = ones(1, M_sub) / M_sub;
    end

    % Ensure non-negative and normalize
    pi_k = max(pi_k, 0);
    if sum(pi_k) > 0
        pi_k = pi_k / sum(pi_k);
    else
        pi_k = ones(1, M_sub) / M_sub;
    end

    visits(:, k) = pi_k(:);
end

% Compute service demands: L(i,k) = visits(i,k) * service_time(i,k)
for i = 1:M_sub
    origIdx = subsetIndices(i);

    for k = 1:K
        % Get service rate from sn
        rate = sn.rates(origIdx, k);

        if isnan(rate) || rate <= 0 || isinf(rate)
            % Disabled or invalid service
            L(i, k) = 0;
        else
            % Service demand = visits * service time
            % For MVA, we normalize visits relative to first station
            serviceTime = 1 / rate;
            L(i, k) = visits(i, k) * serviceTime;
        end
    end
end

% Normalize demands so that total visit to first station is 1
% This is the standard MVA convention
for k = 1:K
    if visits(1, k) > 0
        L(:, k) = L(:, k) / visits(1, k);
    end
end

end
