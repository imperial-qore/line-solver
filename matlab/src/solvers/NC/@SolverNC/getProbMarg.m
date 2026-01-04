function [Pmarg, logPmarg] = getProbMarg(self, node)
% [PMARG, LOGPMARG] = GETPROBMARG(NODE)
%
% Probability distribution for TOTAL queue length at a station.
% Returns P(n total jobs) for n=0,1,...,N, summing over all class combinations.
%
% Compare with getProbAggr: returns probability of a specific per-class
% distribution, e.g., P(2 class-1, 1 class-2).
%
% Input:
%   node - Queue or node object
%
% Output:
%   Pmarg    - Vector where Pmarg(n+1) = P(n total jobs at node)
%   logPmarg - Log probabilities for numerical stability

if GlobalConstants.DummyMode
    Pmarg = NaN;
    logPmarg = NaN;
    return
end

% Get network structure
sn = self.getStruct;

% Store original node parameter and convert to station index
if isa(node, 'Node')
    node_obj = node;
    ist = sn.nodeToStation(node.index);
else
    % If given station index, need to find corresponding node
    ist = node;
    % Find node with this station index
    for i = 1:length(sn.nodetype)
        if sn.nodeToStation(i) == ist
            node_obj = sn.nodes{i};
            break;
        end
    end
end

% Get population and class information
N = sn.njobs;
if ~all(isfinite(N))
    line_error(mfilename, 'getProbMarg not yet implemented for models with open classes.');
end

R = sn.nclasses;  % Number of classes
Ntotal = sum(N);  % Total population

% Initialize output vectors
Pmarg = zeros(1, Ntotal + 1);
logPmarg = -Inf * ones(1, Ntotal + 1);

% For each possible total queue length
for n = 0:Ntotal
    % Enumerate all ways to partition n jobs across R classes
    % subject to constraint: n_r <= N(r) for each class
    partitions = generate_partitions(n, R, N);

    % Sum probabilities over all partitions
    log_prob_n = -Inf;
    for p = 1:size(partitions, 1)
        state_vec = partitions(p, :);

        % Get probability for this specific class distribution
        % Note: getProbAggr already supports LD models via pfqn_ncld
        % Gets unnormalized values (we normalize at the end)
        prob = self.getProbAggr(node_obj, state_vec);

        if prob > 0
            log_p = log(prob);
            % Log-sum-exp trick for numerical stability
            if isinf(log_prob_n)
                log_prob_n = log_p;
            else
                log_prob_n = log_prob_n + log(1 + exp(log_p - log_prob_n));
            end
        end
    end

    % Store results
    if isfinite(log_prob_n)
        Pmarg(n + 1) = exp(log_prob_n);
        logPmarg(n + 1) = log_prob_n;
    else
        Pmarg(n + 1) = 0;
        logPmarg(n + 1) = -Inf;
    end
end

% Normalize to ensure probabilities sum to 1
total_prob = sum(Pmarg);
if total_prob > 0 && abs(total_prob - 1.0) > 1e-10
    Pmarg = Pmarg / total_prob;
    logPmarg = logPmarg - log(total_prob);
end

end

function partitions = generate_partitions(n, R, N_max)
% GENERATE_PARTITIONS Enumerate all ways to partition n into R non-negative integers
%
% @param n Total to partition
% @param R Number of parts
% @param N_max Maximum value for each part (vector of length R)
% @return partitions Matrix where each row is a valid partition

if R == 1
    % Base case: single class
    if n <= N_max(1)
        partitions = n;
    else
        partitions = zeros(0, 1);
    end
    return
end

% Recursive case: try all possible values for first class
partitions = [];
for n1 = 0:min(n, N_max(1))
    % Recursively partition remaining jobs across remaining classes
    sub_partitions = generate_partitions(n - n1, R - 1, N_max(2:end));

    if ~isempty(sub_partitions)
        % Prepend n1 to each sub-partition
        partitions = [partitions; n1 * ones(size(sub_partitions, 1), 1), sub_partitions];
    end
end

end
