function [pi,pis,pi0,scc,isrec] = ctmc_solve_reducible_blkdecomp(Q,pi0,options)
% [PI,PIS,PI0,SCC,ISREC] = CTMC_SOLVE_REDUCIBLE_BLKDECOMP(Q, PI0, OPTIONS)
%
% Compute limiting distribution for a CTMC with reducible generator Q
% using direct block decomposition on the generator matrix.
%
% Algorithm:
%   1. Decompose states into transient and recurrent classes via SCC
%   2. For transient states: solve n * Q_tt = -p0_t for expected sojourn
%   3. Compute hitting probabilities: h = n * Q_ta + p0_r
%   4. For each recurrent class: solve pi_c * Q_cc = 0, scale by hitting prob
%
% Input:
% Q:       infinitesimal generator matrix
% pi0:     initial distribution vector, set to [] if not available
% options: struct where options.tol sets the tolerance
%
% Output:
% pi:    limiting distribution (1 x N)
%        - For an ergodic CTMC, this is the unique limiting distribution.
%        - For a reducible CTMC:
%            - if there is a single transient SCC then this is the limiting
%              distribution when starting uniformly within it.
%            - otherwise pi is the weighted average of pis rows.
% pis:   limiting distribution given initialization in a single SCC (numSCC x N)
% pi0:   starting distribution for each row of pis (numSCC x N)
% scc:   mapping of each state of Q to its SCC index
% isrec: element i is true if SCC i is recurrent

N = size(Q, 1);

if nargin < 3
    options = struct('tol', 1e-12);
end

if nargin < 2
    pin = [];
else
    pin = pi0;
end

% Ensure valid generator
Q = ctmc_makeinfgen(Q);

% Build adjacency from off-diagonal positive entries
Adj = Q;
Adj(1:N+1:end) = 0;
[scc, isrec] = stronglyconncomp(Adj > 0);
numSCC = max(scc);

% Irreducible case: use standard solver
if numSCC == 1
    pi = ctmc_solve(Q, struct('force',true));
    pis = pi;
    pi0 = [];
    return
end

% Build SCC index sets
scc_idx = cell(1, numSCC);
for i = 1:numSCC
    scc_idx{i} = find(scc == i);
end

% Classify SCCs
trans_scc_ids = find(~isrec);
rec_scc_ids = find(isrec);

% Gather ordered state indices
trans_states = [];
for i = trans_scc_ids
    trans_states = [trans_states, scc_idx{i}];
end
trans_states = sort(trans_states);

rec_states = [];
for i = rec_scc_ids
    rec_states = [rec_states, scc_idx{i}];
end
rec_states = sort(rec_states);

nt = length(trans_states);
nr = length(rec_states);

% Extract Q sub-blocks (only when transient states exist)
if nt > 0
    Q_tt = Q(trans_states, trans_states);
    Q_ta = Q(trans_states, rec_states);
end

% Compute per-SCC limiting distributions
pis = zeros(numSCC, N);
pi0 = zeros(numSCC, N);

for s = 1:numSCC
    % Starting distribution: uniform within SCC s
    p0 = zeros(1, N);
    p0(scc_idx{s}) = 1 / length(scc_idx{s});
    pi0(s, :) = p0;

    % Compute absorption probabilities into recurrent states
    hit = zeros(1, nr);
    if nt > 0
        p0_t = p0(trans_states);
        if any(abs(p0_t) > 0)
            % Solve n * Q_tt = -p0_t for expected sojourn in transient states
            % Q_tt is non-singular (Hurwitz) for transient states.
            % Above the dispatch threshold the transient block is what the direct
            % factorization cannot hold; it remains the fallback.
            sojourn = [];
            if nt > 6000
                [xg,gflag] = ctmc_gmres(Q_tt', (-p0_t(:)));
                if gflag == 0
                    sojourn = xg';
                end
            end
            if isempty(sojourn)
                sojourn = (-p0_t) / Q_tt;
            end
            hit = sojourn * Q_ta;
        end
    end

    % Add initial mass already in recurrent states
    hit = hit + p0(rec_states);

    % Solve steady-state per recurrent class, scaled by hitting probability
    for c = rec_scc_ids
        idx_c = scc_idx{c};
        [~, loc] = ismember(idx_c, rec_states);
        reachprob = sum(hit(loc));
        if reachprob < 1e-15
            continue
        end
        if length(idx_c) == 1
            % Absorbing state: hitting probability IS the final probability
            pis(s, idx_c) = reachprob;
        else
            % Solve pi_c * Q_cc = 0 within this recurrent class
            pi_c = ctmc_solve(Q(idx_c, idx_c), struct('force',true));
            pis(s, idx_c) = pi_c * reachprob;
        end
    end
end

% Compute initial SCC probabilities for weighted average
if isempty(pin)
    pinl = ones(1, numSCC);
    % Zero out SCCs containing states with zero column sums (no incoming)
    col_sums = sum(abs(Q), 1);
    for j = find(col_sums < 1e-12)
        pinl(scc(j)) = 0;
    end
    if sum(pinl) > 0
        pinl = pinl / sum(pinl);
    else
        pinl = ones(1, numSCC) / numSCC;
    end
else
    pinl = zeros(1, numSCC);
    for i = 1:numSCC
        pinl(i) = sum(pin(scc_idx{i}));
    end
end

% Weighted average over starting SCCs
pi = zeros(1, N);
for i = 1:numSCC
    if pinl(i) > 0
        pi = pi + pis(i, :) * pinl(i);
    end
end

% Special case: single transient SCC without explicit initial distribution
if isscalar(trans_scc_ids) && isempty(pin)
    pi = pis(trans_scc_ids, :);
end

% Normalize
total = sum(pi);
if total > 0
    pi = pi / total;
end

end
