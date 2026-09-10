function [pi,pis,pi0,scc,isrec] = ctmc_solve_reducible_blkdecomp(Q,pin,options)
% [PI,PIS,PI0,SCC,ISREC] = CTMC_SOLVE_REDUCIBLE_BLKDECOMP(Q, PIN, OPTIONS)
%
% Compute limiting distribution for a CTMC with reducible generator Q
% using direct block decomposition on the generator matrix.
%
% Algorithm:
%   1. Decompose states into transient and recurrent classes via SCC
%   2. For transient states: solve sojourn * Q_tt = -p0_t for expected sojourn
%   3. Compute hitting probabilities: hit = sojourn * Q_ta + p0_r
%   4. For each recurrent class: solve pi_c * Q_cc = 0, scale by hitting prob
%
% Input:
% Q:       infinitesimal generator matrix
% pin:     initial distribution vector, set to [] if not available
% options: struct where options.tol sets the tolerance
%
% Output:
% pi:    limiting distribution (1 x N)
%        - For an ergodic CTMC, this is the unique limiting distribution.
%        - For a reducible CTMC:
%            - if PIN is supplied, this is the exact limiting distribution
%              for that initial vector: each recurrent class (BSCC) carries
%              its absorption probability from PIN and every transient state
%              carries zero.
%            - otherwise, if there is a single transient SCC then this is the
%              limiting distribution when starting uniformly within it, and
%              failing that pi is the weighted average of pis rows.
% pis:   limiting distribution given initialization in a single SCC (numSCC x N)
% pi0:   starting distribution for each row of pis (numSCC x N), so PIN is the
%        vector the caller supplies and PI0 the ones this routine builds
% scc:   mapping of each state of Q to its SCC index
% isrec: element i is true if SCC i is recurrent

N = size(Q, 1);

if nargin < 3
    options = struct('tol', 1e-12);
end

if nargin < 2
    pin = [];
end

% Ensure valid generator
Q = ctmc_makeinfgen(Q);

% Build adjacency from off-diagonal positive entries
Adj = Q;
Adj(1:N+1:end) = 0;
% arc by MAGNITUDE, never by sign: an ME generator embeds genuinely negative off-diagonals -- see _kb/11-conventions-and-gotchas.md
[scc, isrec] = stronglyconncomp(abs(Adj) > GlobalConstants.ArcTol);
numSCC = max(scc);

% Irreducible case: use standard solver
if numSCC == 1
    pi = ctmc_solve(Q, i_solveopts(options));
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
Q_tt = [];
Q_ta = [];
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
    pis(s, :) = absorb_limiting(p0, Q, N, nt, nr, trans_states, rec_states, ...
        scc_idx, rec_scc_ids, Q_tt, Q_ta, options);
end

if isempty(pin)
    % No initial vector: mix the uniform-start rows, weighting the SCCs equally.
    % An SCC holding a state whose column of Q is entirely zero is EXCLUDED: such
    % a state has no transition in and none out, so it is isolated and carries no
    % dynamics to start from. (This is not the same as absorbing, which has
    % incoming transitions and a zero off-diagonal ROW.)
    pinl = ones(1, numSCC);
    col_sums = sum(abs(Q), 1);
    for j = find(col_sums < 1e-12)
        pinl(scc(j)) = 0;
    end
    if sum(pinl) > 0
        pinl = pinl / sum(pinl);
    else
        pinl = ones(1, numSCC) / numSCC;
    end

    pi = zeros(1, N);
    for i = 1:numSCC
        if pinl(i) > 0
            pi = pi + pis(i, :) * pinl(i);
        end
    end

    % Special case: single transient SCC without explicit initial distribution
    if isscalar(trans_scc_ids)
        pi = pis(trans_scc_ids, :);
    end
else
    % An initial vector is available, so the exact absorption probabilities can
    % be computed from it directly. Lumping pin onto its SCCs and mixing the
    % pis rows would instead assume a uniform start within each SCC, which
    % differs from the truth whenever pin puts mass on a transient SCC holding
    % more than one state (states of the same transient SCC reach the recurrent
    % classes with different probabilities).
    p0 = pin(:)';
    total0 = sum(p0);
    if total0 > 0
        p0 = p0 / total0;
    end
    pi = absorb_limiting(p0, Q, N, nt, nr, trans_states, rec_states, ...
        scc_idx, rec_scc_ids, Q_tt, Q_ta, options);
end

% Normalize
total = sum(pi);
if total > 0
    pi = pi / total;
end

end

function piv = absorb_limiting(p0, Q, N, nt, nr, trans_states, rec_states, ...
    scc_idx, rec_scc_ids, Q_tt, Q_ta, options)
% OPTIONS IS A PARAMETER, not a variable this subfunction can see from its
% caller: the inner ctmc_solve below asks for i_solveopts(options), and without
% it in scope MATLAB resolved the name on the PATH instead. The arm is reached
% only by a recurrent class of more than one state carrying non-negligible
% hitting probability, so the defect sat latent until a model with routing
% state (RROBIN) produced one.
% Limiting distribution reached from the initial vector P0: the mass absorbed
% in each recurrent class (BSCC) redistributed over that class according to its
% own stationary vector. Transient states receive zero.
piv = zeros(1, N);

% Absorption probabilities into the recurrent states
hit = zeros(1, nr);
if nt > 0
    p0_t = p0(trans_states);
    if any(abs(p0_t) > 0)
        % Solve sojourn * Q_tt = -p0_t for expected sojourn in transient states.
        % Q_tt is non-singular (Hurwitz) for transient states. Above the
        % dispatch threshold the transient block is what the direct
        % factorization cannot hold; it remains the fallback.
        sojourn = [];
        if nt > 6000
            [xg,gflag] = ctmc_gmres(Q_tt', (-p0_t(:)));
            if gflag == 0
                sojourn = xg';
            else
                % Short-recurrence retry before the cubic factorization, as in ctmc_solve.
                [xb,bflag] = ctmc_bicgstab(Q_tt', (-p0_t(:)));
                if bflag == 0
                    sojourn = xb';
                end
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
        piv(idx_c) = reachprob;
    else
        % Solve pi_c * Q_cc = 0 within this recurrent class
        pi_c = ctmc_solve(Q(idx_c, idx_c), i_solveopts(options));
        piv(idx_c) = pi_c * reachprob;
    end
end

end

function sopts = i_solveopts(options)
% Options for the inner CTMC_SOLVE: FORCE, plus the caller's backend choice.
%
% The struct used to be struct('force',true) alone, which DROPPED
% options.config on the way down -- so a caller asking for a Krylov or GPU
% linear solve (config.linsolver) got the direct factorization and no word of
% it. Only config and verbose are forwarded: the method name is deliberately
% not, since it names the state-space construction of the calling solver and
% means nothing to the linear solve.
sopts = struct('force', true);
if isstruct(options)
    if isfield(options,'config') && isstruct(options.config)
        sopts.config = options.config;
    end
    if isfield(options,'verbose') && ~isempty(options.verbose)
        sopts.verbose = options.verbose;
    end
    if isfield(options,'iter_max') && ~isempty(options.iter_max)
        sopts.iter_max = options.iter_max;
    end
end
end
