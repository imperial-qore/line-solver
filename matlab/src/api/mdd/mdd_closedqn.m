function out = mdd_closedqn(mu, P, servers, N, options)
% OUT = MDD_CLOSEDQN(MU, P, SERVERS, N)
% OUT = MDD_CLOSEDQN(MU, P, SERVERS, N, OPTIONS)
% Solve a single-class closed exponential queueing network whose CTMC state
% space (reachable occupancy vectors) is stored in a Multi-valued Decision
% Diagram (MDD) instead of an explicit state list.
%
% The reachable set is generated with MDD_REACHSET and the generator matrix is
% assembled using the MDD's O(K) state-indexing (MDD/index), so no explicit
% (|S| x width) state matrix is ever materialised -- the diagram is the store.
%
% -- Input
% MU      : 1 x M vector of per-station exponential service rates
% P       : M x M Markovian routing matrix (row-stochastic, irreducible)
% SERVERS : 1 x M number of servers per station (Inf for a delay/IS station)
% N       : closed population (number of jobs)
% OPTIONS : (optional) struct; OPTIONS.verbose prints a storage summary,
%           OPTIONS.mdd reuses an already-built reachable set (as returned in
%           OUT.mdd) instead of regenerating it
% -- Output
% OUT     : struct with fields
%             mdd    - the MDD holding the reachable occupancy set
%             Q      - sparse CTMC generator, rows aligned to MDD/index order
%             pi     - stationary distribution (1 x |S|)
%             states - |S| x M occupancy matrix in MDD/index order
%             QLen   - 1 x M mean number of jobs per station
%             U      - 1 x M utilisation (busy servers / servers; mean busy for IS)
%             X      - 1 x M per-station throughput
%             stats  - MDD storage statistics (see MDD/stats)
%             times  - phase timings in seconds: reach (reachable-set build, 0
%                      when OPTIONS.mdd is supplied), gen (generator assembly),
%                      solve (ctmc_solve), metrics (performance measures)
%
% -- Remarks
% For single-class exponential stations the aggregated (occupancy) chain is
% exact: the rate from n to n-e_i+e_j is mu(i)*min(n_i,servers(i))*P(i,j) for
% n_i>0. This matches LINE's SolverCTMC on the same model.
%
% See also: MDD, mdd_reachset, solver_ctmc, ctmc_solve.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5 || isempty(options)
    options = struct('verbose', false);
end
if ~isfield(options, 'verbose'), options.verbose = false; end

M = numel(mu);
mu = mu(:)';
servers = servers(:)';
domain = (N + 1) * ones(1, M);

% events: a job completes at station i (rate mu(i)*min(n_i,servers(i))) and
% routes to station j with probability P(i,j); self-routing i==j leaves the
% occupancy vector unchanged and is skipped.
[ii, jj, pij] = find(P);
keep = ii ~= jj;
ii = ii(keep); jj = jj(keep); pij = pij(keep);

% all jobs start at station 1; an irreducible routing chain makes every
% composition of N over M stations reachable.
init = zeros(1, M); init(1) = N;

nextfun = @(s) i_next(s, ii, jj);
t0 = tic;
if isfield(options, 'mdd') && ~isempty(options.mdd)
    mdd = options.mdd;             % caller already built the reachable set
    times.reach = 0;
else
    mdd = mdd_reachset(domain, init, nextfun);
    times.reach = toc(t0);
end

t0 = tic;
n = mdd.cardinality();
S = mdd.enumerate();               % n x M occupancy states, in MDD/index order

% assemble the generator directly from the MDD state indexing
E = numel(ii);
ri = zeros(n * E, 1); ci = zeros(n * E, 1); vv = zeros(n * E, 1);
c = 0;
for s = 1:n
    st = S(s, :);
    row = mdd.index(st);           % 0-based
    for a = 1:E
        i = ii(a);
        if st(i) > 0
            rate = mu(i) * min(st(i), servers(i)) * pij(a);
            t = st; t(i) = t(i) - 1; t(jj(a)) = t(jj(a)) + 1;
            c = c + 1;
            ri(c) = row + 1;
            ci(c) = mdd.index(t) + 1;
            vv(c) = rate;
        end
    end
end
Q = sparse(ri(1:c), ci(1:c), vv(1:c), n, n);
Q = ctmc_makeinfgen(Q);
times.gen = toc(t0);

t0 = tic;
if isfield(options, 'ctmcmethod') && ~isempty(options.ctmcmethod)
    p = ctmc_solve(Q, struct('method', options.ctmcmethod, 'verbose', 0));
else
    p = ctmc_solve(Q);
end
p = p(:)';
times.solve = toc(t0);

% performance metrics
t0 = tic;
QLen = p * S;                                  % mean occupancy per station
X = zeros(1, M); U = zeros(1, M);
for i = 1:M
    busy = min(S(:, i), servers(i));           % busy servers in each state
    X(i) = mu(i) * (p * busy);                 % throughput = mean completion rate
    if isinf(servers(i))
        U(i) = p * S(:, i);                    % mean number busy (IS station)
    else
        U(i) = (p * busy) / servers(i);        % server utilisation
    end
end

times.metrics = toc(t0);

out.mdd = mdd;
out.times = times;
out.Q = Q;
out.pi = p;
out.states = S;
out.QLen = QLen;
out.U = U;
out.X = X;
out.stats = mdd.stats();

if options.verbose
    s = out.stats;
    line_printf('\nMDD-stored closed QN: %d stations, N=%d\n', M, N);
    line_printf('  reachable states |S| = %d\n', s.numStates);
    line_printf('  MDD nodes            = %d  (%s per level)\n', ...
        s.numNodes, strtrim(sprintf('%d ', s.nodesPerLevel)));
    line_printf('  storage              = %d ints vs %d explicit (%.1fx)\n', ...
        s.mddInts, s.explicitInts, s.compression);
end
end

function T = i_next(s, ii, jj)
% successor occupancy vectors of s under the routing events (ii -> jj)
E = numel(ii);
T = zeros(E, numel(s));
m = 0;
for a = 1:E
    i = ii(a);
    if s(i) > 0
        t = s; t(i) = t(i) - 1; t(jj(a)) = t(jj(a)) + 1;
        m = m + 1;
        T(m, :) = t;
    end
end
T = T(1:m, :);
end
