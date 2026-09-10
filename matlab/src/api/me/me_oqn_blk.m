function [Q, W, T, U, Ca, Cd, PBa, lambda, iter] = me_oqn_blk(M, lambda0, Ca0, mu, Cs, P, c, N, blockrule, options)
%ME_OQN_BLK Maximum Entropy algorithm for single-class open queueing
%networks with finite buffers, loss and transfer blocking
%
% Extends ME_OQN to open networks in which a station has a finite buffer.
% Two per-station policies are supported:
%
%   loss (blockrule = 0)  - a job that finds the destination full is
%                           discarded. Each station is then a censored
%                           GE/GE/c/0;N queue, and the network is the ME
%                           decomposition of Kouvatsos (1994), Section 4.
%   transfer blocking     - a job that completes service at station i and
%   (blockrule = 1)         finds the destination j full is held in i's
%                           server, which cannot serve anyone else until
%                           j has room (blocking after service, BAS).
%
% Transfer blocking is not work conserving, so a product-form
% approximation cannot be applied to the network as it stands. Following
% Tahilramani, Manjunath and Bose (1999) the network is first made work
% conserving by inserting a GE/GE/inf HOLDING NODE h_ij on every routing
% pair with p_ij > 0 and a finite-buffer destination j. The holding node
% absorbs the blocked job, so station i's server is released; the delay it
% introduces is the residual life of the minimum of the c_j service times
% in progress at j, inflated geometrically because the released job may
% find j full again. Station i's own service time is inflated by the same
% blocking probability so that the jobs queued behind the blocked one
% still see the server as busy. The expanded network is work conserving
% and is solved node by node with the censored ME queue of ME_GEGECN,
% iterating over the blocking probabilities and the first two moments of
% the flows until they converge.
%
% INPUTS:
%   M         - Number of stations
%   lambda0   - External arrival rates [M x 1]
%   Ca0       - External interarrival scv [M x 1] (>= 1 where lambda0 > 0)
%   mu        - Service rates [M x 1]
%   Cs        - Service scv [M x 1] (>= 1 at every finite-buffer station)
%   P         - Routing probabilities [M x M], P(i,j) = p_ij. Row sums
%               below one send the residual flow out of the network
%   c         - Servers per station [M x 1]; Inf marks an infinite server
%   N         - Buffer capacity per station [M x 1] in jobs, in service
%               included; Inf marks an unbounded buffer
%   blockrule - Policy at each finite-buffer station [M x 1]: 0 = loss,
%               1 = transfer blocking (BAS)
%   options   - (optional) struct with fields .tol (default 1e-6),
%               .maxiter (default 1000), .verbose (default false),
%               .damping (default 0.5), applied to the blocking
%               probabilities as in the relaxation scheme of the source
%
% OUTPUTS:
%   Q      - Mean number of jobs at each station [M x 1], including the
%            jobs held blocked in that station's servers
%   W      - Mean response time [M x 1], Q ./ T
%   T      - Throughput (carried flow) [M x 1]
%   U      - Utilization [M x 1], mean fraction of busy servers, the jobs
%            held blocked counted as occupying their server
%   Ca     - Interarrival scv of the offered flow at each station [M x 1]
%   Cd     - Interdeparture scv at each station [M x 1]
%   PBa    - Probability that an arrival at each station finds it full,
%            averaged over the incoming streams [M x 1]
%   lambda - Offered arrival rate at each station [M x 1], attempts
%            included
%   iter   - Number of iterations performed
%
% References:
%   D.D. Kouvatsos, "Entropy Maximisation and Queueing Network Models",
%   Annals of Operations Research, 48:63-126, 1994, Section 4.
%   H. Tahilramani, D. Manjunath, S.K. Bose, "Approximate analysis of open
%   network of GE/GE/m/N queues with transfer blocking", MASCOTS 1999,
%   164-171.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 9 || isempty(blockrule)
    blockrule = zeros(M, 1);
end
if nargin < 10
    options = struct();
end
if ~isfield(options, 'tol'), options.tol = 1e-6; end
if ~isfield(options, 'maxiter'), options.maxiter = 1000; end
if ~isfield(options, 'verbose'), options.verbose = false; end
if ~isfield(options, 'damping'), options.damping = 0.5; end

lambda0 = lambda0(:);
Ca0 = Ca0(:);
mu = mu(:);
Cs = Cs(:);
c = c(:);
N = N(:);
blockrule = blockrule(:);

finiteBuf = isfinite(N) & isfinite(c);
bas = finiteBuf & (blockrule == 1);

for i = 1:M
    if finiteBuf(i) && Cs(i) < 1 - 1e-12
        line_error(mfilename, 'MEM with finite buffers requires a service scv of at least 1 at station %d: the GE distribution is not defined for scv < 1.', i);
    end
    if lambda0(i) > 0 && Ca0(i) < 1 - 1e-12
        line_error(mfilename, 'MEM with finite buffers requires an external interarrival scv of at least 1 at station %d.', i);
    end
end

% see _kb/03-api-layer.md (me_oqn_blk -- GE-type OQN with blocking) for rationale
Pf = P;
muf = mu;
Csf = Cs;
for i = 1:M
    pii = P(i, i);
    if pii > 0
        muf(i) = mu(i) * (1 - pii);
        Csf(i) = pii + (1 - pii) * Cs(i);
        Pf(i, :) = P(i, :) / (1 - pii);
        Pf(i, i) = 0;
    end
end

% see _kb/03-api-layer.md (me_oqn_blk -- GE-type OQN with blocking) for rationale
sigmaS = 2 ./ (Csf + 1);
muRes = c .* muf .* sigmaS;

% Fixed-point state
Ca = ones(M, 1);
Cd = Csf;
PBs = zeros(M, M);   % PB^i_j, blocking of the flow from i into j
PBe = zeros(M, 1);   % PBe_j, blocking of the external flow into j
PBh = zeros(M, M);   % PB^{h_ij}_j, blocking of the flow released by h_ij
PBa = zeros(M, 1);   % aggregate blocking probability at j
Q = zeros(M, 1);
U = zeros(M, 1);
T = zeros(M, 1);
lambda = zeros(M, 1);
Lhold = zeros(M, M); % mean occupancy of each holding node
iter = 0;

for iter = 1:options.maxiter
    Ca_old = Ca;
    PBs_old = PBs;
    PBe_old = PBe;

    % Service inflation at the blocking stations, eqs. (15)-(16): the
    % fraction PBf(i) of the service completions at i is followed by a
    % blocking period during which the server stays unavailable.
    PBf = zeros(M, 1);
    for i = 1:M
        for j = 1:M
            if Pf(i, j) > 0 && bas(j)
                PBf(i) = PBf(i) + Pf(i, j) * PBs(i, j);
            end
        end
    end
    if any(PBf >= 1 - 1e-9)
        line_error(mfilename, 'MEM transfer-blocking fixed point saturates: a station is blocked with probability one. The network has no stable operating point under BAS.');
    end
    muEff = muf .* (1 - PBf);
    CsEff = PBf + Csf .* (1 - PBf);

    % Flow balance on the carried flow, eq. (11). Under loss the fraction
    % PB of a stream is discarded; under transfer blocking every job
    % eventually enters, the delay being charged to the holding node.
    A = zeros(M, M);
    b = zeros(M, 1);
    for j = 1:M
        % see _kb/03-api-layer.md (me_oqn_blk -- GE-type OQN with blocking) for rationale
        b(j) = lambda0(j) * (1 - PBe(j));
        for i = 1:M
            if Pf(i, j) > 0
                if finiteBuf(j) && ~bas(j)
                    A(i, j) = Pf(i, j) * (1 - PBs(i, j));
                else
                    A(i, j) = Pf(i, j);
                end
            end
        end
    end
    T = (eye(M) - A') \ b;
    T(T < 0) = 0;

    % Offered (attempt) rates and the aggregate blocking probability,
    % eq. (13). A blocked job under transfer blocking re-attempts from the
    % holding node, so its stream contributes carried/(1-PB) attempts.
    attExt = zeros(M, 1);
    attInt = zeros(M, M);
    for j = 1:M
        attExt(j) = lambda0(j);
        for i = 1:M
            if Pf(i, j) > 0
                if finiteBuf(j) && bas(j)
                    attInt(i, j) = T(i) * Pf(i, j) / max(1 - PBs(i, j), 1e-12);
                else
                    attInt(i, j) = T(i) * Pf(i, j);
                end
            end
        end
    end
    lambda = attExt + sum(attInt, 1)';
    for j = 1:M
        if lambda(j) > 0
            PBa(j) = (attExt(j) * PBe(j) + sum(attInt(:, j) .* PBs(:, j))) / lambda(j);
        else
            PBa(j) = 0;
        end
    end

    % Interarrival scv of the offered flow, by GE splitting and merging.
    % A stream thinned with probability p has scv 1-p+p*Cd; the merge of
    % GE streams satisfies 1/(Cm+1) = sum_s (lam_s/lam)/(Cs+1).
    CaStreamInt = ones(M, M);
    for j = 1:M
        if lambda(j) <= 0
            continue
        end
        sum_inv = 0;
        if attExt(j) > 0
            sum_inv = sum_inv + (attExt(j) / lambda(j)) / (Ca0(j) + 1);
        end
        for i = 1:M
            if attInt(i, j) > 0
                CaStreamInt(i, j) = 1 - Pf(i, j) + Pf(i, j) * Cd(i);
                sum_inv = sum_inv + (attInt(i, j) / lambda(j)) / (CaStreamInt(i, j) + 1);
            end
        end
        if sum_inv > 0
            Ca(j) = -1 + 1 / sum_inv;
        end
    end

    % Station solution in isolation and the per-stream blocking
    % probabilities, eqs. (12) and (14), i.e. eq. (4.3) of the source
    % evaluated with each stream's own scv.
    PBe_new = zeros(M, 1);
    PBs_new = zeros(M, M);
    PBh_new = zeros(M, M);
    for j = 1:M
        if isinf(c(j))
            % Infinite server: no queueing and no blocking
            Q(j) = 0;
            if muEff(j) > 0
                Q(j) = lambda(j) / muEff(j);
            end
            U(j) = Q(j);
            Cd(j) = Ca(j);
            continue
        end
        if ~finiteBuf(j)
            % Unbounded buffer: the infinite-capacity GE building blocks
            rho = 0;
            if muEff(j) > 0
                rho = lambda(j) / (c(j) * muEff(j));
            end
            if rho >= 1
                Q(j) = Inf;
                U(j) = 1;
                Cd(j) = CsEff(j);
            elseif c(j) == 1
                Q(j) = rho * (Ca(j) + 1) / 2 + rho^2 * (CsEff(j) + Ca(j)) / (2 * (1 - rho));
                U(j) = rho;
                Cd(j) = rho^2 * CsEff(j) + (1 - rho) * Ca(j) + rho * (1 - rho);
            else
                Q(j) = me_gegec_mql(lambda(j), Ca(j), muEff(j), CsEff(j), c(j));
                U(j) = rho;
                Cd(j) = rho^2 * CsEff(j) + (1 - rho) * Ca(j) + rho * (1 - rho);
            end
            continue
        end
        [pj, Lj, Uj] = me_gegecn(lambda(j), Ca(j), muEff(j), CsEff(j), c(j), 0, N(j));
        Q(j) = Lj;
        U(j) = Uj;
        % Interdeparture scv, eq. (4), evaluated at the utilization of the
        % censored queue: with losses the offered load can exceed one
        % while the fraction of busy servers cannot.
        Cd(j) = Uj^2 * CsEff(j) + (1 - Uj) * Ca(j) + Uj * (1 - Uj);
        if attExt(j) > 0
            PBe_new(j) = me_gegecn_pb(pj, 0, N(j), c(j), CsEff(j), Ca0(j));
        end
        for i = 1:M
            if attInt(i, j) > 0
                PBs_new(i, j) = me_gegecn_pb(pj, 0, N(j), c(j), CsEff(j), CaStreamInt(i, j));
                if bas(j)
                    % see _kb/03-api-layer.md (me_oqn_blk -- GE-type OQN with blocking) for rationale
                    q = Pf(i, j) * PBs(i, j);
                    CaH = 1 - q + q * Cd(i);
                    PBh_new(i, j) = me_gegecn_pb(pj, 0, N(j), c(j), CsEff(j), CaH);
                end
            end
        end
    end

    % Relaxation on the blocking probabilities
    w = options.damping;
    PBe = (1 - w) * PBe + w * PBe_new;
    PBs = (1 - w) * PBs + w * PBs_new;
    PBh = (1 - w) * PBh + w * PBh_new;

    delta = max([max(abs(Ca - Ca_old)), max(max(abs(PBs - PBs_old))), max(abs(PBe - PBe_old))]);
    if options.verbose
        fprintf('Iteration %d: max delta = %e\n', int32(iter), delta);
    end
    if delta < options.tol
        break
    end
end

if iter == options.maxiter && delta >= options.tol
    warning('me_oqn_blk:noconverge', 'Did not converge within %d iterations (delta=%e)', int32(options.maxiter), delta);
end

% Holding node occupancy, eq. (7) and the correction of eq. (17). The jobs
% held in h_ij are physically blocked in the servers of station i, so they
% are added back to station i.
Lhold(:) = 0;
for i = 1:M
    for j = 1:M
        if bas(j) && Pf(i, j) > 0 && PBs(i, j) > 0
            rateH = T(i) * Pf(i, j) * PBs(i, j);
            muH = muRes(j) * (1 - PBh(i, j));
            if muH > 0
                Lhold(i, j) = rateH / muH;
            end
        end
    end
end
Q = Q + sum(Lhold, 2);

% see _kb/03-api-layer.md (me_oqn_blk -- GE-type OQN with blocking) for rationale
for i = 1:M
    if muf(i) > 0
        if isinf(c(i))
            U(i) = T(i) / muf(i);
        else
            U(i) = T(i) / (c(i) * muf(i));
        end
    else
        U(i) = 0;
    end
end
W = zeros(M, 1);
for i = 1:M
    if T(i) > 0
        W(i) = Q(i) / T(i);
    end
end
end
