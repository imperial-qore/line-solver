function [QN,UN,RN,TN,CN,XN,runtime] = solver_ba_qrf_analyzer(sn, options)
% SOLVER_BA_QRF_ANALYZER Adapter for QRF library functions within SolverBA
%
% [QN,UN,RN,TN,CN,XN,RUNTIME] = SOLVER_BA_QRF_ANALYZER(SN, OPTIONS)
%
% Bridges the LINE sn struct to QRF (Quadratic Reduction Framework) library
% functions for approximating performance metrics of single-class closed
% queueing networks with PH service.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Tstart = tic;

M = sn.nstations;
K = sn.nclasses;
N = sum(sn.njobs);
S = sn.nservers;
PH = sn.proc;

% THE PREMISES ARE BA_METHOD_REFUSAL'S, the predicate SolverBA.supportsModelMethod
% reports: single-class closed; the alpha-free arms refuse a delay, a
% multiserver or a load-dependent station by naming the two arms that serve
% it (a c>1 station solved as c=1 is not a bound in either direction, measured
% +200% at c=3,N=1 and -10% at c=3,N=3); the load-dependent arms carry
% alpha(i,n) -- the rate law of a delay (n), of a c-server station (min(n,c))
% and of limited load dependence -- and sn_to_qrf_alpha owns the single
% restriction that survives there, exponential service where a station serves
% several jobs at once; a multi-phase law must serve FCFS, the semantics of the
% one phase per station the local state carries; and the blocking arms need
% the tables sn_to_qrf_blocking derives. see _kb/06-solver-catalog.md
refusal = ba_method_refusal(sn, options.method, options);
if ~isempty(refusal)
    line_error(mfilename, '%s', refusal);
end
[alpha_sn, ~, ~, qrfPeak] = sn_to_qrf_alpha(sn);

% Extract MAPs as cell array: MAPs{i} = {D0, D1}
MAPs = cell(M, 1);
K_phases = zeros(M, 1);
for i = 1:M
    if ~isempty(PH{i}{1})
        MAPs{i} = PH{i}{1};
        K_phases(i) = size(MAPs{i}{1}, 1);
    else
        K_phases(i) = 1;
        MAPs{i} = {-1, 1}; % fallback exponential rate 1
    end
end

% Build routing matrix (M x M) from sn.rt (MK x MK)
rt = zeros(M);
for i = 1:M
    for j = 1:M
        rt(i,j) = sn.rt(i, j); % K=1, so indexing is direct
    end
end

% Extract mu and v arrays from MAPs
Kmax = max(K_phases);
mu = zeros(M, Kmax, Kmax);
v = zeros(M, Kmax, Kmax);
for i = 1:M
    D0 = MAPs{i}{1};
    D1 = MAPs{i}{2};
    for h = 1:K_phases(i)
        for k = 1:K_phases(i)
            % Both mu and v are indexed (from phase, to phase), matching
            % qrf_noblo_*'s q{i,j}(k,h) = v{i}(k,h) + r(i,i)*mu{i}(k,h).
            % mu{i}(k,h) is the completion rate from phase k leaving to phase
            % h, i.e. D1(k,h); v{i}(k,h) is the background phase change k -> h,
            % i.e. D0(k,h) off the diagonal. Writing v as D0(h,k) transposes
            % it, which is invisible for a reversible D0 but reverses the phase
            % order of an Erlang (upper bidiagonal D0) and silently changes the
            % model.
            mu(i, h, k) = D1(h, k);
            if h == k
                v(i, h, k) = 0;
            else
                v(i, h, k) = D0(h, k);
            end
        end
    end
end

% Dispatch based on method
switch options.method
    case 'qrf.mmi'
        MR = 1;
        [UN_qrf, QN_qrf] = qrf_noblo_mmi(M, MR, K_phases(:)', N, mu, v, rt);

    case 'qrf.mem'
        [UN_qrf, QN_qrf] = qrf_noblo_mem(MAPs, N, rt);

    case 'qrf.bethe'
        % Same polytope and same phase-1 start as 'qrf.mmi'; the objective is
        % the tree-reweighted (Bethe) free entropy at the uniform spanning-tree
        % weight lambda = 1/M, the largest uniform weight at which the program
        % is convex. See qrf_noblo_bethe.m.
        MR = 1;
        [UN_qrf, QN_qrf] = qrf_noblo_bethe(M, MR, K_phases(:)', N, mu, v, rt);

    case 'qrf.mmi.ld'
        [UN_qrf, QN_qrf, BN_qrf] = qrf_noblo_mmi_ld(MAPs, N, rt, qrf_alpha_of(options, alpha_sn));

    case 'qrf.mmi.linear'
        [UN_qrf, QN_qrf, BN_qrf] = qrf_noblo_mmi_linear(MAPs, N, rt, qrf_alpha_of(options, alpha_sn));

    case 'qrf.bas.mmi'
        % SERVED AGAIN SINCE 2026-09-03. It was removed on 2026-08-01 with
        % qrf_bas_mmi_simple.m, a reduced model whose signature carried no
        % MM/MM1/ZZ/ZM, so it dropped the arrival terms of THM30/THM3 and
        % refused MR > 1. All three reasons the refusal outlived that routine
        % are now gone: this dispatch derives the four tables exactly as
        % 'qrf.bas.mem' and 'qrf.bas.bethe' do; qrf_bas_mmi.m runs, its stale
        % duplicate THM30-group families having been deleted; and the ports'
        % MI block now spans 0..F as this one and the AMPL source do, so all
        % four codebases minimise the SAME functional over the SAME polytope.
        %
        % THAT IS FORMULATION PARITY, NOT VALUE PARITY, AND IT CANNOT BE MORE.
        % The MI objective is CONCAVE on this polytope, so its minimisers are
        % VERTICES and every codebase settles on whichever one its own solver
        % reaches: fmincon takes interior SQP steps here, while the python and
        % C++ ports probe vertices explicitly. Measured on the three-station
        % feeder model at N=3, MATLAB reports U1 = 0.5434 against python's
        % 0.5090 on the same polytope. Use 'qrf.bas.bethe' when the answer has
        % to be a property of the model rather than of the solver -- that is
        % the whole point of the tree-reweighted weight lambda = 1/M.
        if isfield(options.config, 'qrf_params')
            qp = options.config.qrf_params;
        else
            qp = [];
        end
        if isempty(qp)
            [blk, blkMsg] = sn_to_qrf_blocking(sn, options);
            if ~isempty(blkMsg)
                line_error(mfilename, ['qrf.bas.mmi cannot be applied to this model: %s ' ...
                    'Supply options.config.qrf_params explicitly to override the derivation.'], blkMsg);
            end
            qp = blk;
        end
        [UN_qrf, QN_qrf] = qrf_bas_mmi(qp.f, M, qp.MR, qp.MM, qp.MM1, qp.ZZ, qp.ZM, ...
            qp.BB, K_phases(:)', qp.F, N, mu, v, rt);

    case 'qrf.bas.mem'
        if isfield(options.config, 'qrf_params')
            qp = options.config.qrf_params;
        else
            qp = [];
        end
        if isempty(qp)
            % Same tables as 'qrf.bas' over the same BAS polytope, so the same
            % derivation serves both rather than only the LP arm.
            [blk, blkMsg] = sn_to_qrf_blocking(sn, options);
            if ~isempty(blkMsg)
                line_error(mfilename, ['qrf.bas.mem cannot be applied to this model: %s ' ...
                    'Supply options.config.qrf_params explicitly to override the derivation.'], blkMsg);
            end
            qp = blk;
        end
        f = qp.f;
        MR = qp.MR;
        MM = qp.MM;
        MM1 = qp.MM1;
        ZZ = qp.ZZ;
        ZM = qp.ZM;
        BB = qp.BB;
        F = qp.F;
        [UN_qrf, QN_qrf] = qrf_bas_mem(f, M, MR, MM, MM1, ZZ, ZM, BB, K_phases(:)', F, N, mu, v, rt);

    case 'qrf.bas.bethe'
        % Same tables and the same BAS polytope as 'qrf.bas' and 'qrf.bas.mem',
        % so the same derivation serves all three. The objective is the
        % tree-reweighted free entropy of 'qrf.bethe' evaluated over the BAS
        % decision vector; see qrf_bas_bethe.m for the weight lambda = 1/M.
        % All three BAS objectives share one polytope since the stale duplicate
        % THM30-group families were deleted from qrf_bas_mmi.m on 2026-09-03.
        if isfield(options.config, 'qrf_params')
            qp = options.config.qrf_params;
        else
            qp = [];
        end
        if isempty(qp)
            [blk, blkMsg] = sn_to_qrf_blocking(sn, options);
            if ~isempty(blkMsg)
                line_error(mfilename, ['qrf.bas.bethe cannot be applied to this model: %s ' ...
                    'Supply options.config.qrf_params explicitly to override the derivation.'], blkMsg);
            end
            qp = blk;
        end
        [UN_qrf, QN_qrf] = qrf_bas_bethe(qp.f, M, qp.MR, qp.MM, qp.MM1, qp.ZZ, qp.ZM, ...
            qp.BB, K_phases(:)', qp.F, N, mu, v, rt);

    case 'qrf.bas'
        params = sn_to_qrf_params(sn, MAPs, K_phases, N, mu, v, rt, options, true);
        params.verbose = qrf_lp_is_debug(options);
        % Objective queue 1, sense 'max': the analyzer reports an UPPER bound
        % on utilization, the direction the relaxation guarantees. qrf_bas
        % DEFAULTS to 'U1min', whose vertex has U = 0 at the objective station
        % and reads out as a zero throughput. python and the JAR use 'max' too.
        [result] = qrf_bas(params, 1, 'max');
        % qrf_bas returns utilization bounds; derive QN via visit ratios
        [UN_qrf, QN_qrf] = derive_qn_from_bounds(result.U(:)', M, N, S, PH, sn);

    case 'qrf.rsrd'
        % RS-RD carries NO blocking tables: qrf_rsrd reads only M, N, F, K, mu,
        % v, r and the optional alpha, and its PBB constraint sums over every
        % queue that can be full. So it needs no enumeration, and unlike
        % 'qrf.bas' it admits SEVERAL finite-capacity queues.
        params = sn_to_qrf_params(sn, MAPs, K_phases, N, mu, v, rt, options, false);
        params.verbose = qrf_lp_is_debug(options);
        [result] = qrf_rsrd(params, 1, 'max');   % same direction as qrf.bas above
        [UN_qrf, QN_qrf] = derive_qn_from_bounds(result.U(:)', M, N, S, PH, sn);

    otherwise
        line_error(mfilename, 'Unknown QRF method: %s', options.method);
end

% The alpha-free arms return no BN because their alpha is identically 1, and
% there BN = P(n >= 1) = UN_qrf: a single server's departure rate is
% proportional to the probability that it is busy. Setting it here rather than
% in each arm keeps the readout below one formula.
if ~exist('BN_qrf', 'var') || isempty(BN_qrf)
    BN_qrf = UN_qrf;
end

% Normalize QN to population constraint
if sum(QN_qrf) > 0
    QN_qrf = QN_qrf / sum(QN_qrf) * N;
end

% Map 1D QRF results to M x K matrices (K=1)
UN = zeros(M, K);
QN = zeros(M, K);
TN = zeros(M, K);
RN = zeros(M, K);
XN = zeros(1, K);
CN = zeros(1, K);

QN(:, 1) = QN_qrf(:);

% Derive all metrics from the QRF utilization using Little's law and visit
% ratios.
%
% UN_qrf(i) is sum_{ni>=1, ki} p2(i,ni,ki,i,ni,ki), i.e. P(n_i >= 1)
% marginalised over phase. That IS the utilization and it is bounded by 1
% through the ONE constraint. An earlier comment here held that it "is the sum
% of phase effective utilizations, which can exceed 1 for phase-type service"; that
% is false, and what exceeds 1 is the sum ACROSS stations, which is expected.
% Adjudicated against glpsol on the paper's AMPL model at K=2 (2-station
% closed, Exp + Erlang-2): every UN_qrf entry lands inside the glpsol
% [min,max] range for that station's utilization.
%
% Throughput therefore comes from a finite-server station via
% U_i = X * V_i * stime_i, which is exact for a single server. The previous
% inversion X = QN(refstat) / (V*stime_ref) instead assumed R == stime at the
% reference station, which holds only for an infinite server; with a FCFS
% reference station it returned X above the bottleneck capacity and hence
% UN > 1. UN_qrf(i)/(V_i*stime_i) is consistent across stations, so any
% finite-server station determines X.

% Compute visit ratios
if isfield(sn, 'visits') && ~isempty(sn.visits)
    V = cellsum(sn.visits);
else
    [visits] = sn_refresh_visits(sn, sn.chains, sn.rt, sn.rtnodes);
    V = cellsum(visits);
end

refstat = sn.refstat(1);

% Per-station mean service times
stimes = zeros(1, M);
for i = 1:M
    if ~isempty(PH{i}{1})
        stimes(i) = map_mean(PH{i}{1});
    end
end

% System throughput from the ALPHA-WEIGHTED marginal mean BN, the mean number
% of jobs actually in service: E[min(n,c)] at a c-server station, E[n] at a
% delay, P(n >= 1) at a single server. That is what the departure rate is
% proportional to, since alpha(i,n) scales the completion rate, so
% T_i = BN_i / stime_i holds exactly at the relaxed point and XN = T_i / V_i.
%
% The single-server case is the former UN_qrf(i)*S(i)/(V*stime) unchanged, S
% being 1 and BN being UN_qrf there. A delay no longer needs the separate
% refstat fallback below, which computed QN(refstat)/(V*stime): alpha = n makes
% BN = E[n] = QN, so the general formula already IS that fallback.
UN_qrf = UN_qrf(:)';
BN_qrf = BN_qrf(:)';
for i = 1:M
    if stimes(i) > 0 && V(i, 1) > 0 && BN_qrf(i) > 0
        XN(1) = BN_qrf(i) / (V(i, 1) * stimes(i));
        break;
    end
end
if XN(1) == 0
    % No finite-server station carries load: fall back to the reference
    % station, where QN = XN * V * stime holds exactly for a delay.
    if stimes(refstat) > 0 && V(refstat, 1) > 0
        XN(1) = QN(refstat, 1) / (V(refstat, 1) * stimes(refstat));
    end
end

% Derive per-station metrics from XN and visit ratios
for i = 1:M
    if ~isempty(PH{i}{1})
        stime = stimes(i);
        TN(i, 1) = XN(1) * V(i, 1);
        if isinf(S(i))
            % Delay (infinite server): UN = QN by LINE convention
            UN(i, 1) = QN(i, 1);
        else
            % Finite server: the busy fraction of the station's DECLARED peak
            % capacity, BN being the mean number of jobs in service. The
            % normalizer is nservers times the reachable lld peak, LINE's one
            % U = T*S/peak convention (see sn_to_qrf_alpha); at peak 1 it is
            % UN_qrf(i), the utilization the alpha-free arms report directly.
            if stime > 0
                UN(i, 1) = BN_qrf(i) / qrfPeak(i);
            end
        end
        % Response time via Little's law
        if TN(i, 1) > 0
            RN(i, 1) = QN(i, 1) / TN(i, 1);
        end
    end
end

% Cycle time
if XN(1) > 0
    CN(1) = N / XN(1);
end

% Clean NaN values
QN(isnan(QN)) = 0;
UN(isnan(UN)) = 0;
RN(isnan(RN)) = 0;
TN(isnan(TN)) = 0;
XN(isnan(XN)) = 0;
CN(isnan(CN)) = 0;

runtime = toc(Tstart);
end

function alpha = qrf_alpha_of(options, alpha_sn)
% ALPHA = QRF_ALPHA_OF(OPTIONS, ALPHA_SN)
%
% The load-dependent scaling the ld arms run on: derived from the model by
% sn_to_qrf_alpha, with options.config.qrf_alpha kept as an explicit override,
% exactly as options.config.qrf_params overrides the derived blocking tables.
% Reaching for isfield rather than the field directly because reading
% options.config.qrf_alpha throws "Unrecognized field name" whenever config
% carries no such field, which is the common case.
alpha = alpha_sn;
if isfield(options.config, 'qrf_alpha') && ~isempty(options.config.qrf_alpha)
    alpha = options.config.qrf_alpha;
end
end

function bool = qrf_lp_is_debug(options)
% BOOL = QRF_LP_IS_DEBUG(OPTIONS)
%
% Whether the QRF linear programs may print. They are LOUD when they do: the
% constraint builder narrates every block it emits (ZERO, ONE, SYMMETRY,
% MARGINALS, UEFF, THM1..THM4), then reports the variable count, then linprog
% itself is put in 'Display','final' and adds its own exit message -- on
% cqn_bas_blocking that is some thirty lines for a two-station model.
%
% That belongs at DEBUG and nowhere else. The old test was `options.verbose > 0`,
% which is TRUE at VerboseLevel.STD (= 1), so an ordinary getAvgTable() printed
% the whole LP trace. The sibling paths already had it right and are the
% precedent: solver_ba_analyzer hardcodes params.verbose = false for the lr/qr
% bounds, the JAR sets settings.verbose = false, and the python and C++ twins
% never pass a verbose flag at all.
%
% GlobalConstants.Verbose is consulted too, matching solver_ctmc: the global
% level can be raised to DEBUG without it reaching this options struct.
%
% Note QRF_BAS and QRF_RSRD default their own `verbose` field to TRUE when the
% caller sets none, so this must be passed explicitly rather than omitted.
bool = false;
if isfield(options, 'verbose') && options.verbose >= VerboseLevel.DEBUG
    bool = true;
end
if GlobalConstants.Verbose >= VerboseLevel.DEBUG
    bool = true;
end
end

function params = sn_to_qrf_params(sn, MAPs, K_phases, N, mu, v, rt, options, needBlocking)
% Build params struct for qrf_bas / qrf_rsrd from sn struct.
%
% NEEDBLOCKING is true only for 'qrf.bas'. The blocking configuration tables
% (f, MR, BB, MM, MM1, ZZ) are DERIVED from sn by sn_to_qrf_blocking rather
% than demanded from the caller: the model already fixes every one of them.
% options.config.qrf_params remains an explicit override, so hand-built tables
% and the goldens produced with them keep working unchanged.
M = sn.nstations;
params = struct();
params.M = M;
params.N = N;
params.K = K_phases(:);
params.r = rt;

% Convert mu/v from 3D arrays to cell arrays expected by qrf_bas/qrf_rsrd
% Each mu{i}, v{i} must be a [Ki x Ki] matrix (even for Ki=1)
params.mu = cell(M, 1);
params.v = cell(M, 1);
for i = 1:M
    Ki = K_phases(i);
    params.mu{i} = reshape(mu(i, 1:Ki, 1:Ki), Ki, Ki);
    params.v{i} = reshape(v(i, 1:Ki, 1:Ki), Ki, Ki);
end

% Capacity. Read through sn_to_qrf_capacity, which decides BINDING with
% sn_get_buffer_size rather than testing sn.cap for finiteness: refreshCapacity
% derives a finite classcap at every station of every closed model, so raw
% sn.cap reports a buffer even where none can ever refuse a job. It also folds
% setClassCapacity and the reachable population in, which sn.cap alone does not.
[Fcap, ~, capMsg] = sn_to_qrf_capacity(sn);
if ~isempty(capMsg)
    line_error(mfilename, 'The ''%s'' method cannot be applied: %s', options.method, capMsg);
end
params.F = Fcap;

% Blocking tables. Derived from sn by sn_to_qrf_blocking -- the model fixes
% every one of them -- with options.config.qrf_params kept as an explicit
% override for hand-built tables. Reaching for isfield rather than the field
% directly because reading options.config.qrf_params throws "Unrecognized field
% name" whenever config carries no such field, which is the common case.
if isfield(options.config, 'qrf_params')
    qp = options.config.qrf_params;
else
    qp = [];
end

if needBlocking && isempty(qp)
    % Derive. This replaces a refusal: the earlier code demanded the tables
    % from the caller because the fallback it had removed substituted
    % no-blocking defaults (MR=1), which measured against SolverCTMC on
    % sanity_CQN_rm_{fcfs,ps}_1class sat 4.16667 from exact where properly
    % parameterised tables sit at 0.133333, i.e. 31x closer. Refusing was right
    % as long as the tables had to be invented; they do not, so they are built.
    [blk, blkMsg] = sn_to_qrf_blocking(sn, options);
    if ~isempty(blkMsg)
        line_error(mfilename, ['The ''%s'' method cannot be applied to this model: %s ' ...
            'Supply options.config.qrf_params explicitly to override the derivation.'], ...
            options.method, blkMsg);
    end
    params.f = blk.f;
    params.F = blk.F;
    params.MR = blk.MR;
    params.BB = blk.BB;
    params.MM = blk.MM;
    params.ZZ = blk.ZZ;
    params.ZM = blk.ZM;
    params.MM1 = blk.MM1;
    line_debug(options, sprintf(['QRF: derived BAS blocking tables -- f=%d, MR=%d, ZM=%d, ' ...
        'blockers=[%s]'], blk.f, blk.MR, blk.ZM, num2str(blk.blockers)));
elseif ~isempty(qp)
    % Explicit override. Validated only for 'qrf.bas': 'qrf.rsrd' reads none of
    % these fields, so demanding them there would reject a well-formed call.
    if needBlocking
        required_bas = {'f', 'MR', 'BB', 'MM', 'MM1', 'ZZ', 'ZM'};
        for idx = 1:length(required_bas)
            if ~isfield(qp, required_bas{idx})
                line_error(mfilename, 'qrf_params must contain field ''%s'' for qrf.bas method.', required_bas{idx});
            end
        end
    end
    if isfield(qp, 'f'), params.f = qp.f; end
    if isfield(qp, 'MR'), params.MR = qp.MR; end
    if isfield(qp, 'BB'), params.BB = qp.BB; end
    if isfield(qp, 'MM'), params.MM = qp.MM; end
    if isfield(qp, 'ZZ'), params.ZZ = qp.ZZ; end
    if isfield(qp, 'ZM'), params.ZM = qp.ZM; end
    if isfield(qp, 'MM1'), params.MM1 = qp.MM1; end
    if isfield(qp, 'F'), params.F = qp.F(:); end
end

% Load-dependent alpha
if isfield(options.config, 'qrf_alpha') && ~isempty(options.config.qrf_alpha)
    alpha_mat = options.config.qrf_alpha;
    params.alpha = cell(M, 1);
    for i = 1:M
        params.alpha{i} = alpha_mat(i, :)';
    end
end
end

function [UN_qrf, QN_qrf] = derive_qn_from_bounds(U_bounds, M, N, S, PH, sn)
% Derive QN from utilization bounds using visit ratios and Little's law.
% For bounds methods (qrf.bas, qrf.rsrd), only utilization is returned.
% We compute system throughput from the utilization at queue stations,
% then derive QN at all stations via Little's law.
UN_qrf = U_bounds;

% Compute visit ratios
if isfield(sn, 'visits') && ~isempty(sn.visits)
    V = cellsum(sn.visits);
else
    [visits] = sn_refresh_visits(sn, sn.chains, sn.rt, sn.rtnodes);
    V = cellsum(visits);
end

% Find XN from the first finite-server station with nonzero utilization
XN_est = 0;
for i = 1:M
    if ~isinf(S(i)) && ~isempty(PH{i}{1}) && U_bounds(i) > 0 && V(i,1) > 0
        stime = map_mean(PH{i}{1});
        if stime > 0
            % UN = XN * V(i) * stime / S(i), so XN = UN * S(i) / (V(i) * stime)
            XN_est = U_bounds(i) * S(i) / (V(i, 1) * stime);
            break;
        end
    end
end

% Derive QN at each station: QN(i) = XN * V(i) * RN(i)
% For single-server queues: RN(i) >= stime(i) (at least one service time)
% Use Little's law: QN(i) = TN(i) * RN(i) where TN(i) = XN * V(i)
QN_qrf = zeros(1, M);
for i = 1:M
    if ~isempty(PH{i}{1})
        stime = map_mean(PH{i}{1});
        TN_i = XN_est * V(i, 1);
        if isinf(S(i))
            % Delay: QN = TN * stime
            QN_qrf(i) = TN_i * stime;
            UN_qrf(i) = QN_qrf(i); % LINE convention for INF server
        else
            % Queue: QN = TN * stime / (1 - U) for M/G/1-like estimate
            if U_bounds(i) < 1
                QN_qrf(i) = TN_i * stime / (1 - U_bounds(i));
            else
                QN_qrf(i) = N; % saturated
            end
        end
    end
end

% Rescale to population constraint
if sum(QN_qrf) > 0
    QN_qrf = QN_qrf / sum(QN_qrf) * N;
end
end
