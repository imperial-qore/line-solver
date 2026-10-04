function varargout = solver_tr_lc_analyzer(self, phase, varargin)
% VARARGOUT = SOLVER_TR_LC_ANALYZER(SELF, PHASE, VARARGIN)
%
% LOAD CONCEALMENT as a model transformation, and the first ITERATED strategy
% of TRANSFORMSOLVE.
%
% Birman and Kogan (Stochastic Models 8(3):543-563, 1992), Algorithm 2. The
% saddle point analysis of their Corollary 2 shows that chain l may be solved on
% its own provided every station is slowed by the residual capacity the other
% chains leave it,
%
%   A_i = 1 - sum_{k != l} L(i,k) * X_k,
%
% so that chain l sees the concealed demand L(i,l)/A_i. Sweeping the chains in
% Gauss-Seidel order and iterating to a fixed point is the algorithm.
%
% WHAT THIS ADDS OVER THE KERNEL. PFQN_BKLC solves each single-chain subproblem
% on a DEMAND VECTOR, with the inner solve hard-wired to pfqn_mva or pfqn_bkue.
% Here the subproblem is a real single-class Network, so the inner solve is the
% CALLER'S OWN solver: the same decomposition can be evaluated with CTMC, SSA or
% Fluid, which is what makes the concealment approximation measurable rather
% than merely asserted.
%
% IT IS NOT A STRICTLY BETTER LC, and must not be described as one. The kernel
% sees only L; ModelAdapter.aggregateChains refits the chain service law to two
% moments. On a product-form model the two coincide; off it they are different
% approximations of the same quantity.
%
% THE TOLERANCE IS FIXED AT 1e-10 AND IS NOT options.iter_tol. That is a PARITY
% requirement, not a quality knob: a looser tolerance stops the sweep at a
% different iteration in each codebase, which is how MATLAB and python came to
% report 0.99843 and 0.9835 for the same model.
%
% Phases: 'expand', 'couple' (Gauss-Seidel, one chain at a time), 'converged',
% 'lift'.
%
% See also TRANSFORMSOLVE, PFQN_BKLC, SOLVER_TR_CHAINS_ANALYZER.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

switch phase
    case 'expand'
        [varargout{1}, varargout{2}] = expand_(self, varargin{:});
    case 'couple'
        [varargout{1}, varargout{2}] = couple_(self, varargin{:});
    case 'converged'
        [varargout{1}, varargout{2}] = converged_(self, varargin{:});
    case 'lift'
        [varargout{1}, varargout{2}, varargout{3}, ...
         varargout{4}, varargout{5}, varargout{6}] = lift_(self, varargin{:});
    otherwise
        line_error(mfilename, sprintf('Unknown load-concealment phase: %s', phase));
end
end

% ---------------------------------------------------------------- expand ----

function [submodels, ctx] = expand_(self, sn, options)
% Stage 1 is the chain aggregation, which is single pass; stages 2..R+1 are the
% per-chain concealed solves, which iterate.
[chainModel, alpha, deagg] = ModelAdapter.aggregateChains(self.model);
snChain = chainModel.getStruct();
M = snChain.nstations;

[Lc, ~, ~, ~, Nc] = sn_get_demands_chain(snChain);
% The concealment is over CHAINS, and the aggregated model is supposed to carry
% one class per chain, so the demand matrix's column count and the class count
% must agree. Naming the mismatch beats indexing past the end of Lc later.
R = size(Lc, 2);
if snChain.nclasses ~= R
    line_error(mfilename, sprintf(['the chain-aggregated model has %d classes but %d chains: ' ...
        'load concealment needs one class per chain.'], snChain.nclasses, R));
end
isDelay = reshape(snChain.sched == SchedStrategy.INF, 1, []);
% The concealment slows QUEUEING stations only: a delay station holds no queue,
% so the other chains leave its capacity untouched. Zeroing the delay rows makes
% A come out as exactly 1 there rather than needing a second branch.
L = Lc;
L(isDelay, :) = 0;
Z = sum(Lc(isDelay, :), 1);
N = reshape(Nc, 1, []);

% One single-class model per chain, built ONCE and re-concealed in place on
% every sweep. removeClass looks the class up BY NAME on its copy, so the
% original model's class objects stay valid as the copy shrinks.
origClasses = chainModel.classes;
submodels = cell(1, R);
baseService = cell(1, R);
for l = 1:R
    m = chainModel;
    for k = 1:R
        if k ~= l
            m = ModelAdapter.removeClass(m, origClasses{k});
        end
    end
    submodels{l} = m;
    baseService{l} = cell(1, M);
    for i = 1:M
        % ONLY the queueing stations are ever concealed, so only their service
        % law is retained, and it is read from the CHAIN MODEL, where the class
        % objects still name the chain they belong to.
        if isDelay(i)
            continue
        end
        baseService{l}{i} = chainModel.stations{i}.getService(origClasses{l});
    end
end

% ONE CLASS PER CHAIN is the ordinary case for load concealment, and there the
% chain aggregation is the IDENTITY: aggregateChains returns no deaggregation
% tables, because the chain answer already IS the class answer. The lift then
% just re-indexes chains onto their classes.
ctx = struct('sn_orig', sn, 'alpha', alpha, 'deagg', deagg, ...
    'L', L, 'N', N, 'Z', Z, 'M', M, 'R', R, ...
    'identity', sn.nchains >= sn.nclasses, 'Korig', sn.nclasses, ...
    'isDelay', isDelay, 'baseService', {baseService}, 'origClasses', {origClasses}, ...
    'X', lc_seed(L, N, Z), 'Xold', [], 'tol', 1e-10, 'iterated', true);
ctx.Xold = -ones(1, R);
submodels = conceal_all(ctx, submodels);

line_debug(options, '%s: load concealment, %d chains over %d stations, tol %g.', ...
    self.getName(), R, M, ctx.tol);
end

function X = lc_seed(L, N, Z)
% Step 1 of Algorithm 2: the saddle point utilizations of Corollary 1 seed the
% sweep, with the same fallbacks and the same capacity clamp PFQN_BKLC applies,
% so the elevated path starts the iteration from the same point as the kernel.
R = size(L, 2);
X = zeros(1, R);
try
    [~, ~, X] = pfqn_bk(L, N, Z);
catch
    X = zeros(1, R);
end
X = reshape(X, 1, []);
ok = isfinite(X) & X >= 0;
X(~ok) = 0;
for r = 1:R
    if X(r) == 0 && N(r) > 0
        X(r) = N(r) / (Z(r) + sum(L(:, r)));
    end
end
% A chain cannot draw more than the capacity of its own slowest station
for r = 1:R
    cap = max(L(:, r));
    if cap > 0
        X(r) = min(X(r), 1 / cap);
    end
end
end

% ---------------------------------------------------------------- couple ----

function [submodels, ctx] = couple_(~, ctx, submodels, res, e)
% GAUSS-SEIDEL. Chain e's throughput is published the moment it is known, and
% every remaining subproblem is re-concealed against it, so chain e+1 of this
% same sweep already sees it. Deferring this to the end of the sweep would be
% Jacobi: the same fixed point, reached in a different number of sweeps, which
% is exactly what the cross-codebase sweep counts are pinned on.
ctx.X(e) = chain_tput(res{e});
submodels = conceal_all(ctx, submodels);
end

function x = chain_tput(r)
% The single-class subproblem's system throughput.
x = r.XN;
x = x(:)';
if isempty(x)
    x = 0;
else
    x = x(1);
end
if ~isfinite(x) || x < 0
    x = 0;
end
end

function submodels = conceal_all(ctx, submodels)
for l = 1:ctx.R
    A = 1 - (ctx.L * ctx.X(:) - ctx.L(:, l) * ctx.X(l));
    A = max(A, GlobalConstants.FineTol);
    m = submodels{l};
    for i = 1:ctx.M
        if ctx.isDelay(i)
            continue
        end
        % dist_scale_rate multiplies the RATE by its factor, so a factor of A_i
        % divides the mean by A_i: exactly the concealed demand L(i,l)/A_i.
        m.stations{i}.setService(m.classes{1}, dist_scale_rate(ctx.baseService{l}{i}, A(i)));
    end
    % A full rebuild, not refreshProcesses: the submodel comes from
    % ModelAdapter.removeClass and its struct is not yet complete enough for the
    % incremental process refresh (which reads sn.nodetype).
    m.refreshStruct(true);
    submodels{l} = m;
end
end

% ------------------------------------------------------------- converged ----

function [done, ctx] = converged_(~, ctx, res, it) %#ok<INUSD>
done = max(abs(ctx.X - ctx.Xold)) <= ctx.tol * max(1, max(abs(ctx.X)));
ctx.Xold = ctx.X;
end

% ------------------------------------------------------------------ lift ----

function [QN,UN,RN,TN,CN,XN] = lift_(~, ctx, res)
% Reassemble the per-chain answers into a chain-level table, then deaggregate it
% exactly as SOLVER_TR_CHAINS_ANALYZER does: the concealment changed how the
% chain table was obtained, not what it means.
M = ctx.M;
R = ctx.R;
Qchain = zeros(M, R);
Uchain = zeros(M, R);
Rchain = zeros(M, R);
Tchain = zeros(M, R);
Xchain = zeros(1, R);
for l = 1:R
    Qchain(:, l) = column_of(res{l}.QN, M);
    Uchain(:, l) = column_of(res{l}.UN, M);
    Rchain(:, l) = column_of(res{l}.RN, M);
    Tchain(:, l) = column_of(res{l}.TN, M);
    Xchain(l) = chain_tput(res{l});
end
if ctx.identity
    K = ctx.Korig;
    QN = zeros(M,K); UN = zeros(M,K); RN = zeros(M,K); TN = zeros(M,K); XN = zeros(1,K);
    for c = 1:R
        inchain = ctx.sn_orig.inchain{c};
        k = inchain(1);
        QN(:,k) = Qchain(:,c);
        UN(:,k) = Uchain(:,c);
        RN(:,k) = Rchain(:,c);
        TN(:,k) = Tchain(:,c);
        XN(k) = Xchain(c);
    end
    CN = sum(RN, 1);
    return
end

[QN,UN,RN,TN,CN,XN] = sn_deaggregate_chain_results(ctx.sn_orig, ctx.deagg.Lchain, [], ...
    ctx.deagg.STchain, ctx.deagg.Vchain, ctx.alpha, Qchain, Uchain, Rchain, Tchain, [], Xchain);

if isempty(CN)
    CN = sum(RN, 1);
end
end

function v = column_of(A, M)
v = zeros(M, 1);
if isempty(A)
    return
end
n = min(M, size(A, 1));
v(1:n) = A(1:n, 1);
end
