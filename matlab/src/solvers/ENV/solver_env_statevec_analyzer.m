function varargout = solver_env_statevec_analyzer(self, phase, it, e)
% SOLVER_ENV_STATEVEC_ANALYZER  Full state-vector coupling for SolverENV.
%
% Alternative to the mean-field analyzer (solver_env_meanfield_analyzer). Instead
% of collapsing each stage to marginal mean queue lengths and re-seeding the
% next stage with initFromMarginal, this analyzer carries the entire state
% probability vector across environment switches and propagates it with the
% CTMC transient (matrix-exponential action) on the per-stage generator.
%
% For each stage e with infinitesimal generator Q_e (supplied by a SolverCTMC
% inner solver) and entry distribution pi_enter{e}, one iteration computes the
% transient pi(t) = pi_enter{e} * exp(Q_e t) over the stage time span, then:
%   - the exit distribution toward each destination h, pi_exit{e}{h}, as the
%     expectation of pi(t) at the (random) e->h transition time, weighted by the
%     increments of the e->h transition CDF proc{e}{h};
%   - the sojourn-end distribution pi_timeavg{e}, weighted by the overall stage
%     holding-time CDF holdTime{e}, used in the environment-averaged blend.
% Entry distributions are chained as
%   pi_enter{e} = sum_h probOrig(h,e) * resetStateFun{h,e}(pi_exit{h}{e}),
% renormalised, and iterated to an L1 fixed point. The blend reuses the
% discipline-aware CTMC marginal mapping (solver_ctmc_avg_from_pi).
%
% This weighting scheme is identical to solver_env_meanfield_analyzer; the sole
% difference is that the full joint distribution is propagated and chained
% rather than its marginal means, so the two agree when the marginal collapse
% is exact and differ when inter-class/inter-station correlations matter.
%
% Backend: requires a SolverCTMC inner solver (an explicit enumerated generator
% and state space). MAM/LDQBD backends are not yet supported.
%
% Phase dispatch (called by the SolverENV EnsembleSolver hooks):
%   'pre'       pre_(self,it)            -> []
%   'analyze'   analyze_(self,it,e)      -> [results_e, runtime]
%   'post'      post_(self,it)           -> []
%   'finish'    finish_(self)            -> []
%   'converged' converged_(self,it)      -> bool
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

switch phase
    case 'pre'
        pre_(self, it);
    case 'analyze'
        [varargout{1}, varargout{2}] = analyze_(self, it, e);
    case 'post'
        post_(self, it);
    case 'finish'
        finish_(self);
    case 'converged'
        varargout{1} = converged_(self, it);
    otherwise
        line_error(mfilename, sprintf('Unknown statevec-analyzer phase: %s', phase));
end
end

% -------------------------------------------------------------------------
function pre_(self, it)
% Build the per-stage generator, state space and metric-mapping data once,
% and initialise the per-stage entry distributions.
if it ~= 1
    return
end
E = self.getNumberOfModels;
self.Qgen        = cell(1, E);
self.SS          = cell(1, E);
self.SSaggr      = cell(1, E);
self.statevecData = cell(1, E);
self.piEnter     = cell(1, E);
self.piExitDest  = cell(1, E);
self.piTimeAvg   = cell(1, E);

for e = 1:E
    solver_e = self.solvers{e};
    opts_e = solver_e.getOptions;
    if ~isfield(opts_e, 'timespan') || ~isfinite(opts_e.timespan(2))
        line_error(mfilename, sprintf(['The statevec analyzer requires a finite inner-solver timespan ' ...
            'for stage %d, e.g. CTMC(model,''timespan'',[0,T]).'], e));
    end
    if isa(solver_e, 'SolverCTMC')
        % CTMC backend: explicit enumerated generator + state space.
        [Q, SSp, SSaggr, ~, arvRates, depRates, sn_e] = solver_ctmc(self.sn{e}, opts_e);
        self.SS{e}     = SSp;
        self.SSaggr{e} = SSaggr;
        self.statevecData{e} = struct('backend', 'ctmc', 'arvRates', arvRates, ...
            'depRates', depRates, 'sn', sn_e, 'options', opts_e);
    elseif isa(solver_e, 'SolverMAM')
        % MAM backend: level-dependent QBD blocks flattened to a generator.
        % see _kb/06-solver-catalog.md for rationale
        [~, ~, ~, ~, ~, ~, ~, ld] = solver_mam_ldqbd(self.sn{e}, opts_e);
        [Q, levelOf] = solver_mam_ldqbd_flatten(ld);
        self.SS{e}     = [];
        self.SSaggr{e} = [];
        self.statevecData{e} = struct('backend', 'mam', 'ld', ld, ...
            'levelOf', levelOf, 'options', opts_e);
    else
        line_error(mfilename, sprintf(['The statevec analyzer requires a SolverCTMC or SolverMAM ' ...
            'inner solver, but environment stage %d uses %s.'], e, class(solver_e)));
    end
    self.Qgen{e} = Q;
    % Warm-start each entry distribution from the stage's own stationary
    % distribution (a valid probability vector over its state space).
    pi0 = ctmc_solve_reducible(Q);
    pi0 = pi0(:)';
    pi0(pi0 < 0) = 0;
    if sum(pi0) > 0
        pi0 = pi0 / sum(pi0);
    end
    self.piEnter{e} = pi0;
end
self.piEnterPrev = self.piEnter;
end

% -------------------------------------------------------------------------
function [results_e, runtime] = analyze_(self, it, e)
% Propagate the entry distribution of stage e through its sojourn and store
% the per-destination exit distributions and the sojourn-end distribution.
results_e = struct();
results_e.statevec = struct('ok', false);
T0 = tic;

E = self.getNumberOfModels;
Q  = self.Qgen{e};
data = self.statevecData{e};
opts = data.options;

pi0 = self.piEnter{e}(:)';
t0 = opts.timespan(1);
t1 = opts.timespan(2);

% Deterministic sojourn option (off by default): exit is pi0*exp(Q*d_e), same
% toward every destination. see _kb/06-solver-catalog.md for rationale
if isfield(self.options,'sojourn') && strcmpi(self.options.sojourn,'deterministic')
    d_e = max(map_mean(self.envObj.holdTime{e}), eps);
    % Exact deterministic sojourn via uniformization: the exit is pi0*exp(Q*d_e)
    % and the blend uses the time-average (1/d_e) * \int_0^{d_e} exp(Q t) dt.
    [piAvg, piEx] = ctmc_timeaverage(pi0, Q, d_e);
    piExit_e = cell(1, E);
    for h = 1:E
        if self.E0(e, h) > 0
            piExit_e{h} = piEx;
        else
            piExit_e{h} = [];
        end
    end
    self.piExitDest{e} = piExit_e;
    self.piTimeAvg{e}  = piAvg;
    results_e.statevec.ok = true;
    runtime = toc(T0);
    return
end

% Exponential environment sojourn: exit equals the time-average via the
% resolvent s*pi*(sI-Q)^{-1}. see _kb/06-solver-catalog.md for rationale
expSojourn = true;
for h = 1:E
    if self.E0(e, h) > 0 && ~isa(self.envObj.env{e, h}, 'Exp')
        expSojourn = false;
        break
    end
end
if expSojourn
    s_e = sum(self.E0(e, :));            % sojourn ~ Exp(s_e)
    d = size(Q, 1);
    piRes = s_e * (pi0 / (s_e * speye(d) - sparse(Q)));
    piExit_e = cell(1, E);
    for h = 1:E
        if self.E0(e, h) > 0
            piExit_e{h} = piRes;
        else
            piExit_e{h} = [];
        end
    end
    self.piExitDest{e} = piExit_e;
    self.piTimeAvg{e}  = piRes;
    results_e.statevec.ok = true;
    runtime = toc(T0);
    return
end

% General Markovian (PH/Erlang) sojourn: transient pi(t) = pi0 * exp(Q t) on the
% adaptive ode23 grid, averaged over the random holding-time / transition CDFs.
[pit, t] = ctmc_transient(Q, pi0, t0, t1);
t = t(:);

% Per-destination exit distributions: E[ pi(T_{e->h}) ] weighted by the
% increments of the e->h transition CDF proc{e}{h}.
piExit_e = cell(1, E);
for h = 1:E
    proc_eh = self.envObj.proc{e}{h};
    dF = map_cdf(proc_eh, t(2:end)) - map_cdf(proc_eh, t(1:end-1));
    w = [0; dF(:)];
    sw = sum(w);
    if sw > 0 && all(~isnan(w))
        piExit_e{h} = (w' * pit) / sw;
    else
        piExit_e{h} = [];
    end
end
self.piExitDest{e} = piExit_e;

% Sojourn-end distribution: E[ pi(T_sojourn) ] weighted by holdTime{e}.
holdT = self.envObj.holdTime{e};
dFh = map_cdf(holdT, t(2:end)) - map_cdf(holdT, t(1:end-1));
wh = [0; dFh(:)];
swh = sum(wh);
if swh > 0 && all(~isnan(wh))
    self.piTimeAvg{e} = (wh' * pit) / swh;
else
    self.piTimeAvg{e} = pit(end, :); % degenerate: use the terminal distribution
end

results_e.statevec.ok = true;
runtime = toc(T0);
end

% -------------------------------------------------------------------------
function post_(self, it)
% Chain the entry distributions: carry each stage's exit distributions into
% the stages they feed, weighted by the origin probabilities probOrig.
E = self.getNumberOfModels;
self.piEnterPrev = self.piEnter;
piEnterNew = cell(1, E);

for e = 1:E
    nstates_e = size(self.Qgen{e}, 1);
    acc = zeros(1, nstates_e);
    wsum = 0;
    for h = 1:E
        po = self.envObj.probOrig(h, e);
        if po > 0 && ~isempty(self.piExitDest{h}) && ~isempty(self.piExitDest{h}{e})
            pex = self.resetStateFun{h, e}(self.piExitDest{h}{e});
            pex = pex(:)';
            if numel(pex) ~= nstates_e
                line_error(mfilename, sprintf(['resetStateFun{%d,%d} returned a %d-element vector but ' ...
                    'stage %d has %d states. Supply a resetStateFun{%d,%d} that maps the state space of ' ...
                    'stage %d onto that of stage %d.'], h, e, numel(pex), e, nstates_e, h, e, h, e));
            end
            acc = acc + po * pex;
            wsum = wsum + po;
        end
    end
    if wsum > 0
        acc = acc / wsum;
    else
        acc = self.piEnter{e}; % no inflow this cycle: retain the current estimate
    end
    acc(acc < 0) = 0;
    s = sum(acc);
    if s > 0
        acc = acc / s;
    end
    piEnterNew{e} = acc;
end
self.piEnter = piEnterNew;
end

% -------------------------------------------------------------------------
function bool = converged_(self, it)
% Converged when the max L1 change across all entry distributions over a full
% cycle falls below iter_tol.
bool = false;
if it < 1 || isempty(self.piEnterPrev) || isempty(self.piEnter)
    return
end
E = self.getNumberOfModels;
l1 = 0;
for e = 1:E
    a = self.piEnter{e};
    b = self.piEnterPrev{e};
    if isempty(a) || isempty(b) || numel(a) ~= numel(b)
        return
    end
    l1 = max(l1, sum(abs(a(:) - b(:))));
end
if isnan(l1) || isinf(l1)
    return
end
if l1 < self.options.iter_tol
    bool = true;
    line_debug('ENV statevec converged: iteration %d, max L1 entry change %e < iter_tol %e', ...
        it, l1, self.options.iter_tol);
end
end

% -------------------------------------------------------------------------
function finish_(self)
% Environment-averaged blend: map each stage's sojourn-end distribution to
% marginal metrics and weight by the stage probability probEnv.
E = self.getNumberOfModels;
M = self.ensemble{1}.getNumberOfStations;
K = self.ensemble{1}.getNumberOfClasses;

Qval = zeros(M, K);
Uval = zeros(M, K);
Tval = zeros(M, K);

for e = 1:E
    piF = self.piTimeAvg{e};
    if isempty(piF)
        continue
    end
    data = self.statevecData{e};
    if strcmp(data.backend, 'mam')
        [QN, UN, ~, TN] = solver_mam_ldqbd_avg(data.ld, piF, data.levelOf);
    else
        [QN, UN, ~, TN] = solver_ctmc_avg_from_pi(data.sn, piF, self.SS{e}, ...
            self.SSaggr{e}, data.arvRates, data.depRates);
    end
    Qval = Qval + self.envObj.probEnv(e) * QN;
    Uval = Uval + self.envObj.probEnv(e) * UN;
    Tval = Tval + self.envObj.probEnv(e) * TN;
end

self.result.Avg.Q = Qval;
self.result.Avg.U = Uval;
self.result.Avg.T = Tval;

% Environment-blended cache hit/miss ratios written onto the stage-1 reference
% model. see _kb/06-solver-catalog.md for rationale
aggregateCacheBlend_(self, E, K);
end

% -------------------------------------------------------------------------
function aggregateCacheBlend_(self, E, K)
if strcmp(self.statevecData{1}.backend, 'mam')
    return % MAM backend has no enumerated cache state
end
sn1 = self.statevecData{1}.sn;
if ~isfield(sn1, 'nstateful')
    return
end
cacheStateful = [];
for isf = 1:sn1.nstateful
    ind = sn1.statefulToNode(isf);
    if sn1.nodetype(ind) == NodeType.Cache
        cacheStateful(end+1) = isf; %#ok<AGROW>
    end
end
if isempty(cacheStateful)
    return
end
for isf = cacheStateful
    hitT = zeros(1, K); missT = zeros(1, K);
    for e = 1:E
        piF = self.piTimeAvg{e};
        if isempty(piF); continue; end
        pv = piF(:); pv(pv < 0) = 0;
        if sum(pv) > 0; pv = pv / sum(pv); end
        dr = self.statevecData{e}.depRates;
        sne = self.statevecData{e}.sn;
        np = sne.nodeparam{sne.statefulToNode(isf)};
        w = self.envObj.probEnv(e);
        for k = 1:K
            if length(np.hitclass) >= k
                h = np.hitclass(k); mcl = np.missclass(k);
                if h > 0 && mcl > 0
                    hitT(k)  = hitT(k)  + w * (pv' * dr(:, isf, h));
                    missT(k) = missT(k) + w * (pv' * dr(:, isf, mcl));
                end
            end
        end
    end
    hitprob = NaN(1, K); missprob = NaN(1, K);
    for k = 1:K
        tot = hitT(k) + missT(k);
        if tot > 0
            hitprob(k)  = hitT(k) / tot;
            missprob(k) = missT(k) / tot;
        end
    end
    node = self.ensemble{1}.getNodeByIndex(sn1.statefulToNode(isf));
    node.setResultHitProb(hitprob);
    node.setResultMissProb(missprob);
end
end
