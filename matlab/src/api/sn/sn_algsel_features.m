function f = sn_algsel_features(sn)
% F = SN_ALGSEL_FEATURES(SN)
%
% Structural feature map f(x) used by the learned solver selector
% (chooseSolverTree). Mirrors line-test.git/algsel/features.py field for
% field: the deployed tree is fitted on those features, so any divergence
% here silently moves the decision boundaries.
%
% Every feature must be cheap relative to solving the model (structural reads
% of sn plus O(M*K) arithmetic), otherwise the selector costs more than the
% choice it makes. Feature names are a stable schema: append, never rename.

f = struct();

M = double(sn.nstations);
K = double(sn.nclasses);
f.f_nstations = M;
f.f_nnodes = double(sn.nnodes);
f.f_nclasses = K;
f.f_nchains = double(sn.nchains);
f.f_classes_per_chain = K / max(double(sn.nchains), 1);

njobs = double(sn.njobs(:))';
isOpen = isinf(njobs);
isClosed = isfinite(njobs) & njobs > 0;
f.f_n_open = sum(isOpen);
f.f_n_closed = sum(isClosed);
f.f_frac_open = sum(isOpen) / max(K, 1);
f.f_is_open = double(all(isOpen));
f.f_is_closed = double(all(isClosed));
f.f_is_mixed = double(any(isOpen) && any(isClosed));

totalJobs = 0;
if any(isClosed)
    totalJobs = sum(njobs(isClosed));
end
f.f_total_jobs = totalJobs;
f.f_log_total_jobs = log1p(totalJobs);
f.f_jobs_per_chain = totalJobs / max(double(sn.nchains), 1);
f.f_jobs_per_station = totalJobs / max(M, 1);
if any(isClosed)
    f.f_max_class_jobs = max(njobs(isClosed));
else
    f.f_max_class_jobs = 0;
end

% SchedStrategy.toText is lowercase; features.py keys off the enum name, so
% the comparison is done in upper case in both codebases.
schedNames = cell(1, M);
for i = 1:M
    schedNames{i} = upper(char(SchedStrategy.toText(sn.sched(i))));
end
pool = {'PS','FCFS','INF','LCFSPR','SIRO','HOL','DPS','GPS','LCFS'};
denom = max(numel(schedNames), 1);
for s = 1:numel(pool)
    f.(sprintf('f_frac_sched_%s', lower(pool{s}))) = ...
        sum(strcmp(schedNames, pool{s})) / denom;
end
f.f_n_distinct_sched = numel(unique(schedNames));
hasPrio = false;
for i = 1:M
    if ~isempty(strfind(schedNames{i}, 'PRIO')) || strcmp(schedNames{i}, 'HOL') %#ok<STREMP>
        hasPrio = true;
    end
end
f.f_has_prio_sched = double(hasPrio);

nservers = double(sn.nservers(:))';
finiteSrv = nservers(isfinite(nservers));
if isempty(finiteSrv)
    finiteSrv = 1;
end
f.f_max_servers = max(finiteSrv);
f.f_mean_servers = mean(finiteSrv);
f.f_has_multiserver = double(any(finiteSrv > 1));
f.f_n_delay_stations = sum(strcmp(schedNames, 'INF'));

scv = sn.scv(:)';
scv = scv(isfinite(scv) & scv > 0);
if isempty(scv)
    scv = 1;
end
f.f_scv_min = min(scv);
f.f_scv_max = max(scv);
f.f_scv_mean = mean(scv);
f.f_scv_logmean = mean(log(scv));
f.f_frac_exp = mean(abs(scv - 1) < 1e-9);
f.f_frac_hypo = mean(scv < 1 - 1e-9);
f.f_frac_hyper = mean(scv > 1 + 1e-9);

phases = [];
if isfield(sn, 'phases') && ~isempty(sn.phases)
    phases = double(sn.phases(:))';
    phases = phases(isfinite(phases) & phases > 0);
end
if isempty(phases)
    phases = 1;
end
f.f_max_phases = max(phases);
f.f_mean_phases = mean(phases);
f.f_log_total_phases = sum(log(phases));

% MATLAB sn.visits{c} is STATEFUL-indexed, python's is station-indexed; the
% station rows are pulled through stationToStateful so both sum the same
% (nstations x nclasses) matrix.
rates = sn.rates;
V = zeros(M, K);
if isfield(sn, 'visits') && ~isempty(sn.visits)
    for c = 1:numel(sn.visits)
        Vc = sn.visits{c};
        if isempty(Vc)
            continue
        end
        for ist = 1:M
            isf = sn.stationToStateful(ist);
            if isf >= 1 && isf <= size(Vc, 1)
                row = Vc(isf, :);
                row(~isfinite(row)) = 0;
                V(ist, :) = V(ist, :) + row(1:K);
            end
        end
    end
end
D = zeros(M, K);
pos = rates > 0 & isfinite(rates);
D(pos) = V(pos) ./ rates(pos);
D(~isfinite(D)) = 0;
dst = sum(D, 2)';
dst = dst(isfinite(dst));
if isempty(dst)
    dst = 0;
end
dmax = max(dst);
f.f_demand_max = dmax;
f.f_demand_mean = mean(dst);
if dmax > 0
    f.f_demand_cv = std(dst, 1) / dmax;
else
    f.f_demand_cv = 0;
end
ordered = sort(dst, 'descend');
if numel(ordered) > 1 && ordered(2) > 0
    f.f_bottleneck_ratio = ordered(1) / ordered(2);
else
    f.f_bottleneck_ratio = 1;
end
f.f_log_demand_max = log1p(dmax);

lambda = 0;
srcIdx = find(strcmp(schedNames, 'EXT'), 1);
if ~isempty(srcIdx)
    row = rates(srcIdx, :);
    lambda = sum(row(isfinite(row)));
end
f.f_arrival_rate = lambda;
if lambda > 0
    f.f_rho_max_open = min(lambda * dmax, 10);
else
    f.f_rho_max_open = 0;
end

nodeTypes = double(sn.nodetype(:))';
f.f_n_fork = sum(nodeTypes == NodeType.Fork);
f.f_n_join = sum(nodeTypes == NodeType.Join);
f.f_has_fj = double(f.f_n_fork > 0);
f.f_has_classswitch = double(any(nodeTypes == NodeType.ClassSwitch));
f.f_has_cache = double(any(nodeTypes == NodeType.Cache));
f.f_has_source = double(any(nodeTypes == NodeType.Source));

if isfield(sn, 'rt') && ~isempty(sn.rt)
    rt = full(sn.rt);
    f.f_routing_density = mean(rt(:) > 0);
else
    f.f_routing_density = 0;
end

f.f_has_ld = 0;
if isfield(sn, 'lldscaling') && ~isempty(sn.lldscaling)
    f.f_has_ld = double(any(abs(sn.lldscaling(:) - 1) > 1e-12));
end
f.f_has_finite_cap = 0;
if isfield(sn, 'cap') && ~isempty(sn.cap)
    f.f_has_finite_cap = double(any(isfinite(sn.cap(:))));
end
f.f_has_class_prio = 0;
if isfield(sn, 'classprio') && ~isempty(sn.classprio)
    f.f_has_class_prio = double(numel(unique(sn.classprio(:))) > 1);
end

f.f_product_form = double(sn_has_product_form(sn));

% log state-space bound: multiset of the closed population over stations,
% inflated by the service phases, which is what drives CTMC feasibility.
logSS = 0;
if any(isClosed)
    chains = [];
    if isfield(sn, 'chains') && ~isempty(sn.chains)
        chains = sn.chains;
    end
    for c = 1:double(sn.nchains)
        inC = [];
        if ~isempty(chains)
            inC = find(chains(c, :) > 0);
        end
        pop = 0;
        for r = inC
            if r <= numel(njobs) && isfinite(njobs(r))
                pop = pop + njobs(r);
            end
        end
        if pop > 0
            logSS = logSS + gammaln(pop + M) - gammaln(M) - gammaln(pop + 1);
        end
    end
end
logSS = logSS + sum(log(phases)) * max(totalJobs, 1) / max(M, 1);
f.f_log_state_space = min(logSS, 500);
f.f_ctmc_feasible = double(logSS < log(2e4));
end
