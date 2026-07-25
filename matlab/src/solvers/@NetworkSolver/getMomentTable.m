function [MomentTable, mom] = getMomentTable(self, order)
% GETMOMENTTABLE Exact higher moments of the per-class performance measures.
%
% [MOMENTTABLE, MOM] = GETMOMENTTABLE(SELF) returns a table with one row per
% (Station, JobClass) giving, in addition to the means that getAvgTable
% reports, the second moments of that row's queue length and response time:
%   QLen, QLenVar, QLenSCV, RespT, RespTVar, RespTSCV
%
% [..] = GETMOMENTTABLE(SELF, ORDER) selects which moment orders to report.
% ORDER is a set. A scalar k is read as 1:k, "everything up to order k"; an
% explicit vector selects exactly those orders:
%   1        the means only:            QLen, RespT
%   2        (default) means and second moments, i.e. the columns above
%   3        also adds RespTSkew
%   [1 2]    the same as 2
%   [2 3]    second moments and skewness, without the means
% Order 1 contributes QLen and RespT, order 2 contributes the Var and SCV
% columns, order 3 contributes RespTSkew.
%
% ORDER = 3 also adds QLenSkew, the skewness of the per-class queue length.
% That quantity is reachable because the generating parameter need not scale a
% whole demand column: scaling L(i,r) alone is Theorem 1 of Akyildiz and Strelen
% with the class subset T = {r}, and it generates the moments of n(i,r) itself.
% QLenSkew is available only for closed single-server models, which is the scope
% of pfqn_sens_mom; it is NaN otherwise.
%
% All of it is exact, not simulated and not approximated. The queue-length
% moments come from the product-form identity Cov[n(i,r),n(j,s)] =
% L(j,s) dQ(i,r)/dL(j,s), evaluated by the pfqn_sens_* family; see
% _kb/03-api-layer.md.
%
% RESPONSE-TIME MOMENTS ARE FCFS OR PROCESSOR-SHARING. RespTVar and RespTSCV
% are NaN at any station that is neither, and at an LCFS center in particular,
% because the sojourn-time distribution there is not known in general (Strelen
% 1990, Section 4) and a wrong value is worse than a blank. The mean RespT is
% always reported, since it needs no distributional result.
%
% The FCFS moments come from pfqn_sens_respt and are closed-model only. The
% processor-sharing moments come from Mitra and Morrison (1983) and cover two
% configurations, both requiring exponential single-server service:
%   purely open   -> qsys_mm1_ps, exact, at any PS station whose arrivals are
%                    Poisson, that is, that lies on no routing cycle
%   purely closed -> pfqn_respt_ps_moments, for the terminal-driven system the
%                    paper analyses: one PS station visited once per think
%                    cycle, with delay stations holding the think time
% A PS station outside those configurations keeps RespTVar = NaN.
%
% Scope by model type:
%   closed, single-server            -> pfqn_sens_mva
%   closed, multiserver              -> pfqn_sens_mvaldmx
%   mixed open and closed            -> pfqn_sens_mvaldmx
%   purely open, single-server       -> exact BCMP closed form (below)
%   purely open, multiserver         -> not supported, see the error
%
% MOM is a struct carrying the raw results: .qlen is the underlying
% pfqn_sens_mva / pfqn_sens_mvaldmx struct (with the full covariance matrices,
% not just the diagonal this table shows), .respt is the pfqn_sens_respt struct
% or empty, and .psrespt is the pfqn_respt_ps_moments struct or empty. Use it
% when the per-pair covariances, or the route taken at a PS station, are
% needed.
%
% Per-station TOTAL moments, including the third moment and the skewness, are
% in getMomentStationTable: they are only defined for a station total, because
% the parameter that generates them scales a whole demand column.
%
% See also: getAvgTable, getSensitivityTable, getMomentStationTable.

if nargin < 2 || isempty(order)
    order = 2;
end
order = validateMomentOrder(order, 3);

sn = self.model.getStruct();
R = sn.nclasses;
N = sn.njobs;

[lambda, D, Np, Z, mu, Ssrv, ~] = sn_get_product_form_params(sn);
queueIndices = find(sn.nodetype == NodeType.Queue);
Mq = numel(queueIndices);
Ztot = sum(Z, 1);
isOpen = any(isinf(N));
isClosed = any(isfinite(N) & N > 0);
isMixed = isOpen && isClosed;

mom = struct('qlen', [], 'respt', [], 'qlenmom', [], 'psrespt', []);
QLen = zeros(Mq, R);
QLenVar = zeros(Mq, R);

% ---- queue-length moments -------------------------------------------------
if ~isOpen
    if all(Ssrv == 1)
        mom.qlen = pfqn_sens_mva(D, Np, Ztot);
    else
        mom.qlen = pfqn_sens_mvaldmx(zeros(1,R), D, Np, Ztot, mu, Ssrv);
    end
    QLen = mom.qlen.Q;
    QLenVar = mom.qlen.QVar;
elseif isMixed
    mom.qlen = pfqn_sens_mvaldmx(lambda, D, N, Ztot, mu, Ssrv);
    QLen = mom.qlen.Q;
    QLenVar = mom.qlen.QVar;
else
    % Purely open BCMP single-server: closed-form geometric/multinomial joint
    % law, no lattice recursion. see _kb/06-solver-catalog.md for rationale
    if any(Ssrv > 1)
        line_error(mfilename, 'getMomentTable does not support multiserver stations in a purely open model: the queue-length law is not geometric there. Add a closed class, or use a single-server model.');
    end
    rho = zeros(Mq, R);
    for ist = 1:Mq
        for r = 1:R
            if isinf(N(r))
                rho(ist, r) = lambda(r) * D(ist, r);
            end
        end
    end
    rhoTot = sum(rho, 2);
    for ist = 1:Mq
        ri = rhoTot(ist);
        if ri >= 1
            line_error(mfilename, sprintf('Station %s is unstable (utilization %.4f >= 1); its queue-length moments do not exist.', sn.nodenames{queueIndices(ist)}, ri));
        end
        if ri <= 0
            continue;
        end
        En = ri / (1 - ri);
        Vn = ri / (1 - ri)^2;
        for r = 1:R
            pr = rho(ist, r) / ri;
            QLen(ist, r) = pr * En;
            QLenVar(ist, r) = En * (pr - pr^2) + pr^2 * Vn;
        end
    end
end

% ---- per-class queue-length skewness (closed, single-server only) ---------
% pfqn_sens_mom with groups = 1:R scales one class at a time, which is the
% class-subset parameter T = {r} of Akyildiz and Strelen's Theorem 1, so it
% yields the moments of n(i,r) rather than of the station total.
QLenSkew = nan(Mq, R);
mom.qlenmom = [];
if any(order == 3) && ~isOpen && all(Ssrv == 1)
    mom.qlenmom = pfqn_sens_mom(D, Np, Ztot, ones(1,Mq), 1:R);
    QLenSkew = mom.qlenmom.Skew;
end

% ---- response-time moments (FCFS only) ------------------------------------
RespT = zeros(Mq, R);
RespTVar = nan(Mq, R);
RespTSkew = nan(Mq, R);
% see _kb/06-solver-catalog.md for rationale (pfqn_sens_respt reads mu only at FCFS)
respAvail = ~isOpen && fcfsRatesAreClassIndependent(sn, queueIndices, R);
if respAvail
    Ssvc = ones(Mq, 1);
    Vq = zeros(Mq, R);
    for ist = 1:Mq
        sIdx = sn.nodeToStation(queueIndices(ist));
        if sn.sched(sIdx) == SchedStrategy.FCFS
            st = serviceTimeOf(sn, sIdx, R);
            if st > 0
                Ssvc(ist) = st;   % the true per-visit rate; (4.5) needs it
            end
        end
        for r = 1:R
            Vq(ist, r) = D(ist, r) / Ssvc(ist);
        end
    end
    mom.respt = pfqn_sens_respt(Ssvc, Vq, Np, Ztot, Ssrv(:), max(order));
end
for ist = 1:Mq
    sIdx = sn.nodeToStation(queueIndices(ist));
    for r = 1:R
        if D(ist, r) <= 0
            continue;
        end
        % Mean response time per visit by Little's law at the station; always
        % reported. see _kb/06-solver-catalog.md for rationale
        xr = throughputOf(mom, lambda, N, r);
        Vir = D(ist, r) * sn.rates(sIdx, r);
        if xr > 0 && Vir > 0
            RespT(ist, r) = QLen(ist, r) / (xr * Vir);
        end
    end
    % The variance needs the sojourn-time distribution, which is known only at
    % FCFS centers; elsewhere RespTVar stays NaN.
    if ~isempty(mom.respt) && sn.sched(sIdx) == SchedStrategy.FCFS
        for r = 1:R
            if D(ist, r) > 0
                RespT(ist, r) = mom.respt.W(ist, r);
                if max(order) >= 2
                    RespTVar(ist, r) = mom.respt.WVar(ist, r);
                end
                if max(order) >= 3
                    RespTSkew(ist, r) = mom.respt.WSkew(ist, r);
                end
            end
        end
    end
end

% ---- response-time moments at processor-sharing stations ------------------
% Mitra and Morrison (1983) supply the sojourn-time moments the FCFS block
% above cannot reach: exactly at an open PS station fed by Poisson streams,
% and by expansion (or by exact enumeration when small) in the closed
% terminal-driven system. see _kb/05-solvers-overview.md for the scope
mom.psrespt = [];
if any(order == 2) && ~isMixed && Mq > 0
    psQueues = [];
    for ist = 1:Mq
        if sn.sched(sn.nodeToStation(queueIndices(ist))) == SchedStrategy.PS
            psQueues(end+1) = ist; %#ok<AGROW>
        end
    end
    if isOpen
        for ist = psQueues
            sIdx = sn.nodeToStation(queueIndices(ist));
            visiting = find(D(ist, :) > 0);
            if Ssrv(ist) ~= 1 || ~ratesAreExponential(sn, sIdx, visiting) ...
                    || ~stationIsFeedbackFree(sn, sIdx, R)
                continue
            end
            muSt = ones(1, R);
            lamSt = zeros(1, R);
            for r = visiting
                muSt(r) = sn.rates(sIdx, r);
                lamSt(r) = lambda(r) * D(ist, r) * sn.rates(sIdx, r);
            end
            if sum(lamSt ./ muSt) >= 1
                continue
            end
            [Wps, W2ps] = qsys_mm1_ps(lamSt, muSt);
            for r = visiting
                RespTVar(ist, r) = W2ps(r) - Wps(r)^2;
            end
        end
    elseif numel(psQueues) == 1 && Mq == 1
        % the paper's closed system: terminals in series with one PS CPU
        ist = psQueues(1);
        sIdx = sn.nodeToStation(queueIndices(ist));
        visiting = find(D(ist, :) > 0);
        Vps = D(ist, visiting) .* sn.rates(sIdx, visiting);
        singleVisit = all(abs(Vps - 1) <= GlobalConstants.CoarseTol);
        if Ssrv(ist) == 1 && ratesAreExponential(sn, sIdx, visiting) && singleVisit ...
                && all(Ztot(visiting) > 0)
            Sps = ones(1, R);
            Zps = ones(1, R);
            Nps = zeros(1, R);
            for r = visiting
                Sps(r) = 1 / sn.rates(sIdx, r);
                Zps(r) = Ztot(r);
                Nps(r) = N(r);
            end
            [Wps, W2ps, mom.psrespt] = pfqn_respt_ps_moments(Sps, Nps, Zps);
            for r = visiting
                RespTVar(ist, r) = W2ps(r) - Wps(r)^2;
            end
        end
    end
end

% ---- assemble -------------------------------------------------------------
Station = {}; JobClass = {};
QLenv = []; QLenVarv = []; QLenSCVv = [];
RespTv = []; RespTVarv = []; RespTSCVv = []; RespTSkewv = []; QLenSkewv = [];
for ist = 1:Mq
    node = queueIndices(ist);
    for r = 1:R
        if D(ist, r) <= 0
            continue;   % class r does not visit this station
        end
        Station{end+1, 1} = sn.nodenames{node}; %#ok<AGROW>
        JobClass{end+1, 1} = sn.classnames{r};  %#ok<AGROW>
        QLenv(end+1, 1) = QLen(ist, r);         %#ok<AGROW>
        QLenVarv(end+1, 1) = QLenVar(ist, r);   %#ok<AGROW>
        QLenSCVv(end+1, 1) = scv(QLenVar(ist, r), QLen(ist, r)); %#ok<AGROW>
        RespTv(end+1, 1) = RespT(ist, r);       %#ok<AGROW>
        RespTVarv(end+1, 1) = RespTVar(ist, r); %#ok<AGROW>
        RespTSCVv(end+1, 1) = scv(RespTVar(ist, r), RespT(ist, r)); %#ok<AGROW>
        RespTSkewv(end+1, 1) = RespTSkew(ist, r); %#ok<AGROW>
        QLenSkewv(end+1, 1) = QLenSkew(ist, r);   %#ok<AGROW>
    end
end

vars = {Station, JobClass};
names = {'Station', 'JobClass'};
if any(order == 1)
    vars{end+1} = QLenv;      names{end+1} = 'QLen';
end
if any(order == 2)
    vars{end+1} = QLenVarv;   names{end+1} = 'QLenVar';
    vars{end+1} = QLenSCVv;   names{end+1} = 'QLenSCV';
end
if any(order == 3)
    vars{end+1} = QLenSkewv;  names{end+1} = 'QLenSkew';
end
if any(order == 1)
    vars{end+1} = RespTv;     names{end+1} = 'RespT';
end
if any(order == 2)
    vars{end+1} = RespTVarv;  names{end+1} = 'RespTVar';
    vars{end+1} = RespTSCVv;  names{end+1} = 'RespTSCV';
end
if any(order == 3)
    vars{end+1} = RespTSkewv; names{end+1} = 'RespTSkew';
end
MomentTable = table(vars{:}, 'VariableNames', names);
end

% =========================================================================
function order = validateMomentOrder(order, maxorder)
% ORDER is a set of moment orders. A scalar k is shorthand for 1:k, so that
% getMomentTable(2) means "up to the second moment" and not "the second moment
% alone"; a vector of two or more entries is taken literally.
%
% Consequence of MATLAB's isscalar: a one-element vector IS a scalar, so [2]
% takes the 1:k path and yields [1 2]. "The second moment alone" is therefore
% not expressible, which is deliberate rather than overlooked: a variance with
% no mean beside it is not a useful table, and [2 3] remains available for the
% higher orders without the means.
%
% Non-integers are rejected on BOTH paths. Rounding them silently would accept
% [1 2.5] as [1 3], i.e. answer a question that was not asked.
if ~isnumeric(order) || isempty(order) || any(~isfinite(order(:)))
    line_error(mfilename, sprintf('order must be an integer in 1..%d, or a vector of such integers.', maxorder));
end
if any(order(:) ~= round(order(:))) || any(order(:) < 1) || any(order(:) > maxorder)
    line_error(mfilename, sprintf('order must be an integer in 1..%d, or a vector of such integers.', maxorder));
end
if isscalar(order)
    order = 1:order;
    return;
end
order = unique(order(:)');
end

% =========================================================================
function v = scv(variance, meanv)
% Squared coefficient of variation. NaN when the mean is zero, since the SCV is
% then undefined rather than infinite in any useful sense.
if isnan(variance) || meanv <= 0
    v = NaN;
else
    v = variance / meanv^2;
end
end

% =========================================================================
function ok = fcfsRatesAreClassIndependent(sn, queueIndices, R)
% An FCFS station in a BCMP network must serve every class at the same
% exponential rate; that is the precondition for reading a per-visit rate off
% it. Non-FCFS stations are unconstrained here, see the caller.
ok = true;
for ist = 1:numel(queueIndices)
    sIdx = sn.nodeToStation(queueIndices(ist));
    if sn.sched(sIdx) ~= SchedStrategy.FCFS
        continue;
    end
    ref = -1;
    for r = 1:R
        rate = sn.rates(sIdx, r);
        if ~isfinite(rate) || rate <= 0
            continue;
        end
        if ref < 0
            ref = rate;
        elseif abs(rate - ref) > GlobalConstants.FineTol * max(1, ref)
            ok = false;
            return;
        end
    end
end
end

% =========================================================================
function s = serviceTimeOf(sn, sIdx, R)
% The common service time of a station, i.e. the reciprocal of the class-
% independent rate. Zero if no class is served here.
s = 0;
for r = 1:R
    rate = sn.rates(sIdx, r);
    if isfinite(rate) && rate > 0
        s = 1 / rate;
        return;
    end
end
end

% =========================================================================
function ok = ratesAreExponential(sn, sIdx, classes)
% The PS sojourn-time moments of Mitra and Morrison assume exponential
% service; a phase-type service leaves the mean intact but not the moments.
ok = true;
for r = classes
    if sn.procid(sIdx, r) ~= ProcessType.EXP
        ok = false;
        return;
    end
end
end

% =========================================================================
function ok = stationIsFeedbackFree(sn, sIdx, R)
% An open PS station sees Poisson arrivals only when no job can return to it,
% which is Melamed's condition: the station must lie on no routing cycle. The
% test aggregates the classes, so a cycle closed through a class switch also
% disqualifies the station.
M = sn.nstations;
adj = false(M, M);
for i = 1:M
    nt = sn.nodetype(sn.stationToNode(i));
    if nt == NodeType.Source || nt == NodeType.Sink
        continue   % jobs never leave a sink, and the rt arc back to the source is bookkeeping
    end
    for j = 1:M
        blk = sn.rt((i-1)*R + (1:R), (j-1)*R + (1:R));
        adj(i, j) = any(blk(:) > 0);
    end
end
reach = adj(sIdx, :);
frontier = find(reach);
while ~isempty(frontier)
    nxt = any(adj(frontier, :), 1) & ~reach;
    reach = reach | nxt;
    frontier = find(nxt);
end
ok = ~reach(sIdx);
end

% =========================================================================
function x = throughputOf(mom, lambda, N, r)
if isinf(N(r))
    x = lambda(r);
elseif ~isempty(mom.qlen)
    x = mom.qlen.X(r);
else
    x = 0;
end
end
