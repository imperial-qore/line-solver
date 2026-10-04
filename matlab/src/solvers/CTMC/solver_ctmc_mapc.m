function sn = solver_ctmc_mapc(sn)
% SN = SOLVER_CTMC_MAPC(SN)
%
% Rewrite the MAP (and MMPP2) service of every multiserver FCFS station into
% the pair form the CTMC needs to follow the law JMT and both LDES engines
% sample: ONE station-level sampler per (station, class), whose every draw
% starts in the phase the previous draw ENDED in, draws chained in service-start
% order and idle periods included.
%
% With c > 1 servers the next start can happen while earlier draws are still in
% progress, so the landing phase of a draw must be known when it starts. With
% V = (-D0)\D1, a draw started in h ends in j with probability V(h,j). A busy
% server is therefore a PAIR (i,j): current phase i and predetermined landing j,
% which moves i->k at D0(i,k)V(k,j)/V(i,j) and completes at D1(i,j)/V(i,j) (the
% Doob transform of D0 conditioned on ending in j). The per-class memory local
% variable keeps its meaning of carried phase h in 1..p, but it is now the
% landing phase of the most recently STARTED draw: a start from h enters pair
% (h,j) with probability V(h,j) and sets h := j, while phase moves and
% completions leave it alone. For a renewal MAP (V(i,j) = alpha(j)) the chain
% reduces in law to PH/c; with c = 1 the station is left untouched.
%
% The pair law is stored as the lifted MAP D0p (conditioned moves) and
% D1p((i,j),(j,j')) = D1(i,j)/V(i,j)*V(j,j'), which is equivalent in law to the
% original MAP. The bookkeeping that State.afterEventStation and
% State.fromMarginal read is sn.ctmcmapc{ist,r}: p, pairs (T x 2), V and done.
% The rewrite is idempotent and rebuilds sn.state of the rewritten stations.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;
R = sn.nclasses;
if ~isfield(sn,'ctmcmapc') || isempty(sn.ctmcmapc)
    sn.ctmcmapc = cell(M, R);
end
snOld = sn;
changed = false(1, M);
for ist = 1:M
    ind = sn.stationToNode(ist);
    if sn.nodetype(ind) == NodeType.Source || sn.sched(ist) ~= SchedStrategy.FCFS ...
            || ~isfinite(sn.nservers(ist)) || sn.nservers(ist) <= 1
        continue
    end
    for r = 1:R
        if ~isempty(sn.ctmcmapc{ist,r}) || ~(sn.procid(ist,r) == ProcessType.MAP || sn.procid(ist,r) == ProcessType.MMPP2)
            continue
        end
        D0 = sn.proc{ist}{r}{1};
        D1 = sn.proc{ist}{r}{2};
        p = size(D0,1);
        V = (-D0) \ D1;
        V(abs(V) < 1e-14) = 0;
        [jj, ii] = find(V' > 0); % pairs (i,j) with V(i,j) > 0, i-major
        pairs = [ii, jj];
        T = size(pairs,1);
        tp = zeros(p,p);
        for t = 1:T
            tp(pairs(t,1), pairs(t,2)) = t;
        end
        done = zeros(T,1);
        D0p = zeros(T);
        D1p = zeros(T);
        for t = 1:T
            i = pairs(t,1);
            j = pairs(t,2);
            done(t) = D1(i,j) / V(i,j);
            D0p(t,t) = D0(i,i);
            for k = setdiff(find(D0(i,:) ~= 0), i)
                if tp(k,j) > 0
                    D0p(t,tp(k,j)) = D0(i,k) * V(k,j) / V(i,j);
                end
            end
            for jn = 1:p
                if tp(j,jn) > 0
                    D1p(t,tp(j,jn)) = done(t) * V(j,jn);
                end
            end
        end
        sn.proc{ist}{r} = {D0p, D1p};
        sn.pie{ist}{r} = map_pie({D0p, D1p});
        sn.mu{ist}{r} = -diag(D0p);
        sn.phi{ist}{r} = done ./ (-diag(D0p));
        sn.phases(ist,r) = T;
        sn.ctmcmapc{ist,r} = struct('p', p, 'pairs', pairs, 'V', V, 'done', done);
        changed(ist) = true;
    end
end
if ~any(changed)
    return
end
sn.phasessz = max(sn.phases, ones(size(sn.phases)));
sn.phasessz(sn.nodeToStation(sn.nodetype == NodeType.Join),:) = sn.phases(sn.nodeToStation(sn.nodetype == NodeType.Join),:);
if isfield(sn,'markidx') && ~isempty(sn.markidx)
    sn.phasessz(sn.markidx > 1) = 1;
end
sn.phaseshift = [zeros(size(sn.phases,1),1), cumsum(sn.phasessz,2)];
if ~isfield(sn,'state') || isempty(sn.state)
    return
end
for ist = find(changed)
    isf = sn.stationToStateful(ist);
    if numel(sn.state) < isf || isempty(sn.state{isf})
        continue
    end
    sn.state{isf} = sub_rebuild_state(snOld, sn, ist, sn.state{isf});
end
end

function st = sub_rebuild_state(snOld, sn, ist, old)
% A job in service in MAP phase i is placed in the first pair (i,j); the buffer
% and the local variables, the carried phase included, are unchanged.
R = sn.nclasses;
ind = sn.stationToNode(ist);
V = sum(snOld.nvars(ind,:));
Kold = snOld.phasessz(ist,:);
st = zeros(size(old,1), size(old,2) - sum(Kold) + sum(sn.phasessz(ist,:)));
for row = 1:size(old,1)
    oldrow = old(row,:);
    buf = oldrow(1:(end-sum(Kold)-V));
    srvOld = oldrow((end-sum(Kold)-V+1):(end-V));
    vars = oldrow((end-V+1):end);
    srv = zeros(1, sum(sn.phasessz(ist,:)));
    for r = 1:R
        blk = srvOld((snOld.phaseshift(ist,r)+1):(snOld.phaseshift(ist,r)+Kold(r)));
        if isempty(sn.ctmcmapc{ist,r})
            srv((sn.phaseshift(ist,r)+1):(sn.phaseshift(ist,r)+sn.phasessz(ist,r))) = blk;
        else
            pairs = sn.ctmcmapc{ist,r}.pairs;
            for i = find(blk > 0)
                t = find(pairs(:,1) == i, 1);
                srv(sn.phaseshift(ist,r) + t) = srv(sn.phaseshift(ist,r) + t) + blk(i);
            end
        end
    end
    st(row,:) = [buf, srv, vars];
end
end
