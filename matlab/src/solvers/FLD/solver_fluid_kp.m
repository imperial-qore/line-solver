function [QN,UN,RN,TN,xvec_it,QNt,UNt,TNt,xvec_t,t,iters,runtime,QVart,Sigmat] = solver_fluid_kp(sn, options)
% [QN,UN,RN,TN,XVEC_IT,QNT,UNT,TNT,XVEC_T,T,ITERS,RUNTIME,QVART,SIGMAT] = SOLVER_FLUID_KP(SN, OPTIONS)
%
% Fluid and diffusion limits of the (MAP_t/Ph_t/inf)^N network of Y. M. Ko and
% J. Pender, "Diffusion limits for the (MAP_t/Ph_t/inf)^N queueing network",
% Oper. Res. Lett. 45 (2017) 248-253.
%
% The mean and the covariance of the limit are integrated jointly:
%
%   dq/dt     = F(t,q) = A f(t,q)
%   dSigma/dt = J Sigma + Sigma J' + G,   J = A df/dq,  G = A diag(f) A'
%
% with A the jump matrix whose column e is the jump vector of event e and f the
% event rate vector. G is exactly dH dH' of Theorem 3.3, each independent Poisson
% term contributing l_e l_e' f_e. Where f is affine in q -- infinite-server
% stations and the arrival phase process -- J does not depend on q and both
% equations close exactly, so for the (MAP_t/Ph_t/inf)^N case the mean and the
% covariance are exact rather than asymptotic. Finite-server stations are
% admitted through the usual fluid min(x,c) term, where the covariance degrades
% to a linear-noise approximation and a warning is emitted.
%
% This method does NOT reuse the closing ODE. That formulation routes a departure
% from the source to the destination station and returns mass through the
% STATIONARY arrival-instant vector pie, replacing the D1' operator by the
% rank-one map pie*(D1*e)', i.e. by the PH renewal process with representation
% (pie, D0). Its stationary arrival rate is exact but its autocorrelation is
% gone, and a non-renewal arrival stream is the entire point of a MAP.
%
% State layout, station-major, arrival phases before service phases:
%   u-block  one per (EXT station, class): arrival MAP phase occupancy, sum 1
%   x-block  one per (queueing station, class): fluid count in each service phase
%
% MAPt, PHt and NHPP all reach the integrator as a piecewise-constant (D0, D1)
% schedule; an NHPP is the one-phase MAPt case, lifted in LOCAL_SCHEDULE.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

runtime_start = tic;
M = sn.nstations;
K = sn.nclasses;
rt = sn.rt;

if any(isfinite(sn.njobs))
    line_error(mfilename, ['the ''kp'' method analyses the open (MAP_t/Ph_t/inf)^N network ' ...
        'of Ko and Pender (2017); a closed class has no arrival process to modulate. ' ...
        'Use ''closing'' or ''matrix'' for closed models.']);
end

% ---- blocks -------------------------------------------------------------
ublocks = zeros(0,4); % [station class offset nphases]
xblocks = zeros(0,4);
off = 0;
for ist = 1:M
    isExt = sn.sched(ist) == SchedStrategy.EXT;
    for r = 1:K
        h = sn.phases(ist,r);
        if h <= 0 || ~isfinite(sn.rates(ist,r)) || sn.rates(ist,r) <= 0
            continue % disabled class at this station
        end
        if isExt
            ublocks(end+1,:) = [ist r off h]; %#ok<AGROW>
        else
            xblocks(end+1,:) = [ist r off h]; %#ok<AGROW>
        end
        off = off + h;
    end
end
dim = off;
if isempty(ublocks)
    line_error(mfilename, 'the ''kp'' method needs at least one Source with an arrival process');
end

linearModel = true;
for b = 1:size(xblocks,1)
    ist = xblocks(b,1);
    if sn.sched(ist) ~= SchedStrategy.INF && isfinite(sn.nservers(ist))
        linearModel = false;
    end
end
if ~linearModel
    line_warning(mfilename, ['a finite-server station makes the rate functions nonlinear, so the ' ...
        'covariance is a linear-noise approximation rather than the exact second moment; ' ...
        'it is exact for infinite-server stations.']);
end

% ---- schedules ----------------------------------------------------------
sched = struct('station',{},'class',{},'breakpoints',{},'nseg',{},'cyclic',{}, ...
    'segD0',{},'segD1',{});
for ist = 1:M
    for r = 1:K
        [isSched,~,~,bp,sD0,sD1,cyc] = local_schedule(sn, ist, r);
        if isSched
            % assign the cell fields directly: struct() REPLICATES over a cell
            % argument, so neither struct(...,'segD0',sD0) nor {sD0} nor {{sD0}}
            % stores the n-segment cell itself
            entry = struct('station',ist,'class',r,'breakpoints',bp, ...
                'nseg',numel(sD0),'cyclic',cyc,'segD0',[],'segD1',[]);
            entry.segD0 = sD0;
            entry.segD1 = sD1;
            if isempty(sched)
                sched = entry;
            else
                sched(end+1) = entry; %#ok<AGROW>
            end
        end
    end
end

% nominal (D0,D1) per station-class, used where there is no schedule
nomD0 = cell(M,K); nomD1 = cell(M,K); nomPie = cell(M,K);
for ist = 1:M
    for r = 1:K
        if sn.phases(ist,r) <= 0, continue; end
        [isSched, D0bar, D1bar] = local_schedule(sn, ist, r);
        if isSched
            nomD0{ist,r} = D0bar;
            nomD1{ist,r} = D1bar;
        else
            pr = sn.proc{ist}{r};
            if iscell(pr) && numel(pr) >= 2 && all(size(pr{1}) == size(pr{2}))
                nomD0{ist,r} = pr{1};
                nomD1{ist,r} = pr{2};
            else
                lam = sn.rates(ist,r);
                nomD0{ist,r} = -lam; nomD1{ist,r} = lam;
            end
        end
        nomPie{ist,r} = map_pie({nomD0{ist,r}, nomD1{ist,r}});
    end
end

% ---- events -------------------------------------------------------------
% desc rows: [kind i c k j n l ip], kind 1=A0 2=A1 3=S 4=D 5=R
[A, desc] = local_events(sn, ublocks, xblocks, dim, K, rt);

% ---- horizon ------------------------------------------------------------
t0 = options.timespan(1);
tend = options.timespan(2);
unbounded = ~isfinite(tend);
if ~isfinite(t0), t0 = 0; end
period = 0;
for e = 1:numel(sched)
    if sched(e).cyclic
        period = max(period, sched(e).breakpoints(end) - sched(e).breakpoints(1));
    end
end
if unbounded
    finiteRates = sn.rates(isfinite(sn.rates) & sn.rates > 0);
    if isempty(finiteRates), slow = 1; else, slow = min(finiteRates); end
    tend = t0 + max(10, 30/slow);
    if period > 0
        tend = max(tend, t0 + 10*period);
    end
end

% ---- initial condition --------------------------------------------------
q0 = zeros(dim,1);
Sigma0 = zeros(dim,dim);
for b = 1:size(ublocks,1)
    ist = ublocks(b,1); r = ublocks(b,2); o = ublocks(b,3); h = ublocks(b,4);
    [D0t, D1t] = local_pair_at(sn, sched, nomD0, nomD1, ist, r, t0);
    Q = D0t + D1t;
    theta = [Q'; ones(1,h)] \ [zeros(h,1); 1];
    theta = max(theta, 0); theta = theta / sum(theta);
    q0(o+1:o+h) = theta;
    % The arrival phase is drawn from the stationary vector rather than known,
    % so its indicator carries covariance diag(u0) - u0*u0'. Zero here would
    % assert a known initial phase and understate the variance early on.
    Sigma0(o+1:o+h, o+1:o+h) = diag(theta) - theta*theta';
end
% NOT options.init_sol: solver_fluid_analyzer fills that from
% solver_fluid_initsol, which is laid out for the CLOSING state vector. This
% method has its own layout (u-blocks before x-blocks, no returning mass), so
% consuming it silently zeroed the source phase mass and the whole network.
%
% A WRONG-SIZED SEED IS REFUSED, not ignored. Dropping it would integrate from
% the default initial condition under the caller's name and return a plausible
% trajectory for a model the caller did not ask about.
if isfield(options,'config') && isfield(options.config,'kp_init_sol') ...
        && ~isempty(options.config.kp_init_sol)
    if numel(options.config.kp_init_sol) ~= dim
        line_error(mfilename, sprintf(['options.config.kp_init_sol has %d entries but the ''kp'' ' ...
            'state vector of this model has %d, laid out station-major over the (station,class) ' ...
            'blocks. It is NOT laid out ' ...
            'like options.init_sol.'], numel(options.config.kp_init_sol), dim));
    end
    q0 = options.config.kp_init_sol(:);
end
% Companion seed for the covariance. A caller that carries a DISTRIBUTION across
% a handoff -- the ENV blend is the case this exists for -- supplies the second
% moment beside the mean, so the next stage does not restart from a point mass
% it never had. Same layout as kp_init_sol, i.e. dim-by-dim over the u-blocks
% and x-blocks in that order.
if isfield(options,'config') && isfield(options.config,'init_cov') ...
        && ~isempty(options.config.init_cov)
    initCov = options.config.init_cov;
    if ~isequal(size(initCov), [dim dim])
        line_error(mfilename, sprintf(['options.config.init_cov is %dx%d but the ''kp'' state ' ...
            'vector of this model has %d entries, so the covariance must be %dx%d.'], ...
            size(initCov,1), size(initCov,2), dim, dim, dim));
    end
    % Loose enough for the rounding of a covariance that was itself integrated,
    % tight enough to catch a matrix that is simply not one.
    if norm(initCov - initCov', 'fro') > 1e-6 * max(1, norm(initCov, 'fro'))
        line_error(mfilename, 'options.config.init_cov must be symmetric.');
    end
    Sigma0 = initCov;
end

% ---- integrate ----------------------------------------------------------
odeh = @(tt, z) local_rhs(tt, z, A, desc, sn, sched, nomD0, nomD1, nomPie, ...
    ublocks, xblocks, dim, K, rt);
tol = options.tol;
maxstep = (tend - t0)/10;
if period > 0
    widths = [];
    for e = 1:numel(sched)
        widths = [widths, diff(sched(e).breakpoints)]; %#ok<AGROW>
    end
    maxstep = min(maxstep, min(widths)/4);
end
odeopt = odeset('AbsTol', tol*1e-3, 'RelTol', tol, 'MaxStep', maxstep);
z0 = [q0; Sigma0(:)];
solStruct = ode15s(odeh, [t0 tend], z0, odeopt);
t = solStruct.x(:);
% The period average below is a QUADRATURE, and the integrator picks its grid for
% accuracy of the SOLUTION: at a tight tol it takes long steps through the smooth
% stretches and the trapezoid over them loses more than the integration gained.
% Refine with a uniform mesh over the averaging window, and bracket each schedule
% boundary so that no trapezoid interval straddles the jump in the arrival rate.
if unbounded && period > 0
    lo = max(t0, tend - period);
    refine = linspace(lo, tend, 2001).';
    bounds = lo;
    for e = 1:numel(sched)
        bp = sched(e).breakpoints;
        per = bp(end) - bp(1);
        if sched(e).cyclic && per > 0
            kmax = ceil((tend - lo)/per) + 2;
            for kk = -1:kmax
                bounds = [bounds; (bp(:) + kk*per)]; %#ok<AGROW>
            end
        else
            bounds = [bounds; bp(:)]; %#ok<AGROW>
        end
    end
    bounds = bounds(bounds > lo & bounds < tend);
    eps_b = max(1e-9, 1e-7*(tend - lo));
    brackets = [bounds - eps_b; bounds; bounds + eps_b];
    brackets = brackets(brackets > lo & brackets < tend);
    t = unique([t; refine; brackets]);
end
Z = deval(solStruct, t).';
t = t(:);
Qtraj = Z(:, 1:dim).';
Straj = reshape(Z(:, dim+1:end).', dim, dim, numel(t));

% ---- metrics ------------------------------------------------------------
% A cyclic schedule has no fixed point, so a steady-state request is answered by
% the time average over the last full period of the periodic regime; the value at
% tend would be an arbitrary point of the cycle and source and station throughput
% would disagree.
if unbounded && period > 0
    window = [max(t0, tend - period), tend];
else
    window = [];
end

QN = zeros(M,K); UN = zeros(M,K); RN = zeros(M,K); TN = zeros(M,K);
QNt = cell(M,K); UNt = cell(M,K); TNt = cell(M,K); QVart = cell(M,K);
for ist = 1:M
    for r = 1:K
        QNt{ist,r} = zeros(numel(t),1);
        UNt{ist,r} = zeros(numel(t),1);
        TNt{ist,r} = zeros(numel(t),1);
        QVart{ist,r} = zeros(numel(t),1);
    end
end

for b = 1:size(xblocks,1)
    ist = xblocks(b,1); r = xblocks(b,2); o = xblocks(b,3); h = xblocks(b,4);
    qser = sum(Qtraj(o+1:o+h, :), 1).';
    vser = zeros(numel(t),1);
    tser = zeros(numel(t),1);
    for n = 1:numel(t)
        blk = Straj(o+1:o+h, o+1:o+h, n);
        vser(n) = sum(blk(:));
        [~, D1t] = local_pair_at(sn, sched, nomD0, nomD1, ist, r, t(n));
        scale = local_capacity(sn, xblocks, Qtraj(:,n), ist);
        tser(n) = sum(D1t, 2).' * max(Qtraj(o+1:o+h, n), 0) * scale;
    end
    QNt{ist,r} = qser; QVart{ist,r} = vser; TNt{ist,r} = tser;
    if sn.sched(ist) == SchedStrategy.INF || ~isfinite(sn.nservers(ist))
        UNt{ist,r} = qser;
    else
        UNt{ist,r} = min(qser, sn.nservers(ist)) / sn.nservers(ist);
    end
    QN(ist,r) = local_summarise(qser, t, window);
    UN(ist,r) = local_summarise(UNt{ist,r}, t, window);
    TN(ist,r) = local_summarise(tser, t, window);
    if TN(ist,r) > 0
        RN(ist,r) = QN(ist,r) / TN(ist,r);
    end
end

for b = 1:size(ublocks,1)
    ist = ublocks(b,1); r = ublocks(b,2); o = ublocks(b,3); h = ublocks(b,4);
    aser = zeros(numel(t),1);
    for n = 1:numel(t)
        [~, D1t] = local_pair_at(sn, sched, nomD0, nomD1, ist, r, t(n));
        aser(n) = sum(D1t, 2).' * max(Qtraj(o+1:o+h, n), 0);
    end
    TNt{ist,r} = aser;
    TN(ist,r) = local_summarise(aser, t, window);
end

xvec_it = {Qtraj(:,end).'};
xvec_t = Qtraj.';
iters = 1;
Sigmat = Straj;
runtime = toc(runtime_start);
end

% -------------------------------------------------------------------------
function v = local_summarise(series, t, window)
if isempty(window)
    v = series(end);
    return
end
mask = t >= window(1) & t <= window(2);
if nnz(mask) < 2
    v = series(end);
else
    tt = t(mask); ss = series(mask);
    v = trapz(tt, ss) / (tt(end) - tt(1));
end
end

function scale = local_capacity(sn, xblocks, q, ist)
% Server-capacity factor min(n,c)/n shared by every class at station IST.
if sn.sched(ist) == SchedStrategy.INF || ~isfinite(sn.nservers(ist))
    scale = 1;
    return
end
ni = 0;
for b = 1:size(xblocks,1)
    if xblocks(b,1) == ist
        o = xblocks(b,3); h = xblocks(b,4);
        ni = ni + sum(max(q(o+1:o+h), 0));
    end
end
if ni <= sn.nservers(ist)
    scale = 1;
else
    scale = sn.nservers(ist) / ni;
end
end

function [isSched, D0bar, D1bar, bp, segD0, segD1, cyclic] = local_schedule(sn, ist, r)
% The piecewise-constant (D0, D1) schedule of station IST class R, together with
% its width-weighted nominal. ISSCHED is false when the process carries none.
%
% AN NHPP IS A ONE-PHASE MAPt and is lifted to that form here. Its sn.proc slot
% is {breakpoints, rates, cyclic} rather than {breakpoints, D0, D1, cyclic}, so
% segment k contributes the 1-by-1 pair (-lambda_k, lambda_k). Reading it as its
% mean rate instead would integrate a CONSTANT arrival stream and return a flat
% trajectory under a method whose whole subject is the time-varying limit.
isSched = false;
D0bar = []; D1bar = [];
bp = []; segD0 = {}; segD1 = {}; cyclic = false;
pid = sn.procid(ist,r);
if pid == ProcessType.MAPT || pid == ProcessType.PHT
    [D0bar, D1bar, bp, segD0, segD1, cyclic] = sn_schedule_nominal(sn, ist, r);
    isSched = true;
elseif pid == ProcessType.NHPP
    slot = sn.proc{ist}{r};
    bp = slot{1}(:).';
    segRate = slot{2}(:).';
    cyclic = logical(slot{3});
    n = numel(segRate);
    segD0 = cell(1, n);
    segD1 = cell(1, n);
    for k = 1:n
        segD0{k} = -segRate(k);
        segD1{k} = segRate(k);
    end
    % Same width-weighted nominal SN_SCHEDULE_NOMINAL returns for a MAPt, so a
    % consumer that has no schedule entry sees the time-average arrival rate.
    widths = diff(bp);
    D1bar = sum(segRate .* widths) / sum(widths);
    D0bar = -D1bar;
    isSched = true;
end
end

function [D0t, D1t] = local_pair_at(sn, sched, nomD0, nomD1, ist, r, t)
% (D0, D1) in force at time T; the nominal when the process has no schedule.
for e = 1:numel(sched)
    if sched(e).station == ist && sched(e).class == r
        bp = sched(e).breakpoints;
        T = bp(end) - bp(1);
        offset = t - bp(1);
        if sched(e).cyclic
            offset = mod(offset, T);
        elseif offset < 0 || offset >= T
            D0t = zeros(size(sched(e).segD0{1}));
            D1t = zeros(size(sched(e).segD1{1}));
            return
        end
        pos = bp(1) + offset;
        idx = find(pos < bp(2:end), 1, 'first');
        if isempty(idx), idx = sched(e).nseg; end
        D0t = sched(e).segD0{idx};
        D1t = sched(e).segD1{idx};
        return
    end
end
D0t = nomD0{ist,r};
D1t = nomD1{ist,r};
end

function [A, desc] = local_events(sn, ublocks, xblocks, dim, K, rt)
% Enumerate the five event families of Ko-Pender (3.1)-(3.2).
% desc rows: [kind i c k j n l ip], kind 1=A0 2=A1 3=S 4=D 5=R
jumps = {};
desc = zeros(0,8);
    function push(pairs, d)
        col = zeros(dim,1);
        for z = 1:size(pairs,1)
            col(pairs(z,1)) = col(pairs(z,1)) + pairs(z,2);
        end
        jumps{end+1} = col; %#ok<AGROW>
        desc(end+1,:) = d; %#ok<AGROW>
    end
% (A0) arrival-MAP phase change without an arrival
for b = 1:size(ublocks,1)
    ist = ublocks(b,1); r = ublocks(b,2); o = ublocks(b,3); h = ublocks(b,4);
    for k = 1:h
        for j = 1:h
            if k ~= j
                push([o+k -1; o+j 1], [1 ist r k j 0 0 0]);
            end
        end
    end
end
% (A1) arrival-MAP phase change WITH an arrival, into a service phase
for b = 1:size(ublocks,1)
    ist = ublocks(b,1); r = ublocks(b,2); o = ublocks(b,3); h = ublocks(b,4);
    for d = 1:size(xblocks,1)
        n = xblocks(d,1); l = xblocks(d,2); no = xblocks(d,3); nh = xblocks(d,4);
        if rt((ist-1)*K+r, (n-1)*K+l) <= 0, continue; end
        for k = 1:h
            for j = 1:h
                for ip = 1:nh
                    push([o+k -1; o+j 1; no+ip 1], [2 ist r k j n l ip]);
                end
            end
        end
    end
end
% (S) service phase change inside a station
for b = 1:size(xblocks,1)
    ist = xblocks(b,1); r = xblocks(b,2); o = xblocks(b,3); h = xblocks(b,4);
    for p = 1:h
        for q = 1:h
            if p ~= q
                push([o+p -1; o+q 1], [3 ist r p q 0 0 0]);
            end
        end
    end
end
% (D) completion leaving the network and (R) completion routed onward
for b = 1:size(xblocks,1)
    ist = xblocks(b,1); r = xblocks(b,2); o = xblocks(b,3); h = xblocks(b,4);
    pout = local_pout(sn, xblocks, rt, ist, r, K);
    for p = 1:h
        if pout > 0
            push([o+p -1], [4 ist r p 0 0 0 0]);
        end
    end
    for d = 1:size(xblocks,1)
        n = xblocks(d,1); l = xblocks(d,2); no = xblocks(d,3); nh = xblocks(d,4);
        if rt((ist-1)*K+r, (n-1)*K+l) <= 0, continue; end
        for p = 1:h
            for ip = 1:nh
                push([o+p -1; no+ip 1], [5 ist r p 0 n l ip]);
            end
        end
    end
end
A = cell2mat(jumps);
end

function pout = local_pout(sn, xblocks, rt, ist, r, K)
% Probability that a completion at (IST,R) leaves the network: sn.rt is closed
% through the Source, so any destination that is not a service block is an exit.
pout = 0;
for j = 1:sn.nstations
    for l = 1:K
        isBlock = any(xblocks(:,1) == j & xblocks(:,2) == l);
        if ~isBlock
            pout = pout + rt((ist-1)*K+r, (j-1)*K+l);
        end
    end
end
end

function f = local_rates(t, q, desc, sn, sched, nomD0, nomD1, nomPie, ublocks, xblocks, K, rt)
f = zeros(size(desc,1),1);
for e = 1:size(desc,1)
    kind = desc(e,1); ist = desc(e,2); r = desc(e,3);
    k = desc(e,4); j = desc(e,5); n = desc(e,6); l = desc(e,7); ip = desc(e,8);
    [D0t, D1t] = local_pair_at(sn, sched, nomD0, nomD1, ist, r, t);
    switch kind
        case 1 % A0
            o = ublocks(ublocks(:,1)==ist & ublocks(:,2)==r, 3);
            f(e) = D0t(k,j) * max(q(o+k), 0);
        case 2 % A1
            o = ublocks(ublocks(:,1)==ist & ublocks(:,2)==r, 3);
            beta = nomPie{n,l};
            p = rt((ist-1)*K+r, (n-1)*K+l);
            f(e) = D1t(k,j) * p * beta(ip) * max(q(o+k), 0);
        case 3 % S
            o = xblocks(xblocks(:,1)==ist & xblocks(:,2)==r, 3);
            f(e) = D0t(k,j) * max(q(o+k), 0) * local_capacity(sn, xblocks, q, ist);
        case 4 % D
            o = xblocks(xblocks(:,1)==ist & xblocks(:,2)==r, 3);
            pout = local_pout(sn, xblocks, rt, ist, r, K);
            f(e) = sum(D1t(k,:)) * pout * max(q(o+k), 0) * local_capacity(sn, xblocks, q, ist);
        case 5 % R
            o = xblocks(xblocks(:,1)==ist & xblocks(:,2)==r, 3);
            beta = nomPie{n,l};
            p = rt((ist-1)*K+r, (n-1)*K+l);
            f(e) = sum(D1t(k,:)) * p * beta(ip) * max(q(o+k), 0) * local_capacity(sn, xblocks, q, ist);
    end
end
end

function dz = local_rhs(t, z, A, desc, sn, sched, nomD0, nomD1, nomPie, ublocks, xblocks, dim, K, rt)
q = z(1:dim);
Sigma = reshape(z(dim+1:end), dim, dim);
f = local_rates(t, q, desc, sn, sched, nomD0, nomD1, nomPie, ublocks, xblocks, K, rt);
dq = A * f;
% Jacobian by central differences on the assembled rates, so every capacity term
% is differentiated consistently with the drift actually integrated.
J = zeros(dim, dim);
hstep = 1e-6 * max(1, max(abs(q)));
for m = 1:dim
    qp = q; qp(m) = qp(m) + hstep;
    qm = q; qm(m) = qm(m) - hstep;
    fp = local_rates(t, qp, desc, sn, sched, nomD0, nomD1, nomPie, ublocks, xblocks, K, rt);
    fm = local_rates(t, qm, desc, sn, sched, nomD0, nomD1, nomPie, ublocks, xblocks, K, rt);
    J(:,m) = A * (fp - fm) / (2*hstep);
end
G = A * diag(f) * A.';
dS = J*Sigma + Sigma*J.' + G;
dz = [dq; dS(:)];
end
