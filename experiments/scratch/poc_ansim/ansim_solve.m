function [QN, elapsed, iters, ok, msg] = ansim_solve(cfg, mode)
% ANSIM_SOLVE  Analytical-simulation surrogate (ansim v2) for the two-class
% Delay+Queue priority model.
%
% One tagged job of the class being updated cycles through the network while
% every other job is PINNED as a self-looping customer at the station where the
% current queue-length estimate puts it. The tagged job's per-station occupancy,
% scaled by the class population, is the new estimate; iterate to a fixed point.
%
% MODE selects how the surrogate is solved:
%   'ctmcP' -- exact CTMC of the surrogate, KEEPING PSPRIO. The frozen field must
%              then be integral, so it is rounded (largest remainder).
%   'mn'    -- fluid minnormal on FRACTIONAL field populations. The fluid solver
%              cannot represent preemption, so the surrogate degrades to DPS
%              weights and inherits the dps-gap.
%   'live'  -- as 'ctmcP', except that every class STRICTLY HIGHER in priority
%              than the tagged one is kept LIVE inside the surrogate (a real
%              closed class cycling delay<->queue) instead of pinned. Under
%              preemption a pinned high-priority job never leaves the server, so
%              a frozen field starves the tagged low-priority job outright the
%              moment its rounded population at the queue reaches 1. This is the
%              design rule of the conditioning study: keep live exactly those
%              state components on which the tagged job's service rate depends
%              non-smoothly, freeze the rest into means.
%
% Returns the M-by-K mean queue-length matrix.

M = 2;
K = numel(cfg.N);
ok = true; msg = '';

QN = zeros(M,K);
QN(1,:) = cfg.N;             % start with the whole population thinking

t0 = tic;
iters = 0;
for it = 1:cfg.maxIter
    iters = it;
    QNprev = QN;
    for r = 1:K
        % Field seen by the tagged class-r job: its own class contributes
        % N_r - 1 jobs (the arrival-theorem analogue), every other class N_s.
        field = QN;
        if cfg.N(r) > 0
            field(:,r) = QN(:,r) * (cfg.N(r)-1) / cfg.N(r);
        end
        [q, sok, smsg] = ansim_tagged(cfg, r, field, mode);
        if ~sok
            QN = nan(M,K); elapsed = toc(t0); ok = false; msg = smsg; return;
        end
        QN(:,r) = cfg.N(r) * q(:);
    end
    if norm(QN(:)-QNprev(:), 1) < cfg.tol
        break;
    end
end
elapsed = toc(t0);
end

% ------------------------------------------------------------------------

function [q, ok, msg] = ansim_tagged(cfg, r, field, mode)
% Occupancy of the single tagged class-r job across the two stations.
M = 2;
K = numel(cfg.N);
q = zeros(M,1); ok = true; msg = '';

% Which classes stay dynamic inside the surrogate rather than being pinned.
live = false(1,K);
if strcmp(mode,'live')
    live = cfg.prio < cfg.prio(r);       % strictly more urgent than the tagged job
end

switch mode
    case {'ctmcP','live'}
        pop = zeros(M,K);
        for s = 1:K
            pop(:,s) = ansim_roundpreserve(field(:,s), cfg.N(s) - (s==r));
        end
    case 'mn'
        pop = field;                     % fractional field is the point of 'mn'
    otherwise
        error('unknown ansim mode %s', mode);
end

m = Network('ansim_surrogate');
nd{1} = Delay(m, 'Think');
switch mode
    case {'ctmcP','live'}, nd{2} = Queue(m, 'Q', SchedStrategy.PSPRIO);
    case 'mn',             nd{2} = Queue(m, 'Q', SchedStrategy.DPS);
end

isdps = strcmp(mode,'mn');
w = ansim_dpsweights(cfg);

% Pinned field: one self-looping class per (station, class) pair. A live class
% contributes no pinned jobs at all -- it gets a cycling closed class below.
jc = cell(M,K);
lv = cell(1,K);
for s = 1:K
    if live(s)
        lv{s} = ClosedClass(m, sprintf('L%d',s), cfg.N(s), nd{1}, cfg.prio(s));
        continue;
    end
    for i = 1:M
        if isdps, p = 0; else, p = cfg.prio(s); end
        jc{i,s} = SelfLoopingClass(m, sprintf('F%d%d',i,s), pop(i,s), nd{i}, p);
    end
end
% The one free job.
if isdps, pmov = 0; else, pmov = cfg.prio(r); end
mov = ClosedClass(m, 'Moving', 1, nd{1}, pmov);

% A pinned class is served only at the station it is pinned to.
for i = 1:M
    for s = 1:K
        if live(s), continue; end
        for j = 1:M
            if j == i
                if j == 1
                    nd{j}.setService(jc{i,s}, Exp(1/cfg.Z(s)));
                elseif isdps
                    nd{j}.setService(jc{i,s}, ansim_svc(cfg.D(s), cfg.law), w(s));
                else
                    nd{j}.setService(jc{i,s}, ansim_svc(cfg.D(s), cfg.law));
                end
            else
                nd{j}.setService(jc{i,s}, Disabled());
            end
        end
    end
end
for s = 1:K
    if ~live(s), continue; end
    nd{1}.setService(lv{s}, Exp(1/cfg.Z(s)));
    if isdps
        nd{2}.setService(lv{s}, ansim_svc(cfg.D(s), cfg.law), w(s));
    else
        nd{2}.setService(lv{s}, ansim_svc(cfg.D(s), cfg.law));
    end
end
nd{1}.setService(mov, Exp(1/cfg.Z(r)));
if isdps
    nd{2}.setService(mov, ansim_svc(cfg.D(r), cfg.law), w(r));
else
    nd{2}.setService(mov, ansim_svc(cfg.D(r), cfg.law));
end

P = m.initRoutingMatrix;
for c = 1:m.getNumberOfClasses()
    P{c} = Network.serialRouting(nd);
end
m.link(P);

movIdx = mov.index;
try
    switch mode
        case {'ctmcP','live'}
            Q = SolverCTMC(m, 'verbose', false, 'force', true).getAvgQLen();
        case 'mn'
            Q = SolverFLD(m, 'method', 'minnormal', 'verbose', false, 'force', true).getAvgQLen();
    end
catch ME
    ok = false; msg = ME.message; return;
end

q = Q(:, movIdx);
s = sum(q);
if ~isfinite(s) || s <= 0
    ok = false; msg = 'tagged occupancy did not sum to a positive value';
    return;
end
q = q / s;   % the tagged job is somewhere: renormalize away solver drift
end
