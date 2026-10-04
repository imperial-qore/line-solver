function bg = mam_bgchain_ctmc(Nb, STb, Pb, sched, nservers, cshare, supp, options)
% BG = MAM_BGCHAIN_CTMC(NB, STB, PB, SCHED, NSERVERS, CSHARE, SUPP, OPTIONS)
%
% Builds and solves the BACKGROUND MODULATING CHAIN of a mixed network: the
% continuous-time Markov chain of the closed-class population vector, with the
% open classes present only through their mean occupancy QOPEN.
%
% The chain carries B = 1 or 2 background classes. B = 1 is the single closed
% class of the model; B = 2 is the tagged/aggregate pair built by
% SOLVER_MAM_BGCHAIN when the model has several closed chains (class 1 is the
% tagged chain r, class 2 the flow-equivalent aggregate of the other R-1).
%
% Inputs
%   NB       (1 x B)      population of each background class
%   STB      (Mc x B)     mean service time of each background class per station
%   PB       cell(1 x B)  (Mc x Mc) row-stochastic routing matrix per class
%   SCHED    (Mc x 1)     SchedStrategy id of each station
%   NSERVERS (Mc x 1)     number of servers of each station
%   CSHARE   (Mc x Nmax+1) CSHARE(i,e+1) is the mean number of servers the
%                          closed jobs of station i hold when they number e,
%                          i.e. E[min(e+k,c) * e/(e+k)] over the open occupancy
%                          k -- min(e,c) when the station carries no open work
%   SUPP     (Mc x B)     logical; SUPP(i,b) is true when station i is on the
%                         route of background class b
%   OPTIONS               solver options (reads options.config.bgstates_max)
%
% SUPP IS NOT AN OPTIMIZATION. A closed chain that never visits a station cannot
% hold jobs there, and enumerating the population vector over the union of every
% chain's stations puts probability on configurations it can never reach. Worse,
% those configurations ABSORB: the chain's routing matrix has a zero row at a
% station it does not visit, which row-normalizes to a self-loop, so a job
% placed there never leaves. The generator turns reducible and the chain's
% population conservation silently fails -- measured, half the mass of a one-job
% chain sat in states it could not reach. Each class is therefore enumerated
% over ITS OWN stations only, which also shrinks the state space.
%
% THE BLOCK OF ONE BACKGROUND CLASS IS STATE.SPACECLOSEDSINGLE, the same lattice
% primitive SOLVER_CTMC enumerates a closed population over, and the joint space
% is their cartesian product in class-major order (STATE.CARTESIAN's convention:
% the first class is the slowest index). Only the SUPPORT differs -- the columns
% are this class's own stations rather than every station -- so the enumeration
% order and the row count are the reference's, not this file's. Do not restore a
% private composition routine here: MAM_BGCHAIN_STATES predicts the row count of
% this same primitive in closed form, and a second enumeration would be free to
% drift from it.
%
% A station holding n_{i,1} + n_{i,2} = e closed jobs serves background class b
% at rate
%
%   INF                : n_{i,b} / STB(i,b)
%   any other discipline: CSHARE(i,e+1) * (n_{i,b}/e) / STB(i,b)
%
% The second line splits the capacity the closed jobs hold over the background
% classes in proportion to their counts, which is exact under PS and is the
% standard random-order surrogate under FCFS.
%
% The open classes enter ONLY through CSHARE, which is what makes this a
% MODULATING chain rather than a joint model. Note that the exchanged quantity
% is the SHARE, already averaged over the open occupancy, and not the mean open
% occupancy itself: e/(e+k) is convex in k, so rebuilding the share from a mean
% k would bias the closed service rate downwards by Jensen's inequality, and the
% closed throughput with it. SOLVER_MAM_BGCHAIN closes the loop by reading
% CSHARE back out of the modulated QBDs, where the same expectation is taken
% against the joint law of (level, environment).
%
% Outputs (struct BG)
%   .space   (nstates x Mc x B)  population vector of every state
%   .pi      (1 x nstates)       stationary distribution
%   .Q       (nstates x nstates) generator
%   .QLen    (Mc x B)            mean queue length per station per class
%   .Tput    (Mc x B)            mean completion rate per station per class
%   .Ubusy   (Mc x B)            mean fraction of the servers held by the class
%   .totocc  (nstates x Mc)      total closed occupancy per station per state
%
% See also SOLVER_MAM_BGCHAIN, MAM_BGCHAIN_ENV.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Mc = size(STb, 1);
B  = numel(Nb);

if isfield(options, 'config') && isfield(options.config, 'bgstates_max')
    bgstates_max = options.config.bgstates_max;
else
    bgstates_max = 20000;
end

%% Per-class population state spaces, over each class's OWN stations
if nargin < 7 || isempty(supp)
    supp = true(Mc, B);
end
sp = cell(1, B);
idxb = cell(1, B);
nst = zeros(1, B);
for b = 1:B
    idxb{b} = find(supp(:, b))';
    if isempty(idxb{b})
        if Nb(b) > 0
            line_error(mfilename, sprintf(['Background class %d holds %d jobs but visits no ' ...
                'station.'], b, Nb(b)));
        end
        idxb{b} = 1;   % an empty class still needs one slot to be indexed by
    end
    comp = State.spaceClosedSingle(numel(idxb{b}), round(Nb(b)));
    sp{b} = zeros(size(comp, 1), Mc);
    sp{b}(:, idxb{b}) = comp;
    nst(b) = size(comp, 1);
end
nstates = prod(nst);
if nstates > bgstates_max
    line_error(mfilename, sprintf(['The background chain of this model has %d states, above the ' ...
        'limit of %d. The chain enumerates the closed-class population vector over the %d stations ' ...
        'the closed classes visit, so its size grows as nchoosek(N+Mc-1,Mc-1) per class, and it ' ...
        'carries %d classes. Lower options.config.bgaggr to aggregate more of the closed chains ' ...
        'into fewer classes, raise options.config.bgstates_max to solve it anyway, or reduce the ' ...
        'closed populations.'], nstates, bgstates_max, Mc, B));
end

%% Per-class transition target index: tgt{b}(s, (i-1)*Mc+j) is the state
% reached from s when one class-b job moves from station i to station j, or 0
% when station i holds no class-b job in s.
tgt = cell(1, B);
for b = 1:B
    ib = idxb{b};
    mb = numel(ib);
    comp = sp{b}(:, ib);
    radix = (Nb(b) + 1) .^ (0:mb-1);
    keys = comp * radix';
    [keysSorted, ord] = sort(keys);
    tgt{b} = zeros(nst(b), Mc * Mc);
    for ii = 1:mb
        i = ib(ii);
        movable = comp(:, ii) > 0;
        if ~any(movable)
            continue;
        end
        for jj = 1:mb
            if jj == ii
                continue;
            end
            j = ib(jj);
            cand = comp(movable, :);
            cand(:, ii) = cand(:, ii) - 1;
            cand(:, jj) = cand(:, jj) + 1;
            [tf, loc] = ismember(cand * radix', keysSorted);
            col = zeros(nst(b), 1);
            col(movable) = ord(max(loc, 1)) .* tf;
            tgt{b}(:, (i-1)*Mc + j) = col;
        end
    end
end

%% Joint state space, laid out class-major: idx = (s1-1)*nst(2) + s2
space = zeros(nstates, Mc, B);
subidx = zeros(nstates, B);
for s = 1:nstates
    rem = s - 1;
    for b = B:-1:1
        subidx(s, b) = mod(rem, nst(b)) + 1;
        rem = floor(rem / nst(b));
    end
    for b = 1:B
        space(s, :, b) = sp{b}(subidx(s, b), :);
    end
end
strideb = ones(1, B);
for b = 1:B-1
    strideb(b) = prod(nst(b+1:end));
end

%% Service rates: mu(i,b) is the rate of one class-b job in service
mu = zeros(Mc, B);
for b = 1:B
    for i = 1:Mc
        if isfinite(STb(i,b)) && STb(i,b) > 0
            mu(i,b) = 1 / STb(i,b);
        end
    end
end

%% Generator
rateFull = zeros(nstates, Mc, B);   % completion rate of class b at station i
capBusy  = zeros(nstates, Mc);      % servers busy (closed jobs only)
nnzmax = nstates * Mc * B * Mc;
rowsI = zeros(nnzmax, 1); rowsJ = zeros(nnzmax, 1); rowsV = zeros(nnzmax, 1);
nrec = 0;
for s = 1:nstates
    n = reshape(space(s, :, :), Mc, B);
    for i = 1:Mc
        eclosed = sum(n(i, :));
        if eclosed == 0
            continue;
        end
        if sched(i) == SchedStrategy.INF
            held = eclosed;
        else
            held = cshare(i, min(eclosed, size(cshare, 2) - 1) + 1);
        end
        if held <= 0
            continue;
        end
        capBusy(s, i) = held;
        for b = 1:B
            if n(i, b) == 0 || mu(i, b) == 0
                continue;
            end
            r = held * (n(i, b) / eclosed) * mu(i, b);
            rateFull(s, i, b) = r;
            for j = 1:Mc
                if j == i || Pb{b}(i, j) <= 0
                    continue;
                end
                tsub = tgt{b}(subidx(s, b), (i-1)*Mc + j);
                if tsub == 0
                    continue;
                end
                sdest = s + (tsub - subidx(s, b)) * strideb(b);
                nrec = nrec + 1;
                rowsI(nrec) = s;
                rowsJ(nrec) = sdest;
                rowsV(nrec) = r * Pb{b}(i, j);
            end
        end
    end
end

Q = sparse(rowsI(1:nrec), rowsJ(1:nrec), rowsV(1:nrec), nstates, nstates);
Q = ctmc_makeinfgen(Q);
if nstates == 1
    pi = 1;
else
    pi = ctmc_solve(Q, options);
end
pi = max(pi, 0);
pi = pi / sum(pi);

%% Aggregate measures
QLen  = zeros(Mc, B);
Tput  = zeros(Mc, B);
Ubusy = zeros(Mc, B);
for b = 1:B
    QLen(:, b) = (pi * space(:, :, b))';
    Tput(:, b) = (pi * rateFull(:, :, b))';
end
totClosed = sum(space, 3);          % (nstates x Mc)
for i = 1:Mc
    if sched(i) == SchedStrategy.INF
        Ubusy(i, :) = QLen(i, :);
    else
        occ = totClosed(:, i);
        share = zeros(nstates, B);
        nz = occ > 0;
        for b = 1:B
            share(nz, b) = space(nz, i, b) ./ occ(nz);
        end
        for b = 1:B
            Ubusy(i, b) = (pi * (capBusy(:, i) .* share(:, b))) / nservers(i);
        end
    end
end

bg = struct('space', space, 'pi', pi, 'Q', Q, 'QLen', QLen, 'Tput', Tput, ...
    'Ubusy', Ubusy, 'totocc', totClosed, 'nstates', nstates);
end
