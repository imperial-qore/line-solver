function desc = mdd_descriptor(mu, P, servers, N, options)
% DESC = MDD_DESCRIPTOR(MU, P, SERVERS, N)
% DESC = MDD_DESCRIPTOR(MU, P, SERVERS, N, OPTIONS)
% Kronecker rate descriptor of a single-class closed queueing network, for the
% Miner-Ciardo-Donatelli approximate-aggregation solver (mdd_mcd), after
% A.S. Miner, G. Ciardo, S. Donatelli, "Using the exact state space of a
% Markov model to compute approximate stationary measures", SIGMETRICS 2000.
%
% The transition rate matrix is expressed compositionally as
%   R = sum_e ( kron_{k} W_k^e )   restricted to the reachable set,
% with W_k^e[i_k,j_k] = lambda_k^e[i_k] * Prob_k^e(i_k,j_k) (Eq. 1). Each level
% k is a station and each event e is a completion at station a routed to b.
%
% -- Exponential stations
% The local state is the population alone:
%   W_a^e[i,i-1] = mu(a)*min(i,servers(a))*P(a,b)   (i >= 1)   -- departure
%   W_b^e[i,i+1] = 1                                (i <= N-1) -- arrival
%   W_l^e        = I                                (l ~= a,b)
%
% -- Phase-type stations
% The local state is the PAIR (population, phase of the job in service),
% encoded in one level rather than two. Splitting them does not work: on a
% completion routed into station b the phase at b restarts only when b was
% empty, a joint condition on b's two components, which is not a product of
% per-level terms. Merging them keeps every event local. The encoding is
%   index 0                      : station empty
%   index 1 + (n-1)*h + (a-1)    : n jobs present, job in service in phase a
% so the domain is 1 + N*h and h = 1 reproduces the exponential encoding
% index = n exactly. With exit vector t = -D0*1 and entry law pie,
%   departure  (n,a) -> (n-1,b)  at t(a)*P(a,b)*pie(b)   for n >= 2
%              (1,a) -> 0        at t(a)*P(a,b)
%   arrival    0     -> (1,b)    at pie(b)
%              (m,c) -> (m+1,c)  at 1                    for m >= 1
%   internal   (n,a) -> (n,b)    at D0(a,b)              for n >= 1, a ~= b
%
% -- Input
% MU      : 1 x K station service rates. Entry i is ignored when station i is
%           given a phase-type law through OPTIONS.proc.
% P       : K x K routing matrix (row-stochastic)
% SERVERS : 1 x K servers per station (Inf for delay/IS)
% N       : closed population
% OPTIONS : struct, fields
%             proc - 1 x K cell; entry i empty for an exponential station, or
%                    a Markovian distribution object, or {D0,D1}
%             pie  - 1 x K cell of entry laws; taken from the object, or the
%                    stationary entry law of {D0,D1}, when omitted
%             sched- 1 x K cell of SchedStrategy values. Optional, and only
%                    consulted to REJECT a phase-type law at a
%                    preemptive-resume station (see the restriction below).
%                    Supply it whenever the stations are not all
%                    non-preemptive, because the descriptor otherwise has no
%                    way to detect that case.
% -- Output
% DESC    : struct with fields
%             K, N, domain, mu, servers, P
%             nphases  - 1 x K phases per station (1 when exponential)
%             valuemap - 1 x K cell, valuemap{i}(idx+1) = jobs at station i
%             init     - initial local index per station
%             nextfun  - successor function over local indices, for mdd_reachset
%             events   - struct array, one per routing pair, with fields
%                        a, b, lev (levels touched), W (matrices at lev)
%
% -- Restrictions
% A phase-type station must be single-server. With c > 1 or an infinite server
% the local state would have to count jobs per phase rather than name one
% phase, which is a different and much larger encoding.
%
% A phase-type station must also be NON-preemptive. The composite level names
% the phase of the one job in service and restarts it at pie when the next job
% starts; under preemptive resume an arrival suspends that job and its phase
% has to be remembered, so the local state would need a stack of phases. This
% matters for LCFSPR, which is BCMP type 2 and stays product-form under general
% service: that insensitivity is real but is NOT reachable through this
% encoding. Exponential service is unaffected, preemption being immaterial by
% memorylessness. Pass OPTIONS.sched to have the case rejected rather than
% silently modelled as non-preemptive.
%
% See also: mdd_mcd, mdd_reachset, mdd_closedqn.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5 || isempty(options), options = struct(); end
if ~isfield(options, 'proc'), options.proc = {}; end
if ~isfield(options, 'pie'),  options.pie = {}; end
if ~isfield(options, 'sched'), options.sched = {}; end

K = numel(mu);
mu = mu(:)';
servers = servers(:)';

% ---- per-station service law
D0 = cell(1, K); D1 = cell(1, K); pie = cell(1, K); h = ones(1, K);
for i = 1:K
    pr = [];
    if numel(options.proc) >= i, pr = options.proc{i}; end
    if isempty(pr)
        D0{i} = -mu(i); D1{i} = mu(i); pie{i} = 1; h(i) = 1;
        continue
    end
    if iscell(pr)
        D0{i} = full(pr{1}); D1{i} = full(pr{2});
    else
        pc = pr.getProcess();
        D0{i} = full(pc{1}); D1{i} = full(pc{2});
    end
    h(i) = size(D0{i}, 1);
    if numel(options.pie) >= i && ~isempty(options.pie{i})
        pie{i} = options.pie{i}(:)';
    elseif ~iscell(pr) && ismethod(pr, 'getInitProb')
        pie{i} = pr.getInitProb();
        pie{i} = pie{i}(:)';
    else
        pie{i} = i_entrylaw(D1{i}, h(i), i);
    end
    pie{i} = pie{i} / sum(pie{i});
    if h(i) > 1 && servers(i) ~= 1
        line_error(mfilename, sprintf(['station %d has a phase-type service law and %g servers; ' ...
            'a multi-server or delay station would have to count jobs per phase rather than ' ...
            'name the phase of one job in service, which this encoding does not carry'], i, servers(i)));
    end
    % The composite level names ONE in-service phase and restarts it at pie
    % when the next job starts, i.e. NON-preemptive service. Under
    % preemptive-resume an arrival suspends the job in service and its phase
    % must be remembered, so the local state would need a STACK of phases.
    if h(i) > 1 && numel(options.sched) >= i && ~isempty(options.sched{i})
        if ismember(options.sched{i}, [SchedStrategy.LCFSPR, SchedStrategy.FCFSPR, ...
                SchedStrategy.LCFSPRPRIO, SchedStrategy.FCFSPRPRIO])
            line_error(mfilename, sprintf(['station %d combines a phase-type service law with a ' ...
                'preemptive-resume discipline; the suspended jobs'' phases would have to be stacked ' ...
                'in the local state, which this encoding does not carry, and the descriptor would ' ...
                'silently model the non-preemptive chain instead'], i));
        end
        % Under a shared server every job present is in service and holds its
        % own phase, so naming one in-service phase is the wrong state.
        if ismember(options.sched{i}, [SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS])
            line_error(mfilename, sprintf(['station %d combines a phase-type service law with a ' ...
                'shared-server discipline; every job present is in service and holds its own ' ...
                'phase, which this encoding does not carry. Use MDD_PS, whose local state is the ' ...
                'per-phase count vector.'], i));
        end
    end
end

d = 1 + N * h;                              % local domain per station
desc.K = K;
desc.N = N;
desc.domain = d;
desc.mu = mu;
desc.servers = servers;
desc.P = P;
desc.nphases = h;

% ---- index maps
desc.valuemap = cell(1, K);
for i = 1:K
    vm = zeros(1, d(i));
    for n = 1:N
        vm((1 + (n - 1) * h(i) + 1):(1 + n * h(i))) = n;
    end
    desc.valuemap{i} = vm;
end

% ---- initial state: all jobs at station 1, in its entry phase
init = zeros(1, K);
init(1) = i_idx(N, find(pie{1} > 0, 1), h(1));
desc.init = init;
desc.nextfun = @(s) i_next(s, D0, D1, pie, h, N, P, servers, desc.valuemap);

% ---- events
[aa, bb, pab] = find(P);
keep = aa ~= bb;
aa = aa(keep); bb = bb(keep); pab = pab(keep);

events = repmat(struct('a', 0, 'b', 0, 'lev', [], 'W', []), 0, 1);

% internal phase changes, one event per phase-type station
for i = 1:K
    if h(i) == 1, continue; end
    Wi = i_internal(D0{i}, h(i), N, d(i));
    if nnz(Wi) == 0, continue; end
    events(end+1, 1) = struct('a', i, 'b', i, 'lev', i, 'W', {{Wi}}); %#ok<AGROW>
end

% completions routed a -> b
for e = 1:numel(aa)
    a = aa(e); b = bb(e); pr = pab(e);
    Wa = i_departure(D0{a}, D1{a}, pie{a}, h(a), N, d(a), mu(a), servers(a), pr);
    Wb = i_arrival(pie{b}, h(b), N, d(b));
    events(end+1, 1) = struct('a', a, 'b', b, 'lev', [a b], 'W', {{Wa, Wb}}); %#ok<AGROW>
end
desc.events = events;
end

% ------------------------------------------------------------------------
function pie = i_entrylaw(D1i, h, i)
% Entry law carried by a {D0,D1} pair. For a renewal process D1 = t0*pie, so
% every row with a positive exit rate is proportional to pie. Deriving it is
% not optional: defaulting to e_1 silently replaces a hyperexponential (whose
% D0 is diagonal, so a job entering phase 1 can never leave it) by an
% exponential at the phase-1 rate.
t0 = sum(D1i, 2);
live = find(t0 > 0);
if isempty(live)
    pie = [1, zeros(1, h - 1)];
    return
end
rows = D1i(live, :) ./ t0(live);
% A non-renewal MAP restarts in a law that depends on the phase it left from,
% which one entry vector cannot express; the composite level would silently
% model the renewal process instead.
if numel(live) > 1 && max(max(abs(rows - rows(1, :)))) > 1e-9
    line_error(mfilename, sprintf(['station %d carries a service law whose restart distribution ' ...
        'depends on the completing phase (a non-renewal MAP); the local state names one entry ' ...
        'law, so this encoding cannot represent it'], i));
end
pie = rows(1, :);
end

% ------------------------------------------------------------------------
function idx = i_idx(n, a, h)
% local index of (population n, service phase a); 0 when the station is empty
if n == 0, idx = 0; else, idx = 1 + (n - 1) * h + (a - 1); end
end

% ------------------------------------------------------------------------
function [n, a] = i_decode(idx, h)
if idx == 0, n = 0; a = 0; return; end
n = floor((idx - 1) / h) + 1;
a = mod(idx - 1, h) + 1;
end

% ------------------------------------------------------------------------
function W = i_internal(D0i, h, N, d)
% phase changes that do not complete a service, at any population n >= 1
ri = []; ci = []; vv = [];
for n = 1:N
    for a = 1:h
        for b = 1:h
            if a == b || D0i(a, b) == 0, continue; end
            ri(end+1) = i_idx(n, a, h) + 1; %#ok<AGROW>
            ci(end+1) = i_idx(n, b, h) + 1; %#ok<AGROW>
            vv(end+1) = D0i(a, b);          %#ok<AGROW>
        end
    end
end
W = sparse(ri, ci, vv, d, d);
end

% ------------------------------------------------------------------------
function W = i_departure(D0i, D1i, piei, h, N, d, mui, srv, pr)
% completion at this station, routed out with probability pr
ri = []; ci = []; vv = [];
if h == 1
    % exponential: the multi-server and delay rate laws live here
    for n = 1:N
        ri(end+1) = n + 1; ci(end+1) = n;        %#ok<AGROW>
        vv(end+1) = mui * min(n, srv) * pr;      %#ok<AGROW>
    end
else
    t = sum(D1i, 2);                             % exit rate per phase
    for n = 1:N
        for a = 1:h
            if t(a) == 0, continue; end
            if n == 1
                ri(end+1) = i_idx(1, a, h) + 1; ci(end+1) = 1; %#ok<AGROW>
                vv(end+1) = t(a) * pr;                          %#ok<AGROW>
            else
                for b = 1:h
                    if piei(b) == 0, continue; end
                    ri(end+1) = i_idx(n, a, h) + 1;             %#ok<AGROW>
                    ci(end+1) = i_idx(n - 1, b, h) + 1;         %#ok<AGROW>
                    vv(end+1) = t(a) * pr * piei(b);            %#ok<AGROW>
                end
            end
        end
    end
end
W = sparse(ri, ci, vv, d, d);
end

% ------------------------------------------------------------------------
function W = i_arrival(pieb, h, N, d)
% an arrival starts service only when the station was empty
ri = []; ci = []; vv = [];
if h == 1
    for m = 0:(N - 1)
        ri(end+1) = m + 1; ci(end+1) = m + 2; vv(end+1) = 1; %#ok<AGROW>
    end
else
    for b = 1:h
        if pieb(b) == 0, continue; end
        ri(end+1) = 1; ci(end+1) = i_idx(1, b, h) + 1; vv(end+1) = pieb(b); %#ok<AGROW>
    end
    for m = 1:(N - 1)
        for c = 1:h
            ri(end+1) = i_idx(m, c, h) + 1;         %#ok<AGROW>
            ci(end+1) = i_idx(m + 1, c, h) + 1;     %#ok<AGROW>
            vv(end+1) = 1;                          %#ok<AGROW>
        end
    end
end
W = sparse(ri, ci, vv, d, d);
end

% ------------------------------------------------------------------------
function T = i_next(s, D0, D1, pie, h, N, P, servers, valuemap) %#ok<INUSD>
% successor local-index vectors of s, used to generate the reachable set
K = numel(s);
T = zeros(0, K);
for i = 1:K
    [ni, ai] = i_decode(s(i), h(i));
    if ni == 0, continue; end
    % internal phase change
    if h(i) > 1
        for b = 1:h(i)
            if b == ai || D0{i}(ai, b) == 0, continue; end
            t = s; t(i) = i_idx(ni, b, h(i));
            T(end+1, :) = t; %#ok<AGROW>
        end
    end
    % completion routed to j
    if h(i) == 1
        exits = 1;
    else
        exits = sum(D1{i}, 2); exits = exits(ai);
    end
    if exits == 0, continue; end
    for j = 1:K
        if j == i || P(i, j) <= 0, continue; end
        [nj, aj] = i_decode(s(j), h(j));
        if h(i) == 1
            newi = i_idx(ni - 1, 1, 1);
        end
        for bi = 1:h(i)
            if h(i) > 1
                if ni == 1
                    newi = 0;
                elseif pie{i}(bi) == 0
                    continue
                else
                    newi = i_idx(ni - 1, bi, h(i));
                end
            elseif bi > 1
                continue
            end
            for bj = 1:h(j)
                if nj == 0
                    if pie{j}(bj) == 0, continue; end
                    newj = i_idx(1, bj, h(j));
                elseif bj > 1
                    continue
                else
                    newj = i_idx(nj + 1, aj, h(j));
                end
                t = s; t(i) = newi; t(j) = newj;
                T(end+1, :) = t; %#ok<AGROW>
            end
            if ni == 1 && h(i) > 1, break; end
        end
    end
end
end
