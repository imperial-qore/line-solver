function desc = mdd_ps(mu, P, servers, N, options)
% DESC = MDD_PS(MU, P, SERVERS, N, OPTIONS)
% Kronecker rate descriptor of a single-class closed queueing network whose
% stations are processor-sharing or infinite-server, with phase-type service.
%
% Under processor sharing every job at a station is in service at once, each
% holding its own phase, so naming a single in-service phase (what
% MDD_DESCRIPTOR does, which is non-preemptive semantics) cannot represent the
% state. The local state here is instead the PER-PHASE COUNT vector
% v = (v_1,...,v_h), v_a jobs in phase a, with n = sum(v) jobs present. That is
% still a per-station quantity, so every event stays a product of per-level
% terms and the Kronecker form of Eq. 1 survives.
%
% With one server shared by n jobs each job advances at rate 1/n, so from local
% state v with n = sum(v):
%   internal   v -> v - e_a + e_b   at v_a * D0(a,b) / n     (a ~= b)
%   departure  v -> v - e_a         at v_a * t(a) / n * P(i,j)
%   arrival    v -> v + e_b         at pie(b)
% An infinite-server (delay) station is the same without the 1/n scaling. For
% h = 1 the departure rate collapses to n*mu/n = mu at PS and to n*mu at IS,
% reproducing the usual single-server and delay rate laws.
%
% The local domain is the number of compositions of 0..N over h phases,
% C(N+h,h), against 1+N*h for the non-preemptive encoding: the price of
% tracking every job's phase rather than one.
%
% -- Input
% MU      : 1 x K station service rates (ignored where OPTIONS.proc gives a law)
% P       : K x K routing matrix (row-stochastic)
% SERVERS : 1 x K servers per station, 1 (PS) or Inf (IS); no other value has
%           a per-phase-count encoding here
% N       : closed population
% OPTIONS : struct, fields proc and pie as in MDD_DESCRIPTOR
% -- Output
% DESC    : descriptor for MDD_MCD, with K, N, domain, mu, servers, P,
%           nphases, valuemap, init, nextfun, events, and the population
%           invariant through valuemap
%
% See also: mdd_descriptor, mdd_mcd, mdd_reachset.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5 || isempty(options), options = struct(); end
if ~isfield(options, 'proc'), options.proc = {}; end
if ~isfield(options, 'pie'),  options.pie = {}; end

K = numel(mu);
mu = mu(:)'; servers = servers(:)';

D0 = cell(1, K); D1 = cell(1, K); pie = cell(1, K); h = ones(1, K);
for i = 1:K
    if ~ismember(servers(i), [1 Inf])
        line_error(mfilename, sprintf(['station %d has %g servers; only processor sharing (1) ' ...
            'and infinite server (Inf) have a per-phase-count encoding here'], i, servers(i)));
    end
    pr = [];
    if numel(options.proc) >= i, pr = options.proc{i}; end
    if isempty(pr)
        D0{i} = -mu(i); D1{i} = mu(i); pie{i} = 1; h(i) = 1;
        continue
    end
    if iscell(pr)
        D0{i} = full(pr{1}); D1{i} = full(pr{2});
    else
        pc = pr.getProcess(); D0{i} = full(pc{1}); D1{i} = full(pc{2});
    end
    h(i) = size(D0{i}, 1);
    if numel(options.pie) >= i && ~isempty(options.pie{i})
        pie{i} = options.pie{i}(:)';
    elseif ~iscell(pr) && ismethod(pr, 'getInitProb')
        pie{i} = pr.getInitProb(); pie{i} = pie{i}(:)';
    else
        pie{i} = i_entrylaw(D1{i}, h(i), i);
    end
    pie{i} = pie{i} / sum(pie{i});
end

% ---- per-station composition state space
comp = cell(1, K); lut = cell(1, K); d = zeros(1, K);
for i = 1:K
    [comp{i}, lut{i}] = i_compstates(h(i), N);
    d(i) = size(comp{i}, 1);
end

desc.K = K; desc.N = N; desc.domain = d;
desc.mu = mu; desc.servers = servers; desc.P = P; desc.nphases = h;
desc.valuemap = cell(1, K);
for i = 1:K, desc.valuemap{i} = sum(comp{i}, 2)'; end

% ---- initial state: all jobs at station 1, entered in its entry phase
init = zeros(1, K);
for i = 1:K
    v = zeros(1, h(i));
    if i == 1, v(find(pie{1} > 0, 1)) = N; end
    init(i) = i_lookup(lut{i}, v, N) - 1;
end
desc.init = init;
desc.nextfun = @(s) i_next(s, D0, D1, pie, h, comp, lut, N, P, servers);

% ---- events
events = repmat(struct('a', 0, 'b', 0, 'lev', [], 'W', []), 0, 1);
for i = 1:K
    if h(i) == 1, continue; end
    Wi = i_internal(D0{i}, comp{i}, lut{i}, h(i), N, d(i), servers(i));
    if nnz(Wi) == 0, continue; end
    events(end+1, 1) = struct('a', i, 'b', i, 'lev', i, 'W', {{Wi}}); %#ok<AGROW>
end

[aa, bb, pab] = find(P);
keep = aa ~= bb; aa = aa(keep); bb = bb(keep); pab = pab(keep);
for e = 1:numel(aa)
    a = aa(e); b = bb(e); pr = pab(e);
    Wa = i_departure(D1{a}, comp{a}, lut{a}, h(a), N, d(a), servers(a), pr);
    Wb = i_arrival(pie{b}, comp{b}, lut{b}, h(b), N, d(b));
    events(end+1, 1) = struct('a', a, 'b', b, 'lev', [a b], 'W', {{Wa, Wb}}); %#ok<AGROW>
end
desc.events = events;
end

% ------------------------------------------------------------------------
function pie = i_entrylaw(D1i, h, i)
% Entry law carried by a {D0,D1} pair; see the twin in MDD_DESCRIPTOR.
t0 = sum(D1i, 2);
live = find(t0 > 0);
if isempty(live)
    pie = [1, zeros(1, h - 1)];
    return
end
rows = D1i(live, :) ./ t0(live);
if numel(live) > 1 && max(max(abs(rows - rows(1, :)))) > 1e-9
    line_error(mfilename, sprintf(['station %d carries a service law whose restart distribution ' ...
        'depends on the completing phase (a non-renewal MAP); the local state names one entry ' ...
        'law, so this encoding cannot represent it'], i));
end
pie = rows(1, :);
end

% ------------------------------------------------------------------------
function [C, lut] = i_compstates(h, N)
% every per-phase count vector with 0 <= sum <= N, plus a mixed-radix lookup
rows = {};
stack = {zeros(1, h)};
C = i_enum(h, N);
rows = C; %#ok<NASGU>
lut = zeros((N + 1)^h, 1);
for r = 1:size(C, 1)
    lut(i_key(C(r, :), N)) = r;
end
end

% ------------------------------------------------------------------------
function C = i_enum(h, N)
if h == 1
    C = (0:N)';
    return
end
C = zeros(0, h);
sub = i_enum(h - 1, N);
subsum = sum(sub, 2);
for v1 = 0:N
    ok = subsum <= (N - v1);
    C = [C; [v1 * ones(sum(ok), 1), sub(ok, :)]]; %#ok<AGROW>
end
end

% ------------------------------------------------------------------------
function k = i_key(v, N)
k = 1; radix = 1;
for a = 1:numel(v)
    k = k + v(a) * radix;
    radix = radix * (N + 1);
end
end

% ------------------------------------------------------------------------
function idx = i_lookup(lut, v, N)
idx = lut(i_key(v, N));
end

% ------------------------------------------------------------------------
function s = i_share(n, srv)
% service-rate scaling: 1/n shared by n jobs at PS, unscaled at IS
if n == 0, s = 0; elseif isinf(srv), s = 1; else, s = 1 / n; end
end

% ------------------------------------------------------------------------
function W = i_internal(D0i, C, lut, h, N, d, srv)
ri = []; ci = []; vv = [];
for r = 1:d
    v = C(r, :); n = sum(v);
    if n == 0, continue; end
    sc = i_share(n, srv);
    for a = 1:h
        if v(a) == 0, continue; end
        for b = 1:h
            if a == b || D0i(a, b) == 0, continue; end
            w = v; w(a) = w(a) - 1; w(b) = w(b) + 1;
            ri(end+1) = r; ci(end+1) = i_lookup(lut, w, N); %#ok<AGROW>
            vv(end+1) = v(a) * D0i(a, b) * sc;              %#ok<AGROW>
        end
    end
end
W = sparse(ri, ci, vv, d, d);
end

% ------------------------------------------------------------------------
function W = i_departure(D1i, C, lut, h, N, d, srv, pr)
t = sum(D1i, 2);
ri = []; ci = []; vv = [];
for r = 1:d
    v = C(r, :); n = sum(v);
    if n == 0, continue; end
    sc = i_share(n, srv);
    for a = 1:h
        if v(a) == 0 || t(a) == 0, continue; end
        w = v; w(a) = w(a) - 1;
        ri(end+1) = r; ci(end+1) = i_lookup(lut, w, N); %#ok<AGROW>
        vv(end+1) = v(a) * t(a) * sc * pr;              %#ok<AGROW>
    end
end
W = sparse(ri, ci, vv, d, d);
end

% ------------------------------------------------------------------------
function W = i_arrival(pieb, C, lut, h, N, d)
ri = []; ci = []; vv = [];
for r = 1:d
    v = C(r, :);
    if sum(v) >= N, continue; end
    for b = 1:h
        if pieb(b) == 0, continue; end
        w = v; w(b) = w(b) + 1;
        ri(end+1) = r; ci(end+1) = i_lookup(lut, w, N); %#ok<AGROW>
        vv(end+1) = pieb(b);                            %#ok<AGROW>
    end
end
W = sparse(ri, ci, vv, d, d);
end

% ------------------------------------------------------------------------
function T = i_next(s, D0, D1, pie, h, comp, lut, N, P, servers)
K = numel(s);
T = zeros(0, K);
for i = 1:K
    v = comp{i}(s(i) + 1, :); n = sum(v);
    if n == 0, continue; end
    % internal phase moves
    for a = 1:h(i)
        if v(a) == 0, continue; end
        for b = 1:h(i)
            if a == b || D0{i}(a, b) == 0, continue; end
            w = v; w(a) = w(a) - 1; w(b) = w(b) + 1;
            t = s; t(i) = i_lookup(lut{i}, w, N) - 1;
            T(end+1, :) = t; %#ok<AGROW>
        end
    end
    % completions routed to j
    ti = sum(D1{i}, 2);
    for a = 1:h(i)
        if v(a) == 0 || ti(a) == 0, continue; end
        w = v; w(a) = w(a) - 1;
        for j = 1:K
            if j == i || P(i, j) <= 0, continue; end
            vj = comp{j}(s(j) + 1, :);
            if sum(vj) >= N, continue; end
            for b = 1:h(j)
                if pie{j}(b) == 0, continue; end
                wj = vj; wj(b) = wj(b) + 1;
                t = s;
                t(i) = i_lookup(lut{i}, w, N) - 1;
                t(j) = i_lookup(lut{j}, wj, N) - 1;
                T(end+1, :) = t; %#ok<AGROW>
            end
        end
    end
end
end
