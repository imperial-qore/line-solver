function [Q0, Q1, Q2] = ldqbd_mphc(D0, D1, alpha, c, arrRate, sf)
% LDQBD_MPHC Level-dependent QBD blocks of an M/PH/c queue, exact in the phases.
%
% [Q0, Q1, Q2] = LDQBD_MPHC(D0, D1, ALPHA, C, ARRRATE) returns the
% block-tridiagonal generator of a c-server FCFS queue whose service is the
% Markovian process (D0, D1) restarted from ALPHA, with a level-dependent
% arrival rate. Level n is the number of jobs at the station; the coordinate
% INSIDE a level is the MULTISET of the phases the min(n,c) busy servers sit
% in, so the chain is exact for phase-type service at any number of servers.
%
% [Q0, Q1, Q2] = LDQBD_MPHC(..., SF) additionally scales the station's total
% service rate at level n by SF(n) (load dependence). Each busy server then
% runs at SF(n)/min(n,c) of its nominal speed, so the aggregate capacity is
% SF(n) times nominal and SF(n) = min(n,c) reproduces the unscaled queue
% exactly. Pass [] for no scaling.
%
% @par Parameters:
% <table>
% <tr><th>Name<th>Description
% <tr><td>D0<td>service sub-generator (p x p): phase changes without completion
% <tr><td>D1<td>service completion block (p x p); D1 = (-D0*1)*ALPHA for a PH
% <tr><td>ALPHA<td>1 x p vector a server starts each new job in
% <tr><td>C<td>number of identical servers (>= 1; capped at the top level)
% <tr><td>ARRRATE<td>1 x (Nlev+1); ARRRATE(n+1) is the arrival rate out of level n
% <tr><td>SF<td>optional 1 x Nlev total-service-rate multiplier per level
% </table>
%
% @par Returns:
% <table>
% <tr><th>Name<th>Description
% <tr><td>Q0<td>{Nlev x 1}; Q0{n+1} is the upward block, level n -> n+1
% <tr><td>Q1<td>{(Nlev+1) x 1}; Q1{n+1} is the local block of level n
% <tr><td>Q2<td>{Nlev x 1}; Q2{n} is the downward block, level n -> n-1
% </table>
%
% WHY THE MULTISET. The collapsed alternative -- one PH process run at
% min(n,c) times its speed -- gets the aggregate rate right but forgets which
% phase each busy server is in, which is not a detail: it makes the c servers
% behave like one fast server whose remaining work is a single phase-type
% variable. Tracking counts rather than an ordered tuple costs
% NCHOOSEK(min(n,c)+p-1, p-1) states per level instead of p^min(n,c), because
% identical servers are exchangeable.
%
% Level sizes therefore GROW over the boundary levels 0..c and repeat above
% them: level 0 is the single empty configuration, and Q0{1}, Q2{1} are the
% rectangular blocks that join a level to a differently sized neighbour.
% LDQBD, LDQBD_R and LDQBD_PI all accept that heterogeneity.
%
% References:
%   S. Asmussen and J.R. Moller, "Calculation of the steady state waiting time
%   distribution in GI/PH/c and MAP/PH/c queues", Queueing Systems 37(1):9-29,
%   2001.
%   M. F. Neuts, "Matrix-geometric solutions in stochastic models", Johns
%   Hopkins University Press, 1981.
%
% See also ldqbd, ph_multisets, qsys_mapphc, solver_mam_ldqbd
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

D0 = full(D0);
D1 = full(D1);
alpha = full(alpha(:)).';
p = size(D0, 1);
arrRate = arrRate(:).';
Nlev = numel(arrRate) - 1;

if Nlev < 1
    line_error(mfilename, 'ldqbd_mphc needs at least one level above the empty one.');
end
if size(D0, 2) ~= p || any(size(D1) ~= [p p]) || numel(alpha) ~= p
    line_error(mfilename, 'D0, D1 and alpha must all have the same order.');
end
if nargin < 6 || isempty(sf)
    sf = [];
elseif numel(sf) < Nlev
    line_error(mfilename, 'sf must give one total-service-rate factor per level 1..Nlev.');
end

% Servers that can never be busy do not need a coordinate: above the top level
% there is nothing left to serve, so the configuration space stops there.
if ~isfinite(c)
    cmax = Nlev;
else
    cmax = min(max(1, round(c)), Nlev);
end

% One configuration set per busy-server count 0..cmax, plus an index for each.
cfg = cell(cmax + 1, 1);
pos = cell(cmax + 1, 1);
for k = 0:cmax
    cfg{k+1} = ph_multisets(p, k);
    pos{k+1} = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for r = 1:size(cfg{k+1}, 1)
        pos{k+1}(cfgkey(cfg{k+1}(r, :))) = r;
    end
end
nCfg = cellfun(@(x) size(x, 1), cfg);

% The repeating level is the widest one, so it is the size worth guarding.
maxCfg = 2000;
if nCfg(cmax + 1) > maxCfg
    line_error(mfilename, sprintf(['the exact M/PH/c chain needs nchoosek(%d+%d-1,%d-1) = %d ' ...
        'configurations per level for %d servers and %d service phases, above the %d the ' ...
        'level-by-level inverses can carry. Use fewer phases (a lower-order fit), fewer ' ...
        'servers, or SolverCTMC/SolverLDES on this model.'], ...
        cmax, p, p, nCfg(cmax + 1), cmax, p, maxCfg));
end

t = D1 * ones(p, 1);   % completion rate out of each phase, summed over targets

% Structural blocks per busy-server count. LOC holds the within-level phase
% changes (and the full outflow on its diagonal), UP the entry of a newly busy
% server, DN a completion that leaves a server idle.
LOC = cell(cmax + 1, 1);
UP = cell(cmax + 1, 1);
DN = cell(cmax + 1, 1);
for k = 0:cmax
    Ck = cfg{k+1};
    nk = nCfg(k+1);
    Lk = zeros(nk, nk);
    for row = 1:nk
        m = Ck(row, :);
        for i = 1:p
            if m(i) == 0, continue; end
            for j = 1:p
                if j == i, continue; end
                mm = m; mm(i) = mm(i) - 1; mm(j) = mm(j) + 1;
                col = pos{k+1}(cfgkey(mm));
                Lk(row, col) = Lk(row, col) + m(i) * D0(i, j);
            end
            % D0(i,i) is the total outflow of phase i, completions included
            Lk(row, row) = Lk(row, row) + m(i) * D0(i, i);
        end
    end
    LOC{k+1} = Lk;

    if k < cmax
        Uk = zeros(nk, nCfg(k+2));
        for row = 1:nk
            m = Ck(row, :);
            for j = 1:p
                mm = m; mm(j) = mm(j) + 1;
                col = pos{k+2}(cfgkey(mm));
                Uk(row, col) = Uk(row, col) + alpha(j);
            end
        end
        UP{k+1} = Uk;
    end

    if k > 0
        Dk = zeros(nk, nCfg(k));
        for row = 1:nk
            m = Ck(row, :);
            for i = 1:p
                if m(i) == 0, continue; end
                mm = m; mm(i) = mm(i) - 1;
                col = pos{k}(cfgkey(mm));
                Dk(row, col) = Dk(row, col) + m(i) * t(i);
            end
        end
        DN{k+1} = Dk;
    end
end

% A completion at a full server bank takes the next waiting job at once, so the
% server stays busy and only its phase moves: the repeating downward block.
Cc = cfg{cmax+1};
nc = nCfg(cmax+1);
CDEP = zeros(nc, nc);
for row = 1:nc
    m = Cc(row, :);
    for i = 1:p
        if m(i) == 0, continue; end
        for j = 1:p
            mm = m; mm(i) = mm(i) - 1; mm(j) = mm(j) + 1;
            col = pos{cmax+1}(cfgkey(mm));
            CDEP(row, col) = CDEP(row, col) + m(i) * D1(i, j);
        end
    end
end

% Per-server speed. Without load dependence every busy server runs at its
% nominal rate; with it, the aggregate SF(n) is shared over the busy servers.
% SF(n) == min(n,c) is passed through as exactly 1 so the unscaled chain is
% reproduced bit for bit.
speed = ones(1, Nlev);
if ~isempty(sf)
    for n = 1:Nlev
        b = min(n, cmax);
        if sf(n) ~= b
            speed(n) = sf(n) / b;
        end
    end
end

Q0 = cell(Nlev, 1);
Q1 = cell(Nlev + 1, 1);
Q2 = cell(Nlev, 1);

Q1{1} = -arrRate(1);                                  % level 0: arrivals only
for n = 1:Nlev
    b = min(n, cmax);
    Q1{n+1} = speed(n) * LOC{b+1} - arrRate(n+1) * eye(nCfg(b+1));
end

for n = 0:Nlev-1
    b = min(n, cmax);
    if n < cmax
        Q0{n+1} = arrRate(n+1) * UP{n+1};             % a free server takes the job
    else
        Q0{n+1} = arrRate(n+1) * eye(nCfg(b+1));      % the job waits, phases unchanged
    end
end

for n = 1:Nlev
    if n <= cmax
        Q2{n} = speed(n) * DN{n+1};                   % the server falls idle
    else
        Q2{n} = speed(n) * CDEP;                      % the server takes the next job
    end
end

end

function k = cfgkey(m)
% Char key of a configuration row, for the index maps above.
k = sprintf('%d,', m);
end
