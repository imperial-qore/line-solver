%{ @file mg1_dt_queue.m
 %  @brief Discrete-time single-server queue with batch DMAP arrivals
 %
 %  @author LINE Development Team
%}

%{
 % @brief Discrete-time single-server queue with batch DMAP arrivals
 %
 % @details
 % Solves the DBMAP/DMAP/1 queue on a slotted time scale under the late
 % arrival system with delayed access (LAS-DA): within a slot the service
 % completion resolves first, arrivals are appended at the end of the slot and
 % cannot enter service before the next slot, and the level is read after both.
 % This is the convention of Q_DT_MAP_MAP_1, whose QBD blocks this function
 % reproduces for a single arrival per slot, and of the LDES slotted engine.
 %
 % The chain is M/G/1-type because a slot may deliver a batch: with arrival
 % matrices A_k and service pair (S0,S1),
 %   A^(-1) = kron(A_0,S1),  A^(k) = kron(A_k,S0) + kron(A_{k+1},S1),
 %   B^(k)  = kron(A_k,I),
 % the boundary row holding the empty system, where no service runs.
 %
 % @par Syntax:
 % @code
 % [QN, UN, TN, ql] = mg1_dt_queue(ARV, SVC)
 % [QN, UN, TN, ql, DEP] = mg1_dt_queue(ARV, SVC, options)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>ARV<td>Cell {A_0, A_1, ...} of arrival batch matrices
 % <tr><td>SVC<td>Cell {S0, S1} of the service DMAP
 % <tr><td>options<td>(Optional) reads config.space_max and config.dt_maxlevel
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>QN<td>Mean number in system at slot boundaries (LAS-DA)
 % <tr><td>UN<td>Utilization, i.e. fraction of slots with the server busy
 % <tr><td>TN<td>Throughput in jobs per slot
 % <tr><td>ql<td>Queue length pmf, ql(i) = Prob[i-1 in system]
 % <tr><td>DEP<td>Departure process as a DMAP {D0,D1} on the truncated chain
 % </table>
%}
function [QN, UN, TN, ql, DEP] = mg1_dt_queue(ARV, SVC, options)

if nargin < 3
    options = struct();
end
maxNumComp = 500;
if isfield(options, 'config') && isfield(options.config, 'dt_maxlevel') && ~isempty(options.config.dt_maxlevel)
    maxNumComp = options.config.dt_maxlevel;
end

S0 = SVC{1};
S1 = SVC{2};
ms = size(S0, 1);
ma = size(ARV{1}, 1);
K = length(ARV) - 1;
Ims = eye(ms);

lambda = dmap_lambda(ARV);
mu = dmap_lambda(SVC);
if lambda >= mu
    line_error(mfilename, ['The discrete-time load %g of the station is not below one ' ...
        '(%g arrivals per slot against %g completions per busy slot).'], lambda/mu, lambda, mu);
end

% M/G/1-type blocks, level = number in system
m = ma * ms;
Acat = zeros(m, m * (K + 2));
Acat(:, 1:m) = kron(ARV{1}, S1);                       % A^(-1), a departure and no arrival
for k = 0:K
    blk = kron(ARV{k+1}, S0);
    if k + 1 <= K
        blk = blk + kron(ARV{k+2}, S1);
    end
    Acat(:, (k+1)*m + 1 : (k+2)*m) = blk;
end
Bcat = zeros(m, m * (K + 1));
for k = 0:K
    Bcat(:, k*m + 1 : (k+1)*m) = kron(ARV{k+1}, Ims);  % empty system: no service runs
end

G = MG1_CR(Acat);
pivec = MG1_pi(Bcat, Acat, G, 'MaxNumComp', maxNumComp);

nlev = length(pivec) / m;
ql = zeros(1, nlev);
for i = 1:nlev
    ql(i) = sum(pivec((i-1)*m + 1 : i*m));
end
ql = ql / sum(ql);

QN = (0:(nlev-1)) * ql(:);
UN = 1 - ql(1);
TN = lambda;

if nargout < 5
    return;
end

% Departure process of the truncated chain: levels 0..L with arrivals that
% would cross L held at L. Truncation is a level cut, not a rate change, so
% the departure rate it reports stays within the mass left above L.
tolMass = GlobalConstants.Zero;
L = find(cumsum(ql) > 1 - max(tolMass, 1e-10), 1, 'first');
if isempty(L)
    L = nlev;
end
L = max(1, L - 1);

nstates = (L + 1) * m;
D0 = zeros(nstates, nstates);
D1 = zeros(nstates, nstates);
for i = 0:L
    rows = i*m + 1 : (i+1)*m;
    for k = 0:K
        if i == 0
            % empty system: the slot carries no completion
            tgt = min(L, k);
            cols = tgt*m + 1 : (tgt+1)*m;
            D0(rows, cols) = D0(rows, cols) + kron(ARV{k+1}, Ims);
        else
            tgtNo = min(L, i + k);
            colsNo = tgtNo*m + 1 : (tgtNo+1)*m;
            D0(rows, colsNo) = D0(rows, colsNo) + kron(ARV{k+1}, S0);
            tgtDep = min(L, i - 1 + k);
            colsDep = tgtDep*m + 1 : (tgtDep+1)*m;
            D1(rows, colsDep) = D1(rows, colsDep) + kron(ARV{k+1}, S1);
        end
    end
end
DEP = {D0, D1};

end
