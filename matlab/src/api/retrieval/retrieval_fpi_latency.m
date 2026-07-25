%{ @file retrieval_fpi_latency.m
 %  @brief FPI-based approximation of delayed-hit count and expected latency
 %
 %  @author LINE Development Team
%}

%{
 % @brief Approximates the mean delayed-hit count d_i and expected latency Z
 %
 % @details
 % Implements the latency-approximation algorithm of the paper (Appendix
 % "sec:fpi-di", Algorithm "Approximation procedure of the latency Z"), based on
 % Theorem "thm:di":
 %
 %   d_i = phi_i * lambda_i * E0[F_i^2] / (2 E0[F_i])
 %
 % where F_i is the fetch-period duration of item i and E0[.] is the Palm
 % expectation conditioned on a miss. The fetch moments are obtained from a
 % reduced absorbing CTMC of item i's visits to the retrieval stations:
 %
 %   E0[F_i^k] = k! * pi_e^{(i)} (-D0^{(i)})^{-k} e        (eq. moments)
 %
 % Steps (per the algorithm):
 %   1. FPI on the full system -> phi_i, pi_{i,0}            (retrieval_fpi)
 %   2. For each i: FPI without item i -> phi^{(i)}_{s,k};
 %      PS-station occupancy phi~_s^{(i)} = sum_{k!=i} phi^{(i)}_{s,k}
 %   3. Reduced CTMC for item i: one block of PH phases per station; shared
 %      stations (PS/SIRO/FCFS/LCFSPR) have their rates slowed by 1/(1+phi~_s^{(i)})
 %      while IS stations stay independent; routing per R; absorption = return to
 %      cache. D0^{(i)} = transient subgenerator, pi_e^{(i)} = entry distribution
 %      (R(outside->s) * PH-initial).
 %   4. Moments via eq. moments -> d_i via Theorem thm:di.
 %   5. Z = sum_i(phi_i + d_i) / sum_i lambda_i(phi_i + pi_{i,0})   (eq. latency tot)
 %
 % Routing convention: R is (S+1) x (S+1) x n; index 1 = outside (entry on miss /
 % return to cache on completion), indices 2..S+1 = retrieval stations 1..S;
 % R(a,b,i) is the routing probability from a to b for item i (rows sum to 1).
 %
 % @par Syntax:
 % @code
 % [Z, d, phi, pi0] = retrieval_fpi_latency(m, lambda, gamma, alpha, T, R, station_type)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>m<td>Cache list capacities (1 x h)
 % <tr><td>lambda<td>Per-item arrival rates (1 x n)
 % <tr><td>gamma<td>Access factors gamma(i,j)=gamma_{i,j}, (n x h)
 % <tr><td>alpha<td>Cell(1,S); alpha{s}(1,:,i) PH entry vector of item i at station s
 % <tr><td>T<td>Cell(1,S); T{s}(:,:,i) PH subgenerator of item i at station s
 % <tr><td>R<td>(S+1) x (S+1) x n routing matrices (index 1 = outside, 2..S+1 = stations)
 % <tr><td>station_type<td>(1 x S) string/cellstr per station, one of "IS", "PS", "SIRO", "FCFS", "LCFSPR". PS and LCFSPR (symmetric/insensitive BCMP disciplines) admit general phase-type service and class-dependent rates. SIRO and FCFS require exponential (single-phase) service with identical per-class rates; otherwise an error is raised.
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Z<td>Expected latency of the delayed-hit system
 % <tr><td>d<td>Mean number of delayed hits awaiting fetch, per item (1 x n)
 % <tr><td>phi<td>Delayed-hit ratio phi_i per item (1 x n)
 % <tr><td>pi0<td>Miss ratio pi_{i,0} per item (1 x n)
 % </table>
%}
function [Z,d,phi,pi0]=retrieval_fpi_latency(m,lambda,gamma,alpha,T,R,station_type)
n = numel(lambda);
m = m(:).';
S = numel(alpha);
fsz = zeros(1,S);
for s = 1:S
    fsz(s) = size(T{s},1);
end
station_type = string(station_type);
% see _kb/09-ldes-and-cache.md (delayed-hit latency: station-type equivalences)
unsupported = ~((station_type == "IS") | (station_type == "PS") | (station_type == "SIRO") | (station_type == "FCFS") | (station_type == "LCFSPR"));
if any(unsupported)
    bad = strjoin(cellstr(unique(station_type(unsupported))), ', ');
    line_error(mfilename, sprintf(['retrieval_fpi_latency supports only IS, PS, SIRO, FCFS and ' ...
        'LCFSPR retrieval stations; unsupported station type(s): %s'], bad));
end
for s = find(station_type == "SIRO" | station_type == "FCFS")
    if fsz(s) > 1
        line_error(mfilename, sprintf(['retrieval_fpi_latency supports SIRO/FCFS retrieval stations ' ...
            'only with exponential (single-phase) service; station %d has phase-type service.'], s));
    end
end
% SIRO and FCFS reduce to the PS single-exponential sojorn only with
% class-independent service rates: every item must share the same mean service
% time at such a station. (LCFSPR is exempt -- it is insensitive.)
for s = find(station_type == "FCFS" | station_type == "SIRO")
    taus = zeros(1, n);
    for i = 1:n
        taus(i) = -alpha{s}(:,:,i) / T{s}(:,:,i) * ones(fsz(s), 1);
    end
    if max(taus) - min(taus) > 1e-9 * max(taus)
        line_error(mfilename, sprintf(['retrieval_fpi_latency requires class-independent (identical) ' ...
            'mean service rates at SIRO/FCFS station %d.'], s));
    end
end
isIdx = find(station_type == "IS");
psIdx = find(station_type == "PS" | station_type == "SIRO" | station_type == "FCFS" | station_type == "LCFSPR");
r = numel(psIdx);

% --- per-item fetching demands eta_{s,i} = visits_{s,i} * mean PH service time ---
% eta_fpi(i, 1) = sum over IS stations; eta_fpi(i, 1+p) = PS station psIdx(p)
eta_fpi = zeros(n, r+1);
for i = 1:n
    Ri = R(:,:,i);
    a = Ri(1, 2:S+1);              % outside -> station entry probs (1 x S)
    P = Ri(2:S+1, 2:S+1);         % station -> station
    visits = a / (eye(S) - P);    % expected visits per fetch (1 x S)
    tau = zeros(1,S);
    for s = 1:S
        al = alpha{s}(:,:,i);
        Tm = T{s}(:,:,i);
        tau(s) = -al / Tm * ones(fsz(s),1);   % mean PH service time
    end
    eta_s = visits .* tau;        % eta_{s,i} per station
    eta_fpi(i,1) = sum(eta_s(isIdx));
    for p = 1:r
        eta_fpi(i,1+p) = eta_s(psIdx(p));
    end
end

% --- step 1: FPI on the full system ---
[pi0, ~, pdh] = retrieval_fpi(m, lambda, eta_fpi, gamma);
phi = sum(pdh, 1);                % phi_i = sum_s phi_{s,i}

d = zeros(1,n);
for i = 1:n
    % --- step 2: FPI without item i -> PS-station occupancy phi~_s^{(i)} ---
    keep = [1:i-1, i+1:n];
    [~, ~, pdh_i] = retrieval_fpi(m, lambda(keep), eta_fpi(keep,:), gamma(keep,:));
    phitilde = zeros(1,S);                  % occupancy per actual station
    for p = 1:r
        phitilde(psIdx(p)) = sum(pdh_i(1+p, :));
    end

    % --- step 3: reduced absorbing CTMC for item i's fetch ---
    Phi = sum(fsz);
    off = [0, cumsum(fsz(1:end-1))];
    D0 = zeros(Phi, Phi);
    pe = zeros(1, Phi);
    Ri = R(:,:,i);
    for s = 1:S
        rows = off(s)+(1:fsz(s));
        if station_type(s) == "PS" || station_type(s) == "SIRO" || station_type(s) == "FCFS" || station_type(s) == "LCFSPR"
            scale = 1/(1 + phitilde(s));     % PS/SIRO/FCFS/LCFSPR mean-field sharing slowdown
        else
            scale = 1;                       % IS: independent
        end
        blk = scale * T{s}(:,:,i);
        D0(rows, rows) = D0(rows, rows) + blk;
        compl = -blk * ones(fsz(s),1);       % completion-rate vector at station s
        for sp = 1:S
            cols = off(sp)+(1:fsz(sp));
            D0(rows, cols) = D0(rows, cols) + compl * Ri(s+1, sp+1) * alpha{sp}(:,:,i);
        end
        pe(rows) = Ri(1, s+1) * alpha{s}(:,:,i);   % entry distribution
    end

    % --- step 4: moments E0[F_i], E0[F_i^2] (eq. moments) and d_i (thm:di) ---
    e = ones(Phi,1);
    A = -D0;
    M1 = pe * (A \ e);
    M2 = 2 * pe * (A \ (A \ e));
    d(i) = phi(i) * lambda(i) * M2 / (2*M1);
end

% --- step 5: expected latency Z (eq. latency tot) ---
Z = sum(phi + d) / sum(lambda(:).' .* (phi + pi0));
end
