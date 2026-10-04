%{
 % @brief Majumdar-Woodside robust box bounds on throughput for closed
 %        multiclass queueing networks with mixed scheduling disciplines.
 %
 % @details
 % Computes distribution-insensitive (NBUE) upper and lower bounds on the
 % per-class system throughput of a closed multiclass queueing network, as
 % defined by Majumdar and Woodside, "Robust bounds and throughput guarantees
 % for closed multiclass queueing networks", Performance Evaluation 32 (1998)
 % 101-136. The upper bound intersects the no-contention bound (eq. 2) with
 % the utilization-based bound (eq. 3) and is independent of the scheduling
 % discipline. The lower bound is the multiclass throughput guarantee of
 % Theorem 2 (eq. 15): X_c >= N_c / (Z_c + sum_k V_kc (S_kc + d_kc^+)), where
 % the per-visit queueing delay bound d_kc^+ depends on the discipline at
 % station k -- FIFO (Theorem 1 / eq. 6, via Lemma 1), processor sharing
 % (Lemma 2), preemptive priority (Lemma 3) and non-preemptive priority
 % (Lemmas 4-5). The coupled inequalities are resolved by the interval-
 % narrowing fixed point that reproduces the BNR-Prolog robust box bounds;
 % for a single FIFO class it reduces to the Muntz-Wong asymptotic bounds.
 %
 % Bounds are insensitive to the service-time distributions (only NBUE is
 % assumed) and to routing dependencies; only the mean visits V, mean service
 % demands S, populations N, think times Z, per-station disciplines and per-
 % class priorities are required. The think time Z aggregates the pure-delay
 % (infinite-server) stations; only queueing stations are passed in V,S.
 %
 % @par Syntax:
 % @code
 % [Xlo,Xup,Wlo] = pfqn_mwrbb(V,S,N,Z,sched,prio)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>V<td>(K x C) mean visits of class c at queueing station k
 % <tr><td>S<td>(K x C) mean service demand per visit of class c at station k
 % <tr><td>N<td>(1 x C) population of class c
 % <tr><td>Z<td>(1 x C) think time (pure delay) of class c
 % <tr><td>sched<td>(K x 1) discipline code per station: 0=FIFO (default),
 %                  1=PS, 2=non-preemptive priority, 3=preemptive priority,
 %                  4=ABA full-contention (discipline-independent, P_cm=1)
 % <tr><td>prio<td>(1 x C) class priority, lower value = higher priority
 %                 (default: all equal). Used only at priority stations.
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Xlo<td>(1 x C) lower bound on class throughput (Theorem 2)
 % <tr><td>Xup<td>(1 x C) upper bound on class throughput (eqs. 2-3)
 % <tr><td>Wlo<td>(K x C) per-visit residence time at station k for class c
 %               consistent with the lower throughput bound
 % </table>
%}
function [Xlo,Xup,Wlo] = pfqn_mwrbb(V,S,N,Z,sched,prio)
[K,C] = size(V);
N = N(:).';
if nargin < 4 || isempty(Z)
    Z = zeros(1,C);
end
Z = Z(:).';
if nargin < 5 || isempty(sched)
    sched = zeros(K,1);
end
sched = sched(:);
if nargin < 6 || isempty(prio)
    prio = zeros(1,C);
end
prio = prio(:).';

% no-contention upper bound on the cycle rate f_c = X_c/N_c (eqs. 1-2)
fup = zeros(1,C);
for c = 1:C
    fup(c) = 1 / (Z(c) + sum(V(:,c).*S(:,c)));
end
flo = zeros(1,C);

maxiter = 20000;
tol = 1e-13;
for it = 1:maxiter
    fup_old = fup; flo_old = flo;

    % utilization-based narrowing of the upper bounds (eq. 3)
    for c = 1:C
        cap = fup(c);
        for k = 1:K
            other = 0;
            for m = 1:C
                if m ~= c
                    other = other + N(m)*V(k,m)*S(k,m)*flo(m);
                end
            end
            denomk = N(c)*V(k,c)*S(k,c);
            if denomk > 0
                cap = min(cap, (1 - other)/denomk);
            end
        end
        fup(c) = min(fup(c), max(cap,0));
    end

    % lower-bound narrowing (Theorem 2, eq. 15). The higher-priority delay at
    % a priority station carries a 1/f_c factor and is isolated algebraically
    % as Bh so that f_c = (1 - Bh) / DEN.
    for c = 1:C
        [DEN,Bh] = mwrbb_denom(c,V,S,N,Z,fup,flo,sched,prio);
        val = (1 - Bh)/DEN;
        if val < 0, val = 0; end
        flo(c) = max(flo(c), val);
    end

    if max(abs(fup-fup_old)) < tol && max(abs(flo-flo_old)) < tol
        break;
    end
end

Xlo = N .* flo;
Xup = N .* fup;

% per-visit residence at converged rates, consistent with the lower bound
Wlo = zeros(K,C);
for c = 1:C
    for k = 1:K
        if V(k,c) == 0
            continue;
        end
        Wlo(k,c) = mwrbb_residence(k,c,V,S,N,fup,flo,sched,prio);
    end
end
end

% ------------------------------------------------------------------------
% Lower-bound denominator DEN and isolated higher-priority work Bh for class c.
function [DEN,Bh] = mwrbb_denom(c,V,S,N,Z,fup,flo,sched,prio)
[K,C] = size(V);
DEN = Z(c);
Bh = 0;
fc = flo(c);
for k = 1:K
    Vkc = V(k,c);
    if Vkc == 0
        continue;
    end
    d = sched(k);
    if d == 3 || d == 2   % priority stations contribute to Bh
        for m = 1:C
            if prio(m) < prio(c)
                Bh = Bh + N(m)*fup(m)*V(k,m)*S(k,m);
            end
        end
    end
    W = mwrbb_station_wrest(k,c,V,S,N,fup,fc,sched,prio);
    DEN = DEN + Vkc*W;
end
end

% ------------------------------------------------------------------------
% Per-visit residence at station k for class c EXCLUDING the isolated higher-
% priority (1/f_c) term. Includes own service and all bounded delay terms.
function W = mwrbb_station_wrest(k,c,V,S,N,fup,fc,sched,prio)
C = size(V,2);
Vkc = V(k,c);
Skc = S(k,c);
d = sched(k);
if d == 0            % FIFO (Theorem 1 / Lemma 1)
    s = 0;
    for m = 1:C
        if fc*Vkc == 0, pcm = 1; else pcm = min(1,(fup(m)*V(k,m))/(fc*Vkc)); end
        s = s + N(m)*S(k,m)*pcm;
    end
    W = s;           % own service is the m=c term (= N_c S_kc)
elseif d == 1        % processor sharing (Lemma 2)
    dp = 0;
    for m = 1:C
        Ncont = N(m); if m == c, Ncont = N(c)-1; end
        if fc*Vkc == 0, term = Skc; else term = min(Skc,(fup(m)*V(k,m)*S(k,m))/(fc*Vkc)); end
        dp = dp + Ncont*term;
    end
    W = Skc + dp;
elseif d == 4        % ABA full-contention (discipline-independent, P_cm = 1)
    s = 0;
    for m = 1:C
        s = s + N(m)*S(k,m);   % wait behind full service of all customers
    end
    W = s;
else                 % preemptive (3) or non-preemptive (2) priority
    dp = 0;
    for m = 1:C
        if prio(m) == prio(c)          % equal priority (includes c)
            Ncont = N(m); if m == c, Ncont = N(c)-1; end
            if fc*Vkc == 0, pcm = 1; else pcm = min(1,(fup(m)*V(k,m))/(fc*Vkc)); end
            dp = dp + Ncont*S(k,m)*pcm;
        end
        % higher-priority (prio<prio(c)) handled via Bh in mwrbb_denom
    end
    if d == 2                          % non-preemptive: lower-priority water-filling
        L = find(prio > prio(c));
        if ~isempty(L)
            [~,ord] = sort(S(k,L),'descend'); Ls = L(ord);
            budget = 1;
            for ii = 1:numel(Ls)
                l = Ls(ii);
                if fc*Vkc == 0, capr = Inf; else capr = (fup(l)*V(k,l))/(fc*Vkc); end
                if N(l) <= 0, al = 0; else al = min(budget/N(l), capr); end
                if al < 0, al = 0; end
                dp = dp + N(l)*al*S(k,l);
                budget = budget - N(l)*al;
                if budget < 0, budget = 0; end
            end
        end
    end
    W = Skc + dp;
end
end

% ------------------------------------------------------------------------
% Full per-visit residence (including higher-priority delay) for reporting Q.
function W = mwrbb_residence(k,c,V,S,N,fup,flo,sched,prio)
C = size(V,2);
fc = flo(c);
Vkc = V(k,c);
W = mwrbb_station_wrest(k,c,V,S,N,fup,fc,sched,prio);
d = sched(k);
if (d == 2 || d == 3) && fc*Vkc > 0
    for m = 1:C
        if prio(m) < prio(c)
            W = W + N(m)*fup(m)*V(k,m)*S(k,m)/(fc*Vkc);
        end
    end
end
end
