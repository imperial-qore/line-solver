%{
%{
 % @file pfqn_hst.m
 % @brief Operational sensitivity of throughput to homogeneous-service-time
 %        (HST) violations, and the constrained worst case (Suri 1983).
%}
%}

function sens = pfqn_hst(L,N,Z,ist)
%{
%{
 % @brief Robustness certificate for a single-class closed product-form
 %        solution: how far the predicted throughput can move when the
 %        homogeneous-service-time assumption fails at one station.
 %
 %        The HST assumption states that the mean service time at station i does
 %        not depend on the queue length there. Suri (1983) perturbs it to
 %        S_i(n) = S_i (1 + a_n), one relative deviation a_n per queue-length
 %        level n, and shows (eq. 3.11) that to first order
 %          [(1/X0) dX0/da_n] = c_n = P(n_i >= n+1)/u_i - P(n_i >= n),
 %        with u_i = L_i X0 the station utilization and the marginals taken from
 %        the product-form solution, P(n_i >= n) = L_i^n G(N-n)/G(N). The naive
 %        certificate |dX0/X0| <= (sum_n |c_n|) d follows from |a_n| <= d alone,
 %        and by Lemma 3.1 that total equals Q_i(N) - Q_i(N-1).
 %
 %        That bound is loose because the deviations are not free: an
 %        operationally consistent perturbation must leave the observed mean
 %        service time unchanged, sum_n p_n a_n = 0 with p_n = P(n_i = n). The
 %        constrained problem (P1)
 %          max |sum_n c_n a_n|  s.t.  |a_n| <= d,  sum_n p_n a_n = 0
 %        is a one-constraint linear program, solved here exactly: its optimum
 %        sets a_n = +/-d according to whether the ratio c_n/p_n exceeds a
 %        threshold, with at most one fractional coordinate. On the paper's
 %        Figure 1 system it collapses 0.831 d to 0.102 d.
 %
 %        Reference: R. Suri, "Robustness of Queuing Network Formulas",
 %        JACM 30(3):564-594, 1983 (eq. 3.11, Lemma 3.1, problem (P1)).
 % @fn pfqn_hst(L, N, Z, ist)
 % @param L Service demand vector (M x 1) of the queueing stations.
 % @param N Population (nonnegative integer scalar).
 % @param Z Think time (scalar, default 0).
 % @param ist Station index the HST perturbation is applied to (default: the
 %        bottleneck, argmax L).
 % @return sens Struct with fields:
 %        station  the station index analysed;
 %        X        product-form throughput;
 %        U        utilization of that station;
 %        Q        mean queue length there;
 %        Pgeq     P(n_i >= k), k = 0..N;
 %        p        P(n_i = k), k = 0..N;
 %        c        sensitivity coefficients c_k, k = 1..N (eq. 3.11);
 %        total    sum_k |c_k|, the unconstrained certificate per unit d;
 %        worst    the (P1) optimum per unit d;
 %        astar    the worst-case deviation profile a_k/d, k = 1..N.
%}
%}
if nargin < 3 || isempty(Z), Z = 0; end
L = L(:);
M = numel(L);
if nargin < 4 || isempty(ist)
    [~,ist] = max(L);
end
if numel(N) > 1 || size(L,2) > 1
    line_error(mfilename,'pfqn_hst is a single-class method.');
end
if ist < 1 || ist > M
    line_error(mfilename,sprintf('Station index %d is out of range (the model has %d queueing stations).',ist,M));
end
if N < 1 || abs(N-round(N)) > 0
    line_error(mfilename,'pfqn_hst requires an integer population of at least one job.');
end
if L(ist) <= 0
    line_error(mfilename,sprintf('Station %d has zero demand, so its queue-length marginals are degenerate.',ist));
end

[~,~,lg] = pfqn_rgf(L,N,Z);          % lg(k+1) = log G(k), k = 0..N

X = exp(lg(N) - lg(N+1));
y = L(ist);
Pgeq = exp((0:N)*log(y) + lg(N+1:-1:1) - lg(N+1));
p = Pgeq - [Pgeq(2:end), 0];
U = y*X;
Q = sum(Pgeq(2:end));

% eq. (3.11): c_n = P(>= n+1)/u - P(>= n), n = 1..N
c = zeros(1,N);
for n = 1:N
    if n+1 <= N
        Pn1 = Pgeq(n+2);
    else
        Pn1 = 0;
    end
    c(n) = Pn1/U - Pgeq(n+1);
end
total = sum(abs(c));

% (P1): one equality constraint plus a box. At the optimum a_n = sign(c_n -
% lambda p_n) d, so sorting by the ratio c_n/p_n and sweeping the split point
% enumerates every candidate lambda; the constraint fixes the single fractional
% coordinate at the split.
pp = p(2:end);                          % p_n for n = 1..N
[~,ord] = sort(c./max(pp,realmin), 'descend');
best = 0; astar = zeros(1,N);
for k = 0:N
    a = -ones(1,N);
    a(ord(1:k)) = 1;
    for piv = 1:N
        aa = a;
        rest = pp*aa' - pp(piv)*aa(piv);
        if pp(piv) <= 0, continue; end
        v = -rest/pp(piv);
        if v < -1 || v > 1, continue; end
        aa(piv) = v;
        obj = c*aa';
        if abs(obj) > abs(best)
            best = obj; astar = aa;
        end
    end
end
if best < 0
    best = -best; astar = -astar;       % the feasible set is symmetric
end

sens = struct('station',ist,'X',X,'U',U,'Q',Q,'Pgeq',Pgeq,'p',p, ...
    'c',c,'total',total,'worst',best,'astar',astar);
end
