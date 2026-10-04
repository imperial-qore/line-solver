%{
%{
 % @file pfqn_looping.m
 % @brief Eager Looping approximate MVA.
%}
%}

function [Xlo,Xup,QN,RN,it]=pfqn_looping(L,N,Z,tol,maxiter)
%{
%{
 % @brief Eager Looping approximate MVA.
 %
 % D. L. Eager, "Bounding Algorithms for Queueing Network Models of Computer
 % Systems", Ph.D. thesis, Tech. Rept. CSRG-156, University of Toronto, 1984.
 % Looping supplies the initial pessimistic and optimistic estimates that the
 % multiple-class performance bound hierarchy (pfqn_pbh) starts from, so it
 % carries a pair of bounds rather than a single fixed point.
 %
 % A HEAP H_j is the class-j congestion that the current queue-length LOWER
 % BOUNDS have not yet accounted for; it is charged back at the pessimistic
 % inflation factor V_c = max_k D_ck or the optimistic one L_c = min_k D_ck,
 % which are the largest and smallest delays one customer can inflict. The
 % whole bracket rests on Q_jk(N - 1_c) being a LOWER bound: the residence
 % (2.20) is then a lower bound, and the two heap charges bracket the
 % congestion it has not accounted for.
 %
 % The queue-length lower bounds are seeded from Little's law at the station,
 %
 %   Q_jk(N - 1_c) = X_j(N - 1_c) R_jk(N - 1_c) >= n_j D_jk / (Z_j + U_j),
 %
 % where R_jk >= D_jk (an arrival that queues behind nobody) and U_j is any
 % UPPER bound on R_j(N - 1_c): the level-0 PBH bound B_j, or the pessimistic
 % R_j(N) of the current iterate, whichever is smaller -- response time is
 % nondecreasing in the population, so R_j(N - 1_c) <= R_j(N). Each iterate
 % lowers R^(pess), which raises the seed, which lowers R^(pess) again, so the
 % refinement is monotone and every iterate is a bound.
 %
 % It previously refined the queue lengths through the convolution identity of
 % Zahorjan (1980), Q_jk(N - 1_c) = [X_j(N - 1_c)/X_j(N)] Q_jk(N), applied with
 % one class-level ratio at every station. That ratio is not a per-station
 % under-estimate (the true per-station ratios differ, so the stations whose
 % ratio is smallest come out too high), the queue lengths stopped being lower
 % bounds, and the heaps then clamped to zero and let R^(opt) collapse onto
 % R^(pess) -- a bracket of zero width around a point that is not exact. On a
 % 2-station closed model it reported X <= 0.731 against an exact 1.111.
 %
 % The level-0 multiple-class PBH bounds on the mean response time are
 %
 %   J_j(n) = sum_k D_jk,   B_j(n) = sum_k D_jk + (sum(n) - 1) max_k D_jk,
 %
 % i.e. an arriving customer queues behind nobody, respectively behind every
 % other customer in the network at its own worst centre.
 %
 % Looping is a BOUNDING algorithm, not a point estimator: it returns the
 % bracket X^(pess) of eq. (2.23) and its optimistic counterpart
 % N_c/(Z_c + R_c^(opt)) built on eq. (2.24). It is reached through SolverBA
 % as 'looping.lower'/'looping.upper'.
 %
 % @fn pfqn_looping(L, N, Z, tol, maxiter)
 % @param L Service demand matrix (stations x classes).
 % @param N Population vector.
 % @param Z Think time vector.
 % @param tol Convergence tolerance (default: 1e-6).
 % @param maxiter Maximum number of iterations (default: 1000).
 % @return Xlo Pessimistic (lower) throughput bound, one entry per class.
 % @return Xup Optimistic (upper) throughput bound, one entry per class.
 % @return QN Mean queue lengths on the pessimistic side, eq. (2.25).
 % @return RN Residence times, eq. (2.20).
 % @return it Number of iterations performed.
%}
%}

[M,R]=size(L);
if nargin<3 || isempty(Z)
    Z = zeros(1,R);
end
Z = sum(Z,1);
if nargin<4 || isempty(tol)
    tol = 1e-6;
end
if nargin<5 || isempty(maxiter)
    maxiter = 1000;
end

Dtot = sum(L,1);              % sum_k D_jk
Vpess = max(L,[],1);          % V_c = max_k D_ck
Lopt  = min(L,[],1);          % L_c = min_k D_ck
Ntot = sum(N);

% level-0 multiple-class PBH bounds on R_j, at N and at N - 1_c
Jbnd = Dtot;                                            % J_j, any population
Bm   = Dtot + max(Ntot-2,0)*Vpess;                      % B_j(N - 1_c)

beta = eye(R);                % beta_{c,j} = 1 iff j == c

% Population of class j at the stations when one class-c job is removed, at
% its optimistic and pessimistic extremes: n_j R_j/(Z_j + R_j) with R_j at the
% level-0 lower bound J_j and upper bound B_j. These cap what the queue-length
% lower bounds can legitimately account for, and the remainder is the heap.
popt = zeros(R,R);
ppess = zeros(R,R);
for c=1:R
    for j=1:R
        nj = max(N(j)-beta(c,j),0);
        if nj > 0 && Z(j)+Jbnd(j) > 0
            popt(j,c) = Jbnd(j)/(Z(j)+Jbnd(j))*nj;
        end
        if nj > 0 && Z(j)+Bm(j) > 0
            ppess(j,c) = Bm(j)/(Z(j)+Bm(j))*nj;
        end
    end
end

% Zero is the one queue-length LOWER bound available before any iterate, and
% the heaps then carry the whole population. Seeding Qm with the even split
% N/M instead, as this once did, starts the iteration from a quantity that is
% not a bound at all: the split exceeds the true queue at every station below
% average, so the very first residence (2.20) is already not a lower bound.
Qm = zeros(M,R,R);            % Qm(k,j,c) = Q_jk(N - 1_c), a LOWER bound
Hopt = popt;                  % Hopt(j,c) = H_j^(opt)(N - 1_c)
Hpess = ppess;
Ub = Bm;                      % upper bound on R_j(N - 1_c)

QN = zeros(M,R);
RN = zeros(M,R);
XN = zeros(1,R);
Rc = zeros(1,R);
Rpess = zeros(1,R);
Ropt = zeros(1,R);
for it=1:maxiter
    QN_1 = QN;
    % (2.19)-(2.23)
    for c=1:R
        if N(c) == 0
            RN(:,c) = 0; XN(c) = 0; Rc(c) = 0; Rpess(c) = 0;
            continue
        end
        Qk = sum(Qm(:,:,c),2);                          % (2.19)
        RN(:,c) = L(:,c) .* (1 + Qk);                   % (2.20)
        Rc(c) = sum(RN(:,c));                           % (2.21)
        Rpess(c) = Rc(c) + Vpess(c)*sum(Hpess(:,c));    % (2.22)
        XN(c) = N(c)/(Z(c)+Rpess(c));                   % (2.23)
    end
    % (2.24), needs X^(pess) of every class
    for c=1:R
        if N(c) == 0
            Ropt(c) = 0;
            continue
        end
        sat = -Inf;
        for ist=1:M
            den = 1 - (sum(XN.*L(ist,:)) - XN(c)*L(ist,c));
            if den > 0
                sat = max(sat, L(ist,c)*N(c)/den - Z(c));
            end
        end
        heaped = Rc(c) + Lopt(c)*sum(Hopt(:,c));
        Ropt(c) = max([sat, heaped, Dtot(c)]);
        % an optimistic bound can never exceed the pessimistic one
        Ropt(c) = min(Ropt(c), Rpess(c));
    end
    % (2.25)-(2.28)
    for c=1:R
        QN(:,c) = XN(c)*RN(:,c);
    end
    % R_j(N - 1_c) <= R_j(N) <= R_j^(pess), so the pessimistic iterate tightens
    % the upper bound the seed divides by, and never loosens it.
    for j=1:R
        if N(j) > 0 && Rpess(j) > 0 && Rpess(j) < Ub(j)
            Ub(j) = Rpess(j);
        end
    end
    for c=1:R
        for j=1:R
            nj = N(j) - beta(c,j);
            if N(j) <= 0 || nj <= 0 || Z(j)+Ub(j) <= 0
                Qm(:,j,c) = 0;
            else
                Qm(:,j,c) = (nj/(Z(j)+Ub(j))) * L(:,j);                            % (2.26)
            end
            qsum = sum(Qm(:,j,c));
            % The clamp cannot fire: Ub >= J = Dtot, so qsum <= popt by
            % construction. It stays as the statement of the invariant the
            % heaps depend on.
            Hopt(j,c) = max(0, popt(j,c) - qsum);                                  % (2.27)
            Hpess(j,c) = max(0, ppess(j,c) - qsum);                                % (2.28)
        end
    end
    nz = N > 0;
    if isempty(find(nz,1)) || (it>1 && max(max(abs(QN(:,nz)-QN_1(:,nz)))) < tol)
        break
    end
end
Xlo = XN;
Xup = zeros(1,R);
for c=1:R
    if N(c) > 0
        Xup(c) = N(c)/(Z(c)+Ropt(c));
    end
end
end
