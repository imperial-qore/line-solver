%{ @file qbd_setupdelayoff_closed.m
 %  @brief Finite-population queue with setup delay and turn-off
 %
 %  @author LINE Development Team
%}

%{
 % @brief Mean queue length and throughput of a CLOSED setup/delay-off queue
 %
 % @details
 % The closed twin of qbd_setupdelayoff. The population N is finite and Z is the
 % complementary delay, the mean time a customer spends away from this station,
 % so the arrival rate is state dependent, lambda(n) = (N-n)/Z, and the level
 % index is bounded by N. That makes the chain a LEVEL-DEPENDENT QBD over
 % finitely many levels, i.e. a finite CTMC, and it is solved exactly rather
 % than by a matrix-geometric tail.
 %
 % THE SEMANTICS ARE THE SIMULATOR'S, not the mean-value shortcut's. When the
 % queue empties the server begins a delay-off period; an arrival DURING it
 % finds the server still warm and resumes without setup (Solver_ssj's
 % cancelDelayoff), and only an arrival after the delay-off has expired pays the
 % setup. That is an M/M/1 with setup time AND close-down time. The
 % per-instance cold-start race p_cold*E[setup] + S this replaces raced the
 % delay-off against the per-instance idle time and carried NO queueing term, so
 % it described a serverless instance pool rather than a single-server vacation
 % queue and left the reported response time byte-identical across a tenfold
 % change in the setup mean.
 %
 % The phase index is overloaded by level, exactly as in the open twin: at level
 % 0 phase 1 is the OFF server and the rest are the delay-off; above level 0 the
 % phases are the setup and the last one is the busy server. Both phases are
 % taken in CANONICAL COXIAN form because an arrival to an off server has to
 % enter the setup at phase 1, which is what that form guarantees for every SCV.
 %
 % @par Syntax:
 % @code
 % [QN, XN] = qbd_setupdelayoff_closed(N, Z, mu, alpharate, alphascv, betarate, betascv)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>N<td>Population of the closed chain
 % <tr><td>Z<td>Complementary delay, the mean time a customer spends away
 % <tr><td>mu<td>Service rate of the station
 % <tr><td>alpharate<td>Rate of the setup phase
 % <tr><td>alphascv<td>Squared coefficient of variation of the setup phase
 % <tr><td>betarate<td>Rate of the delay-off phase
 % <tr><td>betascv<td>Squared coefficient of variation of the delay-off phase
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>QN<td>Mean number of jobs at the station
 % <tr><td>XN<td>Throughput of the station
 % </table>
%}
function [QN, XN] = qbd_setupdelayoff_closed(N, Z, mu, alpharate, alphascv, betarate, betascv)

N = round(N);
if N <= 0 || mu <= 0
    QN = 0; XN = 0;
    return
end
Z = max(Z, GlobalConstants.FineTol);

Ta = coxian_phase(alpharate, alphascv);
na = size(Ta,1);
ta = -sum(Ta,2);
Tb = coxian_phase(betarate, betascv);
nb = size(Tb,1);
tb = -sum(Tb,2);

% Only the REACHABLE states are enumerated. The open twin pads every level to
% na+nb phases and lets the unused ones sit at zero, which a finite chain cannot
% do: an unreachable row is an absorbing row and makes the stationary solve
% singular.
off = 1;                              % level 0, server off
base = 1 + nb;                        % level 0 delay-off occupies 2..1+nb
m = base + N*(na+1);
Q = zeros(m,m);
lam = @(n) (N-n)/Z * (n < N);

if lam(0) > 0
    Q(off, idxSetup(1,1)) = Q(off, idxSetup(1,1)) + lam(0);
end
for j = 1:nb
    for j2 = 1:nb
        if j2 ~= j
            Q(idxDoff(j), idxDoff(j2)) = Q(idxDoff(j), idxDoff(j2)) + Tb(j,j2);
        end
    end
    Q(idxDoff(j), off) = Q(idxDoff(j), off) + tb(j);
    % An arrival during the delay-off cancels it and resumes WITHOUT setup.
    if lam(0) > 0
        Q(idxDoff(j), idxBusy(1)) = Q(idxDoff(j), idxBusy(1)) + lam(0);
    end
end
for n = 1:N
    for i = 1:na
        for i2 = 1:na
            if i2 ~= i
                Q(idxSetup(n,i), idxSetup(n,i2)) = Q(idxSetup(n,i), idxSetup(n,i2)) + Ta(i,i2);
            end
        end
        Q(idxSetup(n,i), idxBusy(n)) = Q(idxSetup(n,i), idxBusy(n)) + ta(i);
        % An arrival during the setup joins the queue and the setup carries on in
        % the SAME phase: the level rises, the phase does not move.
        if n < N && lam(n) > 0
            Q(idxSetup(n,i), idxSetup(n+1,i)) = Q(idxSetup(n,i), idxSetup(n+1,i)) + lam(n);
        end
    end
    if n < N && lam(n) > 0
        Q(idxBusy(n), idxBusy(n+1)) = Q(idxBusy(n), idxBusy(n+1)) + lam(n);
    end
    % A completion that empties the queue starts the delay-off at its phase 1.
    if n-1 >= 1
        Q(idxBusy(n), idxBusy(n-1)) = Q(idxBusy(n), idxBusy(n-1)) + mu;
    else
        Q(idxBusy(n), idxDoff(1)) = Q(idxBusy(n), idxDoff(1)) + mu;
    end
end
for i = 1:m
    Q(i,i) = -sum(Q(i,:));
end

p = ctmc_solve(Q);
p = max(p(:)', 0);
tot = sum(p);
if tot <= 0
    QN = 0; XN = 0;
    return
end
p = p / tot;

QN = 0;
pbusy = 0;
for n = 1:N
    lvl = p(idxBusy(n));
    for i = 1:na
        lvl = lvl + p(idxSetup(n,i));
    end
    QN = QN + n*lvl;
    pbusy = pbusy + p(idxBusy(n));
end
XN = mu * pbusy;

    function k = idxDoff(j)
        k = 1 + j;
    end

    function k = idxSetup(n, i)
        k = base + (n-1)*(na+1) + i;
    end

    function k = idxBusy(n)
        k = base + (n-1)*(na+1) + na + 1;
    end
end

function p = coxian_phase(rate, scv)
% P = COXIAN_PHASE(RATE, SCV) - canonical PH sub-generator of a phase given its
% rate, entered at phase 1 for every SCV. Exponential phases are built from the
% rate directly, which is exact and avoids the mean-based FineTol cutoff in the
% fitters: an Immediate setup has rate GlobalConstants.Immediate, whose mean is
% exactly FineTol, and the round trip turns it into an infinite rate. Same
% helper as qbd_setupdelayoff.m.
if scv == 1.0
    proc = Exp(rate).getProcess;
else
    proc = Coxian.fitMeanAndSCV(1/rate, scv).getProcess;
end
p = proc{1};
end
