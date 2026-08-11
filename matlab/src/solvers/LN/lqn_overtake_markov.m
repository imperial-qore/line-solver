function prOt = lqn_overtake_markov(clientPhases, prVisit, xj, y_aj)
% LQN_OVERTAKE_MARKOV  Overtaking probability via the LQNS Markov phased-server chain.
%
% Layer-1 port of LQNS V6 (lqns/slice.cc setRates + prOvertakingStates,
% lqns/overtake.cc computeOvertaking) for the single-conditioning case
% (calling entry == conditioning entry). Reproduces `lqns -t overtaking`:
% e.g. 31-overtaking -> prOt(server phase 2) = 0.5.
%
% Inputs
%   clientPhases : (maxPhaseA+1) x 5 matrix, row idx = client phase p (p=0..maxPhaseA).
%                  columns [nSlices, service, y_ij, y_ik, t_k] per phase, where
%                  nSlices = 1 + sum(rendezvous calls in phase p),
%                  service = total phase host residence (slice mean = service/nSlices),
%                  y_ij    = calls to the server task, y_ik = calls to other tasks,
%                  t_k     = mean time at other tasks (rendezvousDelay/y_ik).
%                  Row 1 (p=0) is the client think slice (service = thinkTime).
%   prVisit      : client entry visit probability (1 for ref/sole entry).
%   xj           : server residenceTimeForPhase(j) for the server phase j being tested.
%   y_aj         : (maxPhaseA+1) vector, y_aj(1)=total calls, y_aj(i+1)=phase-i calls.
%
% Output
%   prOt         : overtaking probability P(new arrival finds server in phase j).
%
% Reference: Franks & Woodside, "Effectiveness of early replies in client-server
% systems", Perf. Eval. 36 (1999). See [[overtaking-markov-port]].

maxPhaseA = size(clientPhases,1) - 1;
nStates = maxPhaseA + 1;

% --- setRates for each client phase (0..maxPhaseA) ---
a = zeros(1,nStates); b = zeros(1,nStates); c = zeros(1,nStates); d = zeros(1,nStates);
for idx = 1:nStates
    p = idx - 1;
    if p == maxPhaseA, prA = prVisit; else, prA = 1.0; end
    [a(idx),b(idx),c(idx),d(idx)] = local_setRates(xj, prA, ...
        clientPhases(idx,1), clientPhases(idx,2), clientPhases(idx,3), ...
        clientPhases(idx,4), clientPhases(idx,5));
end

% --- prOvertakingStates: PrOT(i+1,r+1,1)=overtaking, (:,:,2)=next ---
PrOT = zeros(nStates, nStates, 2);
for i0 = 0:maxPhaseA
    temp = 1.0;
    for r0 = 0:maxPhaseA
        if r0 == i0, continue; end
        temp = temp * local_prodOfB(b(r0+1), d(r0+1));
    end
    product = 1.0 / local_denominator(b(i0+1), d(i0+1), temp);
    r0 = i0;
    while true
        PrOT(i0+1, r0+1, 1) = c(i0+1) * product;
        PrOT(i0+1, r0+1, 2) = a(i0+1) * product;
        if r0 == 0, r0 = maxPhaseA; else, r0 = r0 - 1; end
        product = product * local_prodOfB(b(r0+1), d(r0+1));
        if r0 == i0, break; end
    end
end

% --- computeOvertaking (entA == entC): nextProb(i)=1 for the call-carrying phase ---
nextProb = zeros(1, maxPhaseA);
if maxPhaseA >= 1, nextProb(maxPhaseA) = 1.0; end
prOt = 0.0;
for i = 1:maxPhaseA
    if clientPhases(i+1,3) == 0, continue; end   % y_ij for phase i
    temp = 1.0;                                   % single-conditioning: (y_ab/y_aj_i)*(y_cd/y_cj0)=1
    acc = 0.0;
    for r = 1:maxPhaseA
        acc = acc + temp * nextProb(r) * PrOT(i+1, r+1, 1);
    end
    prOt = prOt + acc * (y_aj(1) / y_aj(i+1));    % condition
end
end

% ---------------------------------------------------------------------------
function [a,b,c,d] = local_setRates(xj, prA, nSlices, service, y_ij, y_ik, t_k)
y_sum = y_ij + y_ik + 1.0;
if nSlices ~= 0, slice = service / nSlices; else, slice = 0.0; end
temp = xj + t_k;
if ~isfinite(temp), q0=1.0; q3=0.0;
elseif temp ~= 0,   q0=xj/temp; q3=t_k/temp;
else,               q0=0.0; q3=1.0; end
temp = xj + slice;
if ~isfinite(temp), q1=0.0; q5=1.0;
elseif temp ~= 0,   q1=xj/temp; q5=slice/temp;
else,               q1=1.0; q5=0.0; end
q2 = y_ik / y_sum;  q4 = y_ij / y_sum;  q6 = prA / y_sum;
a = q5 + q1*q2*q3;  b = q1*q6;  c = q1*q4;  d = q0*q1*q2;
end

function v = local_prodOfB(b, d), v = b / (1.0 - d); end
function v = local_denominator(b, d, product), v = 1.0 - (b*product + d); end
