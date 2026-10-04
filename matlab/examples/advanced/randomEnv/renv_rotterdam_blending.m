% renv_rotterdam_blending - Exact reproduction of the fyp26 SOQN blending result
%
% Reproduces the 24-environment exponential blending accuracy result of the
% fyp26 Rotterdam container-terminal study (Dhingra 2018), using the same
% level-dependent QBD blocks and the cyclic resolvent blending that the
% ENV state-vector analyzer uses internally.
%
% Per-environment model (PoissonSOQNSpec): a semi-open SOQN with N tokens, a
% Poisson arrival rate lambda_h (hour h of the day), and two flow-equivalent
% server (FES) subnetworks in tandem (upstream S1, downstream S2), each a
% load-dependent station whose throughput curve mu1(n)/mu2(n) is calibrated by
% exact MVA on the underlying closed subnetwork. The LD-QBD state is (n,k):
% level n = jobs upstream of S2 (S1 + external backlog), phase k = S2 occupancy.
%
% The 24 hourly environments are visited cyclically with exponential sojourns
% (mean = hour duration). Each visit applies the resolvent s*pi*(sI-Q)^{-1}
% (exit = time-average, memoryless), entry vectors are chained to an L1 fixed
% point, and the blend weights stages by their time fraction. Mean external
% wait W and external queue length Qex are compared to the Dhingra simulation
% (Table 13).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

clear; warning off;

%% Dhingra 24-hour schedule and simulation ground truth (Table 13)
dailyRateHr = [  6,  30,  40,  62,  76,  79, 119, 164, 152, 130,  79,  70, ...
                57,  57, 113, 130, 162, 202, 148, 118,  92,  62,  36,   8];
dailyDurHr  = [0.51,0.84,0.76,0.98,0.67,2.33,1.80,0.62,0.36,1.30,1.02,0.98, ...
               0.86,1.56,1.30,1.57,0.34,0.22,0.47,1.33,1.56,0.93,0.71,1.11];
lambda = dailyRateHr * 0.5 / 60;     % arrivals per minute  (buildLambdasMin)
durMin = dailyDurHr * 60;            % hour durations in minutes
fracW  = durMin / sum(durMin);       % time-fraction blend weights f_i

simN   = 24:34;
simW   = [154.0661,118.3936,99.962,88.2379,80.9272,75.0887,70.745,67.4326,64.9719,62.7913,61.2534];
simQex = [89.8508,63.3970,49.3650,40.8783,35.3546,30.6543,27.4111,24.7218,22.3622,20.3893,18.9252];

tailFactor = 15;                     % M_trunc = N*(1+tailFactor) (ACC_TAIL)
params = struct('numEntryServers',6,'numStacks',29,'entryServiceTime',6.0, ...
    'travelToStackTime',5.6,'stackServiceTime',6.0,'numExitServers',6, ...
    'travelToExitTime',5.6,'exitServiceTime',6.0);

%% Reproduce the blending for a few token counts and compare to simulation
Nlist = [24 28 34];
fprintf('Rotterdam SOQN exponential blending vs Dhingra simulation (tailFactor=%d)\n', tailFactor);
fprintf('%4s | %10s %10s %7s | %10s %10s %7s\n', 'N', 'W_blend','W_sim','err%','Qex_bl','Qex_sim','err%');
for N = Nlist
    mu1 = fesCurve(N, params, 'up');
    mu2 = fesCurve(N, params, 'down');
    [Wbl, Qexbl] = blendSOQN(N, lambda, durMin, fracW, mu1, mu2, tailFactor);
    j = find(simN == N);
    fprintf('%4d | %10.4f %10.4f %6.2f | %10.4f %10.4f %6.2f\n', ...
        N, Wbl, simW(j), 100*abs(Wbl-simW(j))/simW(j), ...
        Qexbl, simQex(j), 100*abs(Qexbl-simQex(j))/simQex(j));
end

%% ---- local functions ------------------------------------------------------
function mu = fesCurve(N, p, which)
% Load-dependent FES throughput curve mu(l), l=0..N, by exact MVA on the
% closed subnetwork at each population l (Norton flow-equivalent).
mu = zeros(1, N+1);
for l = 1:N
    if strcmp(which, 'up')
        m = Network('S1');
        eg = Queue(m,'EntryGates',SchedStrategy.FCFS); eg.setNumberOfServers(p.numEntryServers);
        tv = Delay(m,'TravelToStack');
        st = cell(1,p.numStacks);
        for i=1:p.numStacks, st{i}=Queue(m,sprintf('Stack%d',i),SchedStrategy.FCFS); st{i}.setNumberOfServers(1); end
        cls = ClosedClass(m,'Trucks',l,eg);
        eg.setService(cls,Exp(1/p.entryServiceTime));
        tv.setService(cls,Exp(1/p.travelToStackTime));
        for i=1:p.numStacks, st{i}.setService(cls,Exp(1/p.stackServiceTime)); end
        ns = 2 + p.numStacks; P = m.initRoutingMatrix();
        P{1}(1,2) = 1;                          % EntryGates -> TravelToStack
        P{1}(2,3:ns) = 1/p.numStacks;           % TravelToStack -> Stack_i
        P{1}(3:ns,1) = 1;                        % Stack_i -> EntryGates
        m.link(P);
        T = MVA(m,'method','exact','verbose',false).getAvgTable;
        mu(l+1) = T.Tput(1);                     % throughput at EntryGates
    else
        m = Network('S2');
        tv = Delay(m,'TravelToExit');
        xg = Queue(m,'ExitGates',SchedStrategy.FCFS); xg.setNumberOfServers(p.numExitServers);
        cls = ClosedClass(m,'Trucks',l,tv);
        tv.setService(cls,Exp(1/p.travelToExitTime));
        xg.setService(cls,Exp(1/p.exitServiceTime));
        P = m.initRoutingMatrix(); P{1}(1,2)=1; P{1}(2,1)=1; m.link(P);
        T = MVA(m,'method','exact','verbose',false).getAvgTable;
        mu(l+1) = T.Tput(2);                     % throughput at ExitGates
    end
end
end

function Q = soqnGenerator(N, lam, mu1, mu2, tailFactor)
% Flat generator of the Poisson SOQN LD-QBD. State (n,k), flat index n*M+k+1.
% level n=0..Mtr (S1+backlog), phase k=0..N (S2 occupancy).
M = N + 1; Mtr = N + tailFactor*N; dim = (Mtr+1)*M;
ri = zeros(0,1); ci = zeros(0,1); vi = zeros(0,1);
for n = 0:Mtr
    for k = 0:N
        row = n*M + k + 1;
        j   = min(n, N-k);
        m1  = mu1(j+1);            % mu1.at(min(n,N-k))
        m2  = mu2(k+1);            % mu2.at(k)
        atTop = (n == Mtr);
        lamEff = lam * (~atTop);   % no arrivals leave the truncated chain at top
        % diagonal
        ri(end+1,1)=row; ci(end+1,1)=row; vi(end+1,1)=-(lamEff + m1 + m2); %#ok<AGROW>
        % S2 completion within level: (n,k)->(n,k-1) at mu2(k)
        if k>=1
            ri(end+1,1)=row; ci(end+1,1)=n*M+(k-1)+1; vi(end+1,1)=m2; %#ok<AGROW>
        end
        % arrival up: (n,k)->(n+1,k) at lam
        if n<Mtr
            ri(end+1,1)=row; ci(end+1,1)=(n+1)*M+k+1; vi(end+1,1)=lam; %#ok<AGROW>
        end
        % S1 completion down: (n,k)->(n-1,k+1) at mu1(min(n,N-k))
        if n>=1 && k<=N-1
            ri(end+1,1)=row; ci(end+1,1)=(n-1)*M+(k+1)+1; vi(end+1,1)=m1; %#ok<AGROW>
        end
    end
end
Q = sparse(ri, ci, vi, dim, dim);
end

function [W, Qex] = blendSOQN(N, lambda, durMin, fracW, mu1, mu2, tailFactor)
% Cyclic exponential blending over the 24 hourly SOQN environments, via the
% resolvent s*pi*(sI-Q)^{-1}, to an L1 fixed point; then blend and extract W,Qex.
K = numel(lambda);
M = N + 1; Mtr = N + tailFactor*N; dim = (Mtr+1)*M;
Q = cell(1,K);
for h = 1:K, Q{h} = soqnGenerator(N, lambda(h), mu1, mu2, tailFactor); end

piEnter = cell(1,K);
pi0 = sparse(1, dim); pi0(1) = 1;        % empty system
piEnter{1} = pi0;
for h = 2:K, piEnter{h} = pi0; end

Ssp = speye(dim);
for iter = 1:200
    prev = piEnter;
    for h = 1:K
        s = 1/durMin(h);
        piEnter{mod(h,K)+1} = s * (piEnter{h} / (s*Ssp - Q{h}));   % resolvent = exit
    end
    l1 = 0;
    for h = 1:K, l1 = max(l1, sum(abs(piEnter{h}-prev{h}))); end
    if l1 < 1e-10, break; end
end

% Blend (exp sojourn: exit = time-average), weighted by time fraction
piAvg = sparse(1, dim);
for h = 1:K
    s = 1/durMin(h);
    piEx = s * (piEnter{h} / (s*Ssp - Q{h}));
    piAvg = piAvg + fracW(h) * piEx;
end
piAvg = full(piAvg); piAvg(piAvg<0)=0; piAvg = piAvg/sum(piAvg);

% Metrics (extractMetricsPoissonSOQN)
QLex=0; QL1=0; QL2=0; tput=0;
for n = 0:Mtr
    for k = 0:N
        p = piAvg(n*M+k+1);
        if p==0, continue; end
        extQ = max(0, n+k-N); s1 = min(n, N-k);
        QLex = QLex + extQ*p; QL1 = QL1 + s1*p; QL2 = QL2 + k*p;
        if k>=1, tput = tput + mu2(k+1)*p; end
    end
end
Qex = QLex;
W   = (QLex + QL1 + QL2) / tput;     % external wait + internal time (Little)
end
