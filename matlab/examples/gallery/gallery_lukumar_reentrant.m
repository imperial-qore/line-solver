function model = gallery_lukumar_reentrant(schedStrategy)
% GALLERY_LUKUMAR_REENTRANT Lu-Kumar (Rybko-Stolyar) unstable network example
%
% This example implements the classic Lu-Kumar/Rybko-Stolyar network,
% which demonstrates that queueing networks can be unstable even when
% the utilization at each station is strictly less than 1.
%
% Network structure (two chains sharing two stations):
%   - 2 stations (servers)
%   - 4 job classes in 2 chains
%   - Chain 1: Source -> Class1@Station1 -> Class2@Station2 -> Sink
%   - Chain 2: Source -> Class3@Station2 -> Class4@Station1 -> Sink
%
% Under FCFS scheduling: The network is stable when all station
% utilizations are < 1.
%
% Under the "bad" priority policy (HOL with class 1 > class 4 at
% station 1, and class 3 > class 2 at station 2): The network becomes
% unstable due to the "virtual bottleneck" phenomenon, even with
% utilization < 1 at each station.
%
% Reference:
%   Lu, S.H. and Kumar, P.R. (1991), "Distributed Scheduling Based on Due
%   Dates and Buffer Priorities", IEEE Transactions on Automatic Control,
%   Vol. 36, No. 12, pp. 1406-1416.
%
%   Rybko, A.N. and Stolyar, A.L. (1992), "Ergodicity of stochastic
%   processes describing the operation of open queueing networks",
%   Problems of Information Transmission, Vol. 28, pp. 199-220.
%
% Usage:
%   model = gallery_lukumar_reentrant();          % FCFS scheduling (stable)
%   model = gallery_lukumar_reentrant('FCFS');    % FCFS scheduling (stable)
%   model = gallery_lukumar_reentrant('HOL');     % HOL priority (unstable)
%   model = gallery_lukumar_reentrant('PS');      % Processor sharing (stable)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 1
    schedStrategy = 'FCFS';
end

% Convert string to scheduling strategy
switch upper(schedStrategy)
    case 'FCFS'
        sched1 = SchedStrategy.FCFS;
        sched2 = SchedStrategy.FCFS;
    case 'HOL'
        sched1 = SchedStrategy.HOL;
        sched2 = SchedStrategy.HOL;
    case 'PS'
        sched1 = SchedStrategy.PS;
        sched2 = SchedStrategy.PS;
    otherwise
        sched1 = SchedStrategy.FCFS;
        sched2 = SchedStrategy.FCFS;
end

model = Network('Lu-Kumar-Reentrant');

%% Block 1: Nodes
source = Source(model, 'Source');
station1 = Queue(model, 'Station1', sched1);
station2 = Queue(model, 'Station2', sched2);
sink = Sink(model, 'Sink');

%% Block 2: Job Classes
% For HOL scheduling, priority is specified as the third parameter
% Higher value = higher priority
%
% "Bad" priority policy for instability demonstration (FBFS):
%   Station 1: Class 1 (priority 1) > Class 4 (priority 0)
%   Station 2: Class 3 (priority 1) > Class 2 (priority 0)
%
% This prioritizes first-visit jobs at each station, creating a
% "virtual bottleneck" that leads to instability.

class1 = OpenClass(model, 'Class1', 1);  % Higher priority at Station 1 (Chain 1, first visit)
class2 = OpenClass(model, 'Class2', 0);  % Lower priority at Station 2 (Chain 1, second visit)
class3 = OpenClass(model, 'Class3', 1);  % Higher priority at Station 2 (Chain 2, first visit)
class4 = OpenClass(model, 'Class4', 0);  % Lower priority at Station 1 (Chain 2, second visit)

%% Block 3: Arrival Process
% External arrivals to Class 1 (Chain 1) and Class 3 (Chain 2)
arrivalRate = 0.08;  % arrival rate lambda for each chain
source.setArrival(class1, Exp(arrivalRate));  % Chain 1 arrivals
source.setArrival(class2, Disabled());
source.setArrival(class3, Exp(arrivalRate));  % Chain 2 arrivals
source.setArrival(class4, Disabled());

%% Block 4: Service Times
% ASYMMETRIC service times - classic Kumar-Seidman instability configuration
% First-visit classes (Class1, Class3) have SLOW service (m=10)
% Second-visit classes (Class2, Class4) have FAST service (m=1)
%
% This asymmetry creates the "virtual bottleneck" that causes instability
% under FBFS priority policy, even with physical utilization < 100%.
%
% Mean service times:
%   m1 = 10.0 (Class 1 at Station 1) - slow first visit
%   m2 = 1.0  (Class 2 at Station 2) - fast second visit
%   m3 = 10.0 (Class 3 at Station 2) - slow first visit
%   m4 = 1.0  (Class 4 at Station 1) - fast second visit
%
% Station utilizations (with lambda=0.08 per chain):
%   rho1 = lambda * (m1 + m4) = 0.08 * (10 + 1) = 0.88 < 1
%   rho2 = lambda * (m2 + m3) = 0.08 * (1 + 10) = 0.88 < 1

m1 = 10.0;  % Mean service time for Class 1 at Station 1 (slow)
m2 = 1.0;   % Mean service time for Class 2 at Station 2 (fast)
m3 = 10.0;  % Mean service time for Class 3 at Station 2 (slow)
m4 = 1.0;   % Mean service time for Class 4 at Station 1 (fast)

% Station 1 services (Class 1 and Class 4)
station1.setService(class1, Exp(1/m1));
station1.setService(class2, Disabled());
station1.setService(class3, Disabled());
station1.setService(class4, Exp(1/m4));

% Station 2 services (Class 2 and Class 3)
station2.setService(class1, Disabled());
station2.setService(class2, Exp(1/m2));
station2.setService(class3, Exp(1/m3));
station2.setService(class4, Disabled());

%% Block 5: Routing (Two-Chain Topology)
% Chain 1: Source -> Class1@Station1 -> Class2@Station2 -> Sink
% Chain 2: Source -> Class3@Station2 -> Class4@Station1 -> Sink
%
% This is the classic Lu-Kumar/Rybko-Stolyar network topology where
% two chains share two stations, visiting them in opposite order.

P = model.initRoutingMatrix();

% Chain 1 routing
P{class1, class1}(source, station1) = 1;     % Source -> Class1@Station1
P{class1, class2}(station1, station2) = 1;   % Class1@Station1 -> Class2@Station2
P{class2, class2}(station2, sink) = 1;       % Class2@Station2 -> Sink

% Chain 2 routing
P{class3, class3}(source, station2) = 1;     % Source -> Class3@Station2
P{class3, class4}(station2, station1) = 1;   % Class3@Station2 -> Class4@Station1
P{class4, class4}(station1, sink) = 1;       % Class4@Station1 -> Sink

model.link(P);

end
