%% Open Tandem Queue with Phase-Type Service, solved by the agent-based solver
%
% This example demonstrates that the agent-based methods represent a phase-type
% service law exactly rather than collapsing it to its mean rate. Each
% component of the agent decomposition is a QBD whose level is the queue length
% and whose phase is the pair (arrival phase, service phase), so the first
% station -- an isolated M/PH/1, since it sees the Poisson source directly --
% comes out at the Pollaczek-Khinchine mean whatever the reversed-rate
% iteration does. Both stations here have the same mean service time and differ
% only in their variability, which is exactly what an M/M/1 reading cannot see.
%
% References: Casale and Harrison, "AutoCAT: Automated Product-Form Solution of
%             Stochastic Models", Stochastic Models 27, 2013
%             Neuts, "Matrix-Geometric Solutions in Stochastic Models", 1981
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

clear; clc;

%% Parameters
lambda = 0.5;   % Arrival rate
meanS = 1.0;    % Mean service time at both queues
scv1 = 0.5;     % Erlang-2 service at queue 1 (less variable than exponential)
scv2 = 4.0;     % HyperExp service at queue 2 (more variable)

%% Create model: Source -> Queue1 -> Queue2 -> Sink
model = Network('Tandem-MPH1');

source = Source(model, 'Source');
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS);
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');

oclass = OpenClass(model, 'Class1');
source.setArrival(oclass, Exp(lambda));
queue1.setService(oclass, Erlang.fitMeanAndSCV(meanS, scv1));
queue2.setService(oclass, HyperExp.fitMeanAndSCV(meanS, scv2));

model.link(Network.serialRouting({source, queue1, queue2, sink}));

%% The exact M/G/1 mean at Queue1, which sees the Poisson source directly
rho = lambda * meanS;
pk1 = rho + rho^2 * (1 + scv1) / (2 * (1 - rho));
fprintf('=== Open Tandem with Phase-Type Service ===\n\n');
fprintf('Queue1 is an isolated M/Er2/1, so its exact mean queue length is the\n');
fprintf('Pollaczek-Khinchine value %.6f. The M/M/1 reading would be %.6f.\n\n', ...
    pk1, rho / (1 - rho));

%% Solve with the agent-based methods
avgTableINAP = AG(model, 'method', 'inap').getAvgTable()

% 'inapinf' additionally drops the maxStates truncation, solving each open
% component on its infinite state space through Neuts' rate matrix R. That
% matters most at Queue2, whose service law has the heavier tail.
avgTableINAPINF = AG(model, 'method', 'inapinf').getAvgTable()
