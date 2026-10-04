%% G-Network (Gelenbe Network) with Negative Customers
%
% This example demonstrates SolverAG, the agent-based solver, on a G-network with negative customers.
% Negative customers (signals) remove jobs from queues when they arrive,
% modeling job cancellations or service interrupts.
%
% Network topology:
%   - Source generates positive customers (Positive) and negative signals (Negative)
%   - Positive customers flow: Source -> Queue1 -> Queue2 -> Sink
%   - Negative signals flow:   Source -> Queue1, where they are CONSUMED,
%     each removing one job from Queue1
%
% Queue1 is where the signal acts, and the reported figures are Gelenbe's
% product form exactly:
%   rho1 = lambda_pos / (mu1 + lambda_neg) = 1.0/2.3 = 0.43478, QLen1 = 0.76923
%   Queue1 throughput = mu1*rho1 = 0.86957 -- the balance is killed, not served
%   rho2 = 0.86957/mu2 = 0.28986, QLen2 = 0.40816
% Queue2 is affected only INDIRECTLY, through Queue1's reduced departure rate;
% no signal reaches it.
%
% Reference: Gelenbe, E. (1991). "Product-form queueing networks with
%            negative and positive customers", Journal of Applied Probability
%
% Copyright (c) 2012-2025, Imperial College London
% All rights reserved.

clear; clc;

%% Parameters
lambda_pos = 1.0;   % Positive customer arrival rate
lambda_neg = 0.3;   % Negative signal arrival rate
mu1 = 2.0;          % Service rate at Queue1
mu2 = 3.0;          % Service rate at Queue2

%% Create model
model = Network('GNetwork-Example');

source = Source(model, 'Source');
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS);
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');

% Positive customer class (normal jobs)
posClass = OpenClass(model, 'Positive');
source.setArrival(posClass, Exp(lambda_pos));
queue1.setService(posClass, Exp(mu1));
queue2.setService(posClass, Exp(mu2));

% Negative signal class (removes jobs from target queue)
% Using Signal class with SignalType.NEGATIVE for automatic G-network handling
negClass = Signal(model, 'Negative', SignalType.NEGATIVE);
source.setArrival(negClass, Exp(lambda_neg));
% A signal is a TRIGGER, not a job: the agent builder reads only its Source
% arrival rate, so these two service laws are inert and are set purely to keep
% the class fully specified. What the gate does require is that the SOURCE rate
% above be exponential, because the removal is folded into the agent as a
% scalar rate.
queue1.setService(negClass, Exp(mu1));
queue2.setService(negClass, Exp(mu2));

% Set routing matrix
P = model.initRoutingMatrix();
% Positive customers: Source -> Queue1 -> Queue2 -> Sink
P{posClass, posClass}(source, queue1) = 1.0;
P{posClass, posClass}(queue1, queue2) = 1.0;
P{posClass, posClass}(queue2, sink) = 1.0;
% Negative signals: Source -> Queue1, where each one removes a job and is
% itself consumed. The onward entries below are the routing a non-consumed
% class would take; a NEGATIVE signal does not survive its target, which is why
% the results show Queue1 Negative ArvR 0.3 with Tput 0 and no Negative row at
% Queue2 at all.
P{negClass, negClass}(source, queue1) = 1.0;
P{negClass, negClass}(queue1, queue2) = 1.0;
P{negClass, negClass}(queue2, sink) = 1.0;
model.link(P);

%% Solve with AG using the INAP method (finite-state fixed point)
solverAG = AG(model, 'method', 'inap');
avgTableAG = solverAG.getAvgTable();
disp('INAP (inap):');
disp(avgTableAG);

%% Solve with the matrix-geometric INAP variant (inapinf)
% 'inapinf' solves each isolated open component on its infinite state space by
% a matrix-geometric decomposition (no state-space truncation), yielding the
% exact geometric marginal pi_n = (1-rho)*rho^n. See Marin, Rota Bulo, Balsamo,
% "A Numerical Algorithm for the Decomposition of Cooperating Structured Markov
% Processes", MASCOTS 2012.
solverInf = AG(model, 'method', 'inapinf');
avgTableInf = solverInf.getAvgTable();
disp('INAP infinite-state (inapinf):');
disp(avgTableInf);
