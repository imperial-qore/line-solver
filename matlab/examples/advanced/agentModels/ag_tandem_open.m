%% Open Tandem Queue with MAM
%
% This example demonstrates MAM with RCAT methods on an open tandem queueing network
% (M/M/1 -> M/M/1) using the INAP algorithm.
%
% Reference: Marin and Rota-Bulo', "A Mean-Field Analysis of a Class of
%            Interactive Distributed Systems", MASCOTS 2009
%
% Copyright (c) 2012-2025, Imperial College London
% All rights reserved.

clear; clc;

%% Parameters
lambda = 0.5;  % Arrival rate
mu1 = 1.0;     % Service rate at queue 1
mu2 = 1.5;     % Service rate at queue 2

%% Create model: Source -> Queue1 -> Queue2 -> Sink
model = Network('Tandem-MM1');

source = Source(model, 'Source');
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS);
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');

oclass = OpenClass(model, 'Class1');
source.setArrival(oclass, Exp(lambda));
queue1.setService(oclass, Exp(mu1));
queue2.setService(oclass, Exp(mu2));

model.link(Network.serialRouting({source, queue1, queue2, sink}));

%% Solve with MAM using INAP method (default)
solverINAP = MAM(model, 'method', 'inap');
avgTableINAP = solverINAP.getAvgTable()

%% Solve with MAM using exact method (AutoCAT)
solverExact = MAM(model, 'method', 'exact');
avgTableExact = solverExact.getAvgTable()

%% Compare with MVA
solverMVA = MVA(model);
avgTableMVA = solverMVA.getAvgTable()
