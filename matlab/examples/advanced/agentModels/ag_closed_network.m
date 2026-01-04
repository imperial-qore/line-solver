%% Closed Network with MAM
%
% This example demonstrates MAM with RCAT methods on a closed queueing network
% with processor-sharing (PS) queues using the INAP algorithm.
%
% The RCAT (Reversed Compound Agent Theorem) decomposes the network
% into interacting stochastic processes and uses fixed-point iteration
% to compute equilibrium measures.
%
% Copyright (c) 2012-2025, Imperial College London
% All rights reserved.

clear; clc;

%% Parameters
N = 10;       % Number of jobs
mu1 = 2.0;    % Service rate at queue 1
mu2 = 1.0;    % Service rate at queue 2

%% Create model: Queue1 <-> Queue2 (closed loop)
model = Network('Closed-2Q');

queue1 = Queue(model, 'Queue1', SchedStrategy.PS);
queue2 = Queue(model, 'Queue2', SchedStrategy.PS);

cclass = ClosedClass(model, 'Class1', N, queue1);
queue1.setService(cclass, Exp(mu1));
queue2.setService(cclass, Exp(mu2));

model.link(Network.serialRouting({queue1, queue2}));

%% Solve with MAM using INAP method (default)
solverINAP = MAM(model, 'method', 'inap');
avgTableINAP = solverINAP.getAvgTable()

%% Solve with MAM using exact method (AutoCAT)
solverExact = MAM(model, 'method', 'exact');
avgTableExact = solverExact.getAvgTable()

%% Solve with MVA for comparison
solverMVA = MVA(model);
avgTableMVA = solverMVA.getAvgTable()

%% Solve with CTMC for exact results (if state space is small)
if N <= 20
    solverCTMC = CTMC(model);
    avgTableCTMC = solverCTMC.getAvgTable()
end
