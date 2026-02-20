%% Multiclass Closed Network with MAM
%
% This example demonstrates MAM with RCAT methods on a multiclass closed queueing
% network. The RCAT algorithm handles multiple job classes by creating
% separate processes for each (station, class) pair.
%
% Copyright (c) 2012-2025, Imperial College London
% All rights reserved.

clear; clc;

%% Parameters
N1 = 5;       % Number of class 1 jobs
N2 = 3;       % Number of class 2 jobs
mu1_q1 = 2.0; % Service rate for class 1 at queue 1
mu2_q1 = 1.5; % Service rate for class 2 at queue 1
mu1_q2 = 1.0; % Service rate for class 1 at queue 2
mu2_q2 = 0.8; % Service rate for class 2 at queue 2

%% Create model
model = Network('Multiclass-Closed');

queue1 = Queue(model, 'Queue1', SchedStrategy.PS);
queue2 = Queue(model, 'Queue2', SchedStrategy.PS);

class1 = ClosedClass(model, 'Class1', N1, queue1);
class2 = ClosedClass(model, 'Class2', N2, queue1);

queue1.setService(class1, Exp(mu1_q1));
queue1.setService(class2, Exp(mu2_q1));
queue2.setService(class1, Exp(mu1_q2));
queue2.setService(class2, Exp(mu2_q2));

% Simple routing: each class cycles through both queues
P = model.initRoutingMatrix();
P{1,1} = Network.serialRouting({queue1, queue2});
P{2,2} = Network.serialRouting({queue1, queue2});
model.link(P);

%% Solve with MAM using INAP method (default)
solverINAP = MAM(model, 'method', 'inap');
avgTableINAP = solverINAP.getAvgTable()

%% Solve with MVA for comparison
solverMVA = MVA(model);
avgTableMVA = solverMVA.getAvgTable()
