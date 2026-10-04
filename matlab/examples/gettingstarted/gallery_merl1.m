function [model,source,queue,sink,oclass] = gallery_merl1
% GALLERY_MERL1 Create a M/E_r/1 queue model
%
% This function creates a simple M/E_r/1 queueing model:
% - Markovian (exponential) arrivals
% - Erlang service distribution  
% - Single server
%
% Output:
%   model - Network model object
%   source - Source node
%   queue - Queue node with Erlang service
%   sink - Sink node
%   oclass - Open job class

% Create the network model
model = Network('M/E_r/1');

%% Block 1: nodes
source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');

%% Block 2: classes  
oclass = OpenClass(model, 'Class1');

% Set exponential arrivals with rate 1
source.setArrival(oclass, Exp(1));

% Set Erlang service with mean 0.5 and order 2 (E_2 distribution)
% This gives coefficient of variation = 1/sqrt(2) ≈ 0.707
queue.setService(oclass, Erlang.fitMeanAndOrder(0.5, 2));

%% Block 3: topology
model.link(Network.serialRouting(source,queue,sink));

end