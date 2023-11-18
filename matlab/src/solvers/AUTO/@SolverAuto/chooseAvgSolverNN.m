function solver = chooseAvgSolverNN(self)
%CHOOSEAVGSOLVERNN Summary of this function goes here
%   Detailed explanation goes here
this_model = self.model;
sn = this_model.getStruct;
dataVector = zeros(1, 15);

% Station and scheduling information
dataVector(1) = sum(sn.schedid == SchedStrategy.ID_FCFS); % Num FCFS queues
dataVector(2) = sum(sn.schedid == SchedStrategy.ID_PS); % Num PS queues
dataVector(3) = sum(sn.schedid == SchedStrategy.ID_INF); % Num delays
dataVector(4) = sn.nnodes - sn.nstations; % Num CS nodes
dataVector(5) = sum(sn.nservers(~isinf(sn.nservers))); % Num queue servers

% Job information
dataVector(6) = sn.nchains; % Num chains
dataVector(7) = sn.nclosedjobs; % Number of jobs in the system

% Service process information
numexp = 0;
numscvGT1 = 0;
numscvLT1 = 0;

for i = 1 : this_model.getNumberOfStations
    for j = 1 : this_model.getNumberOfClasses        
        switch this_model.stations{i}.serviceProcess{j}.name
            case 'Exp'
                numexp = numexp + 1;
            case 'HyperExp'
                numscvGT1 = numscvGT1 + 1;
            case 'Erlang'
                numscvLT1 = numscvLT1 + 1;
            case {'PH','APH'}
                if this_model.stations{i}.serviceProcess{j}.getSCV>1
                    numscvGT1 = numscvGT1 + 1;                    
                elseif this_model.stations{i}.serviceProcess{j}<1
                    numscvLT1 = numscvLT1 + 1;
                else
                    numexp = numexp + 1;
                end

        end
    end
end

dataVector(8:10) = [numexp numscvGT1 numscvLT1]; % Num of each distribution type
sn.rates(sn.rates>1e3)=1e3;
dataVector(11) = mean(sn.rates, 'all', 'omitnan'); % Avg service rate
dataVector(12) = mean(sn.scv, 'all', 'omitnan'); % Avg SCV
dataVector(13) = mean(sn.phases, 'all', 'omitnan'); % Avg phases

% Misc
dataVector(14) = sum(sn.nodetype == NodeType.Queue) == 1; % If only 1 Queue, special for nc.mmint
dataVector(15) = this_model.hasProductFormSolution; % Check if model has product form solution

load('classifier.mat','methodNames');

datasav = load('corrected_trained_model.mat'); trained_model = datasav.trained_model;
%trained_model = importONNXNetwork('corrected_trained_model.onnx', 'OutputLayerType', 'classification');%, 'Classes', categorical({'1','2','3'}));
%load('dataset_20000_feature.mat');
meanvalue = [0.5001    0.4999    1.0000    0.5024   10.7505    3.1566   74.3889    2.6468    2.6805    2.6726    5.3563    6.4636    7.0142    1.0000    0.5704];%mean(allData);
stdvalue = [0.5000    0.5000         0    0.7110   12.7188    2.0883   64.9605    3.9126    3.9408    3.9389    6.8807    9.0994    8.7297         0    0.4950];%std(allData);
dataVector = (dataVector - meanvalue) ./ (stdvalue + 1e-8);  % data normalization 

%dataVector
pred = predict(trained_model, dataVector);

chosenMethod = find(pred==max(pred));

chosenMethod = chosenMethod(1);
if chosenMethod == 2
    chosenMethod = 3;
elseif chosenMethod == 3
    chosenMethod = 8;
end
%if chosenMethod ~=3
%chosenMethod
%end
solver = Solver.load(methodNames(chosenMethod), this_model);
end

