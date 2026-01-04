function solver = chooseAvgSolverAI(self)
% This function takes as input a QN model defined in LINE and returns
% a Solver object with the predicted method loaded
this_model = self.model;
sn = this_model.getStruct;
dataVector = zeros(1, 15);

% Station and scheduling information
dataVector(1) = sum(sn.sched == SchedStrategy.FCFS); % Num FCFS queues
dataVector(2) = sum(sn.sched == SchedStrategy.PS); % Num PS queues
dataVector(3) = sum(sn.sched == SchedStrategy.INF); % Num delays
dataVector(4) = sn.nnodes - sn.nstations; % Num CS nodes
dataVector(5) = sum(sn.nservers(~isinf(sn.nservers))); % Num queue servers

% Job information
dataVector(6) = sn.nchains; % Num chains
dataVector(7) = sn.nclosedjobs; % Number of jobs in the system

% Service process information
numexp = 0;
numhyperexp = 0;
numerlang = 0;

for i = 1 : this_model.getNumberOfStations
    for j = 1 : this_model.getNumberOfClasses
        switch this_model.stations{i}.serviceProcess{j}.name
            case 'Exp'
                numexp = numexp + 1;
            case 'HyperExp'
                numhyperexp = numhyperexp + 1;
            case 'Erlang'
                numerlang = numerlang + 1;
        end
    end
end

dataVector(8:10) = [numexp numhyperexp numerlang]; % Num of each distribution type
dataVector(11) = mean(sn.rates, 'all', 'omitnan'); % Avg service rate
dataVector(12) = mean(sn.scv, 'all', 'omitnan'); % Avg SCV
dataVector(13) = mean(sn.phases, 'all', 'omitnan'); % Avg phases

% Misc
dataVector(14) = sum(sn.nodetype == NodeType.Queue) == 1; % If only 1 Queue, special for nc.mmint
dataVector(15) = this_model.hasProductFormSolution; % Check if model has product form solution

%Add derived features
dataVector = [dataVector dataVector(:, 1:3) ./ sum(dataVector(:, 1:3), 2)]; % Percent FCFS, PS, Delay
dataVector = [dataVector logical(dataVector(:, 4))]; % Has CS or not
dataVector = [dataVector dataVector(:, 5) ./ sum(dataVector(:, 1:2), 2)]; % Avg svrs per Queue
dataVector = [dataVector dataVector(:, 7) ./ dataVector(:, 6)]; % Num jobs per chain
dataVector = [dataVector dataVector(:, 8:10) ./ sum(dataVector(:, 8:10), 2)]; % Percent distributions

load('classifier.mat', 'classifier', 'methodNames', 'selected');
if isa(classifier, 'cell')
    chosenMethod = predictEnsemble(classifier, dataVector(selected(1:length(dataVector))));
else
    chosenMethod = predict(classifier, dataVector(selected(1:length(dataVector))));
end

solver = LINE.load(methodNames(chosenMethod), this_model);
end

