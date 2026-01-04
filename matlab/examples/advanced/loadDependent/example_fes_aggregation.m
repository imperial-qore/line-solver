clear node jobclass

% Flow-Equivalent Server (FES) Aggregation Example
%
% This example demonstrates how to use ModelAdapter.aggregateFES to replace
% a subset of stations in a closed product-form queueing network with a
% single Flow-Equivalent Server (FES).
%
% The FES has Limited Joint Dependence (LJD) service rates where the rate
% for class-c in state (n1,...,nK) equals the throughput of class-c in an
% isolated subnetwork consisting only of the subset stations.

fprintf('=== Flow-Equivalent Server (FES) Aggregation Example ===\n\n');

%% Create original 4-station tandem network with 2 classes
fprintf('Creating original 4-station network...\n');

N1 = 3; % number of class-1 jobs
N2 = 2; % number of class-2 jobs

model = Network('OriginalModel');

% Create stations
node{1} = Delay(model, 'ThinkTime');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
node{3} = Queue(model, 'Queue2', SchedStrategy.PS);
node{4} = Queue(model, 'Queue3', SchedStrategy.PS);

% Create job classes
jobclass{1} = ClosedClass(model, 'Class1', N1, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', N2, node{1}, 0);

% Set service times
node{1}.setService(jobclass{1}, Exp.fitMean(5.0));
node{1}.setService(jobclass{2}, Exp.fitMean(4.0));

node{2}.setService(jobclass{1}, Exp.fitMean(1.5));
node{2}.setService(jobclass{2}, Exp.fitMean(2.0));

node{3}.setService(jobclass{1}, Exp.fitMean(1.0));
node{3}.setService(jobclass{2}, Exp.fitMean(1.2));

node{4}.setService(jobclass{1}, Exp.fitMean(0.8));
node{4}.setService(jobclass{2}, Exp.fitMean(1.0));

% Set up tandem routing (all jobs visit all stations in order)
P = model.initRoutingMatrix();
P{1,1} = [0, 1, 0, 0;
          0, 0, 1, 0;
          0, 0, 0, 1;
          1, 0, 0, 0];
P{2,2} = P{1,1};
model.link(P);

%% Solve original model with MVA
fprintf('\n--- Solving Original Model ---\n');
solverOriginal = SolverMVA(model);
AvgTableOriginal = solverOriginal.getAvgTable;
fprintf('Original model results:\n');
disp(AvgTableOriginal);

%% Aggregate stations 2 and 3 into a Flow-Equivalent Server
fprintf('\n--- Creating FES Model ---\n');
fprintf('Aggregating Queue1 and Queue2 into a single FES...\n');

stationSubset = {node{2}, node{3}};
options.verbose = true;
options.solver = 'mva';

try
    [fesModel, fesStation, deaggInfo] = ModelAdapter.aggregateFES(model, stationSubset, options);

    fprintf('\nFES model created successfully!\n');
    fprintf('FES station name: %s\n', fesStation.getName());
    fprintf('Number of stations in FES model: %d\n', fesModel.getNumberOfStations());

    %% Solve FES model
    fprintf('\n--- Solving FES Model ---\n');
    solverFES = SolverMVA(fesModel);
    AvgTableFES = solverFES.getAvgTable;
    fprintf('FES model results:\n');
    disp(AvgTableFES);

    %% Compare throughputs
    fprintf('\n--- Throughput Comparison ---\n');

    % Extract original throughputs (at Delay station, which is common)
    for k = 1:length(jobclass)
        className = jobclass{k}.name;

        % Find throughput in original
        for row = 1:height(AvgTableOriginal)
            if strcmp(string(AvgTableOriginal.JobClass(row)), className) && ...
               strcmp(string(AvgTableOriginal.Station(row)), 'ThinkTime')
                tputOrig = AvgTableOriginal.Tput(row);
                break;
            end
        end

        % Find throughput in FES model
        for row = 1:height(AvgTableFES)
            if strcmp(string(AvgTableFES.JobClass(row)), className) && ...
               strcmp(string(AvgTableFES.Station(row)), 'ThinkTime')
                tputFES = AvgTableFES.Tput(row);
                break;
            end
        end

        relError = abs(tputOrig - tputFES) / max(tputOrig, 1e-10) * 100;
        fprintf('%s: Original=%.4f, FES=%.4f, RelError=%.2f%%\n', ...
            className, tputOrig, tputFES, relError);
    end

    %% Examine deaggregation info
    fprintf('\n--- Deaggregation Info ---\n');
    fprintf('Subset station indices: %s\n', mat2str(deaggInfo.subsetIndices));
    fprintf('Complement station indices: %s\n', mat2str(deaggInfo.complementIndices));
    fprintf('Cutoffs used: %s\n', mat2str(deaggInfo.cutoffs));
    fprintf('FES node index: %d\n', deaggInfo.fesNodeIdx);

catch ME
    fprintf('Error during FES aggregation: %s\n', ME.message);
    fprintf('Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
end

fprintf('\n=== Example Complete ===\n');
