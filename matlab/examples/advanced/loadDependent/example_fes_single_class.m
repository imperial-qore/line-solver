clear node jobclass

% Single-class FES Aggregation Example
%
% This tests FES aggregation with a single job class, where it should be exact.

fprintf('=== Single-Class FES Aggregation Example ===\n\n');

%% Create original 4-station tandem network with 1 class
fprintf('Creating original 4-station network...\n');

N1 = 5; % number of jobs

model = Network('OriginalModel');

% Create stations
node{1} = Delay(model, 'ThinkTime');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
node{3} = Queue(model, 'Queue2', SchedStrategy.PS);
node{4} = Queue(model, 'Queue3', SchedStrategy.PS);

% Create job class
jobclass{1} = ClosedClass(model, 'Class1', N1, node{1}, 0);

% Set service times
node{1}.setService(jobclass{1}, Exp.fitMean(5.0));
node{2}.setService(jobclass{1}, Exp.fitMean(1.5));
node{3}.setService(jobclass{1}, Exp.fitMean(1.0));
node{4}.setService(jobclass{1}, Exp.fitMean(0.8));

% Set up tandem routing
P = model.initRoutingMatrix();
P{1,1} = [0, 1, 0, 0;
          0, 0, 1, 0;
          0, 0, 0, 1;
          1, 0, 0, 0];
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

    % Find throughput at ThinkTime station
    tputOrig = 0;
    for row = 1:height(AvgTableOriginal)
        if strcmp(string(AvgTableOriginal.Station(row)), 'ThinkTime')
            tputOrig = AvgTableOriginal.Tput(row);
            break;
        end
    end

    tputFES = 0;
    for row = 1:height(AvgTableFES)
        if strcmp(string(AvgTableFES.Station(row)), 'ThinkTime')
            tputFES = AvgTableFES.Tput(row);
            break;
        end
    end

    relError = abs(tputOrig - tputFES) / max(tputOrig, 1e-10) * 100;
    fprintf('Throughput: Original=%.4f, FES=%.4f, RelError=%.2f%%\n', ...
        tputOrig, tputFES, relError);

    if relError < 1.0
        fprintf('\nSUCCESS: FES aggregation is nearly exact (< 1%% error)\n');
    else
        fprintf('\nWARNING: FES aggregation has significant error\n');
    end

catch ME
    fprintf('Error during FES aggregation: %s\n', ME.message);
    fprintf('Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
end

fprintf('\n=== Example Complete ===\n');
