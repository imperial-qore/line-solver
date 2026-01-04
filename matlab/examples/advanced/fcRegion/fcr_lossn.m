%% Loss Network with Finite Capacity Region (NC Solver)
% This example demonstrates the NC solver's ability to analyze open
% loss networks using the Erlang fixed-point approximation.
% The model has a multiclass Delay node inside an FCR with DROP policy.

clear; clc;

%% Create network
model = Network('FCR Loss Network');

%% Add nodes
source = Source(model, 'Source');
delay = Delay(model, 'Delay');
sink = Sink(model, 'Sink');

%% Add job classes
class1 = OpenClass(model, 'Class1', 0);
class2 = OpenClass(model, 'Class2', 1);

%% Set arrival and service rates
lambda1 = 0.3;  lambda2 = 0.2;  % arrival rates
mu1 = 1.0;      mu2 = 0.8;      % service rates

source.setArrival(class1, Exp(lambda1));
source.setArrival(class2, Exp(lambda2));
delay.setService(class1, Exp(mu1));
delay.setService(class2, Exp(mu2));

%% Create routing matrix
P = model.initRoutingMatrix();
P.set(class1, class1, source, delay, 1.0);
P.set(class1, class1, delay, sink, 1.0);
P.set(class2, class2, source, delay, 1.0);
P.set(class2, class2, delay, sink, 1.0);
model.link(P);

%% Add finite capacity region with constraints
% When region is full, arriving jobs are dropped (lost)
fcr = model.addRegion({delay});
fcr.setGlobalMaxJobs(5);            % Global: max 5 jobs in region
fcr.setClassMaxJobs(class1, 3);     % Class1: max 3 jobs
fcr.setClassMaxJobs(class2, 3);     % Class2: max 3 jobs
fcr.setDropRule(class1, true);      % true = drop jobs
fcr.setDropRule(class2, true);      % true = drop jobs

%% Run NC solver (uses Erlang fixed-point for loss networks)
fprintf('Running NC solver (lossn method)...\n');
solverNC = SolverNC(model);
avgTableNC = solverNC.getAvgTable();

fprintf('\nNC Results:\n');
disp(avgTableNC);

%% Run JMT for comparison
fprintf('Running JMT for comparison...\n');
solverJMT = SolverJMT(model, 'seed', 23000, 'samples', 500000);
avgTableJMT = solverJMT.getAvgTable();

fprintf('\nJMT Results:\n');
disp(avgTableJMT);

%% Compare results
fprintf('\n=== Comparison ===\n');
fprintf('%-15s %12s %12s %12s\n', 'Metric', 'NC', 'JMT', 'Rel Err');
fprintf('%s\n', repmat('-', 1, 55));

% Extract Delay node throughputs
delayRowsNC = strcmp(string(avgTableNC.Station), 'Delay');
delayRowsJMT = strcmp(string(avgTableJMT.Station), 'Delay');

tputNC = avgTableNC.Tput(delayRowsNC);
tputJMT = avgTableJMT.Tput(delayRowsJMT);

for i = 1:length(tputNC)
    relErr = abs(tputNC(i) - tputJMT(i)) / tputJMT(i) * 100;
    classNames = avgTableNC.JobClass(delayRowsNC);
    fprintf('Tput %-10s %12.4f %12.4f %10.2f%%\n', ...
        string(classNames(i)), tputNC(i), tputJMT(i), relErr);
end

fprintf('\nNote: NC uses analytical Erlang fixed-point approximation.\n');
fprintf('JMT uses discrete-event simulation.\n');
