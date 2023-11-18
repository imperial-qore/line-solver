%% This example is temporarily disabled
clear solver AvgTable
fprintf(1,'This example illustrates the solution of a moderately large LQN.\n')

cwd = fileparts(which(mfilename));
model = LayeredNetwork.parseXML([cwd,filesep,'ofbizExample.xml']);

options = SolverLQNS.defaultOptions;
options.keep = true; % uncomment to keep the intermediate XML files generates while translating the model to LQNS

%% Solve with LN without initialization
solver{1} = SolverLQNS(model);
AvgTable{1} = solver{1}.getAvgTable;
AvgTable{1}

%% Solve with LN without initialization
solver{2} = SolverLN(model, @(x) SolverNC(x,'verbose',false));
Tnoinit = tic;
AvgTable{2} = solver{2}.getAvgTable;
AvgTable{2}
Tnoinit = toc(Tnoinit)
