%% This example is temporarily disabled
clear solver AvgTable
fprintf(1,'This example illustrates the solution of a moderately large LQN.\n')

cwd = fileparts(which(mfilename));
model = LayeredNetwork.parseXML([cwd,filesep,'lqn_ofbiz.xml']);

%options = LQNS.defaultOptions;
%options.keep = true; % uncomment to keep the intermediate XML files generates while translating the model to LQNS

%% Solve with LQNS (disabled)
%solver{1} = LQNS(model);
%AvgTable{1} = solver{1}.getAvgTable;
%fprintf(1, '\nLQNS Results:\n');
%disp(AvgTable{1});

%% Solve with LN without initialization
solver{2} = LN(model, @(x) NC(x,'verbose',false));
Tnoinit = tic;
AvgTable{2} = solver{2}.getAvgTable;
fprintf(1, '\nLN(NC) Results:\n');
disp(AvgTable{2});
Tnoinit = toc(Tnoinit)
