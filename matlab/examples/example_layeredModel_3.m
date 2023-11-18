%% This example is temporarily disabled
fprintf(1,'This example illustrates the initialization of SolverLN using the output of SolverLQNS.\n')

cwd = fileparts(which(mfilename));
model = LayeredNetwork.parseXML([cwd,filesep,'example_layeredModel_1.xml']);

options = SolverLQNS.defaultOptions;
options.keep = true; % uncomment to keep the intermediate XML files generates while translating the model to LQNS

%% Solve the model using LQNS
lqnssolver = SolverLQNS(model);
[NodeAvgTable,CallAvgTable] = lqnssolver.getRawAvgTables()

%% Solve with LN without initialization
lnsolver = SolverLN(model, @(x) SolverMVA(x));
Tnoinit = tic;
AvgTable = lnsolver.getAvgTable;
AvgTable
Tnoinit = toc(Tnoinit)

%% Solve with LN with LQNS initialization
% lnsolver = SolverLN(model, @(x) SolverMVA(x));
% Tinit = tic;
% lnsolver.initFromRawAvgTables(NodeAvgTable, CallAvgTable);
% AvgTable = lnsolver.getAvgTable;
% AvgTable
% Tinit = toc(Tinit)

%% We now obtain the CDF of response times
fluidsolver = SolverFluid(model.ensemble{3});

RD = fluidsolver.getCdfRespT
% plot 'AS2=>E2' at 'Clients'
%semilogx(RD{1,4}(:,2), RD{1,4}(:,1))
%ylabel('cumulative distribution')
%xlabel('time - t')
%title('AS2=>E2 class at Clients')