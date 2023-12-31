% EXAMPLE_BPMN_1 exemplifies the use of LINE to analize BPMN models. 
%
% Copyright (c) 2012-2016, Imperial College London 
% All rights reserved.

clear; 
clc;

%% input files: bpmn + extensions
filename = fullfile(pwd,'data','BPMN','bpmn_gateways.bpmn');
extFilename = fullfile(pwd,'data', 'BPMN','bpmn_ext_gateways.xml');

verbose = 1;
%% input files parsing
%% create lqn model from bpmn
myLQN = BPMN2LQN(filename, extFilename, verbose);

%% obtain line performance model from lqn
options = SolverLQNS.defaultOptions;
options.keep = true; % uncomment to keep the intermediate XML files generates while translating the model to LQNS

solver{1} = SolverLQNS(myLQN);
AvgTable{1} = solver{1}.getAvgTable();
AvgTable{1}

useLQNSnaming = true;
AvgTable{2} = solver{1}.getAvgTable(useLQNSnaming);
AvgTable{2}


useLQNSnaming = true;
[AvgTable{3}, CallAvgTable{3}] = solver{1}.getRawAvgTables();
AvgTable{3}
CallAvgTable{3}