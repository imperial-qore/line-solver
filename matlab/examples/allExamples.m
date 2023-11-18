%% run all examples
format compact
warning on backtrace
clc
fprintf(1,'<strong>This script runs all LINE examples.</strong>\n');
fprintf(1,'The current workspace will be cleared and figures will be closed. \n');
fprintf(1,'Please press a key to continue or CTRL-C to terminate.\n');
%pause; clc
clear;
close all;

%% LINE examples
fprintf(1,'\n<strong>RUNNING: example_closedModel_*</strong>');
fprintf(1,'\n\nExample: <strong>example_closedModel_1</strong>\n');
fprintf(1,'This example shows all solvers on a basic single-class closed model.\n')
clear; example_closedModel_1; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_closedModel_2</strong>\n');
fprintf('This example shows a model with a multiclass FCFS station.\n')
clear; example_closedModel_2; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_closedModel_3</strong>\n');
fprintf('This example shows the exact solution of a product-form queueing network.\n')
fprintf(1,'In this example we also calculate performance indexes by chain.\n')
clear; example_closedModel_3; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_closedModel_4</strong>\n');
fprintf(1,'This example shows state space generation for a station.')
clear; example_closedModel_4;  fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_closedModel_5</strong>\n');
fprintf(1,'This example shows a 1-line solution of a cyclic queueing network.\n');
clear; example_closedModel_5; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_closedModel_6</strong>\n');
fprintf(1,'This example shows a model with round-robin scheduling.\n');
clear; example_closedModel_6; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_closedModel_7</strong>\n');
fprintf(1,'This example shows equivalnce of models due to product-form solution (BCMP theorem).\n');
clear; example_closedModel_7; fprintf(1,'Pausing...'); pause(3.0);

%%
fprintf(1,'\n<strong>RUNNING: example_loadDependent_*</strong>');
fprintf(1,'\n\nExample: <strong>example_loadDependent_1</strong>\n');
fprintf(1,'This example shows all solvers on a basic single-class closed model.\n')
clear; example_loadDependent_1; fprintf(1,'Pausing...'); pause(3.0);

%%
fprintf(1,'\n<strong>RUNNING: example_initState_*</strong>');
fprintf(1,'\n\nExample: <strong>example_initState_1</strong>\n');
fprintf(1,'This example shows the execution of the transient solver on a 2-class 2-node class-switching model.')
clear; example_initState_1; fprintf(1,'Pausing...'); pause(3.0); try close(handleFig); end
fprintf(1,'\n\nExample: <strong>example_initState_2</strong>\n');
fprintf(1,'This example shows the execution of the transient solver on a 2-class 2-node class-switching model.')
clear; example_initState_2; fprintf(1,'Pausing...'); pause(3.0); try close(handleFig); end

%%
fprintf(1,'\n<strong>RUNNING: example_openModel_*</strong>');
fprintf(1,'\n\nExample: <strong>example_openModel_1</strong>\n');
clear; example_openModel_1; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_openModel_2</strong>\n');
clear; example_openModel_2; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_openModel_3</strong>\n');
clear; example_openModel_3; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_openModel_4</strong>\n');
clear; example_openModel_4; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_openModel_5</strong>\n');
fprintf(1,' This model examplifies how to specify models with multiple sinks (virtual sinks).\n');
clear; example_openModel_5; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_openModel_5</strong>\n');
fprintf(1,' This model examplifies a large multiclass open model.\n');
clear; example_openModel_6; fprintf(1,'Pausing...'); pause(3.0);

%%
fprintf(1,' This model examplifies how to specify priorities.\n');
clear; example_prio_1; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,' Another model that examplifies how to specify priorities.\n');
clear; example_prio_2; fprintf(1,'Pausing...'); pause(3.0);

%%
fprintf(1,'\n<strong>RUNNING: example_mixedModel_*</strong>');
fprintf(1,'\n\nExample: <strong>example_mixedModel_1</strong>\n');
clear; example_mixedModel_1; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_mixedModel_2</strong>\n');
clear; example_mixedModel_2; fprintf(1,'Pausing...'); pause(3.0);

%%
fprintf(1,'\n<strong>RUNNING: example_stateDependentRouting_*</strong>');
fprintf(1,'This example analyzes round-robin scheduling.\n');
fprintf(1,'\n\nExample: <strong>example_sdRouting_1</strong>\n');
clear; example_stateDependentRouting_1; fprintf(1,'Pausing...'); pause(3.0);

%%
fprintf(1,'\n<strong>RUNNING: example_stateProbabilities</strong>');
fprintf(1,'\n\nExample: <strong>example_stateProbabilities_1</strong>\n');
clear; example_stateProbabilities_1; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_stateProbabilities_2</strong>\n');
clear; example_stateProbabilities_2; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_stateProbabilities_3</strong>\n');
clear; example_stateProbabilities_3; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_stateProbabilities_4</strong>\n');
clear; example_stateProbabilities_4; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_stateProbabilities_5</strong>\n');
clear; example_stateProbabilities_5; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_stateProbabilities_6</strong>\n');
clear; example_stateProbabilities_6; fprintf(1,'Pausing...'); pause(3.0);

%%
fprintf(1,'\n<strong>RUNNING: example_cdfRespT_*</strong>');
fprintf(1,'\n\nExample: <strong>example_cdfRespT_1</strong>\n');
clear; example_cdfRespT_1; fprintf(1,'Pausing...'); pause(3.0); try close(handleFig); end
fprintf(1,'\n\nExample: <strong>example_cdfRespT_2</strong>\n');
clear; example_cdfRespT_2; fprintf(1,'Pausing...'); pause(3.0); try close(handleFig); end
fprintf(1,'\n\nExample: <strong>example_cdfRespT_3</strong>\n');
clear; example_cdfRespT_3; fprintf(1,'Pausing...'); pause(3.0); try close(handleFig); end
fprintf(1,'\n\nExample: <strong>example_cdfRespT_4</strong>\n');
clear; example_cdfRespT_4; fprintf(1,'Pausing...'); pause(3.0); try close(handleFig); end
fprintf(1,'\n\nExample: <strong>example_cdfRespT_5</strong>\n');
clear; example_cdfRespT_5; fprintf(1,'Pausing...'); pause(3.0); try close(handleFig); end

%%
fprintf(1,'\n<strong>RUNNING: example_randomEnvironment_*</strong>');
fprintf(1,'\n\nExample: <strong>example_randomEnvironment_1</strong>\n');
clear; example_randomEnvironment_1; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_randomEnvironment_2</strong>\n');
clear; example_randomEnvironment_2; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_randomEnvironment_3</strong>\n');
clear; example_randomEnvironment_3; fprintf(1,'Pausing...'); pause(3.0);

%%
try % LQNS must be available on the system path
fprintf(1,'\n<strong>RUNNING: example_layeredModel_*</strong>');
fprintf(1,'\n\nExample: <strong>example_layeredModel_1</strong>\n');
    clear; example_layeredModel_1; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_layeredModel_2</strong>\n');
    clear; example_layeredModel_2; fprintf(1,'Pausing...'); pause(3.0);
catch
    warning('LQNS is not available on this computer. Skipping LQN tests.');
end

%%
fprintf(1,'\n<strong>RUNNING: example_misc_*</strong>');
fprintf(1,'\n\nExample: <strong>example_misc_1</strong>\n');
fprintf(1,'This example shows how to solve only for selected performance indexes.\n');
clear; example_misc_1; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_misc_2</strong>\n');
clear; example_misc_2; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_misc_3</strong>\n');
fprintf(1,'This example illustrates the solution of DPS models.\n')
clear; example_misc_3; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_misc_4</strong>\n');
fprintf(1,'This example shows that LINE automatically checks if a solver is feasible for a given model.\n');
fprintf(1,'If not, an empty result set is returned.\n');
clear; example_misc_4; fprintf(1,'Pausing...'); pause(3.0);

%%
fprintf(1,'\n<strong>RUNNING: example_cacheModel_*</strong>');
fprintf(1,'\n\nExample: <strong>example_cacheModel_1</strong>\n');
fprintf('This example shows a small cache model with an open arrival process.\n')
clear; example_cacheModel_1; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_cacheModel_2</strong>\n');
fprintf('This example shows a small cache model with a closed arrival process.\n')
clear; example_cacheModel_2; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_cacheModel_3</strong>\n');
fprintf('This example shows a small cache model within a layered queueing network.\n')
clear; example_cacheModel_3; fprintf(1,'Pausing...'); pause(3.0);

%%
fprintf(1,'\n<strong>RUNNING: example_svcEstimation_*</strong>');
fprintf(1,'\n\nExample: <strong>example_svcEstimation_1</strong>\n');
fprintf('This example shows service demand estimation in a single class model using the utilization-based regression (UBR) method.\n')
clear; example_svcEstimation_1; fprintf(1,'Pausing...'); pause(3.0);
fprintf('This example shows service demand estimation in a multiclass model using the ERPS method.\n')
clear; example_svcEstimation_2; fprintf(1,'Pausing...'); pause(3.0);
fprintf('This example shows service demand estimation in a multiclass model using the utilization-based regression (UBR) method.\n')
clear; example_svcEstimation_3; fprintf(1,'Pausing...'); pause(3.0);
fprintf('This example shows service demand estimation in a multiclass model using the utilization-based optimization (UBO) method.\n')
clear; example_svcEstimation_4; fprintf(1,'Pausing...'); pause(3.0);

%%
fprintf(1,'\n<strong>RUNNING: example_forkJoin_*</strong>');
fprintf(1,'\n\nExample: <strong>example_forkJoin_1</strong>\n');
fprintf(1,'This example shows the simulation of a fork-join open queueing network.\n');
clear; example_forkJoin_1; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_forkJoin_2</strong>\n');
fprintf(1,'This example shows the simulation of a multiclass fork-join open queueing network.\n');
clear; example_forkJoin_2; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_forkJoin_3</strong>\n');
fprintf(1,'This example shows the simulation of nested forks and joins.\n');
clear; example_forkJoin_3; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>example_forkJoin_4</strong>\n');
fprintf(1,'This example shows a model with a fork but without a join.\n');
clear; example_forkJoin_4; fprintf(1,'Pausing...'); pause(3.0);

%%
fprintf(1,'\nExamples completed.\n')
