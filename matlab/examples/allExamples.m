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

%% LINE examples - Basic
fprintf(1,'\n<strong>RUNNING: basic/cacheModel examples</strong>');
fprintf(1,'\n\nExample: <strong>cache_replc_rr</strong>\n');
fprintf('This example shows a small cache model with an open arrival process.\n')
clear; cache_replc_rr; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>cache_replc_fifo</strong>\n');
fprintf('This example shows a small cache model with a closed arrival process.\n')
clear; cache_replc_fifo; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>lcq_singlehost</strong>\n');
fprintf('This example shows a small cache model within a layered queueing network.\n')
clear; lcq_singlehost; fprintf(1,'Pausing...'); pause(3.0);

%%
fprintf(1,'\n<strong>RUNNING: basic/closedQN examples</strong>');
fprintf(1,'\n\nExample: <strong>cqn_repairmen</strong>\n');
fprintf(1,'This example shows all solvers on a basic single-class closed model.\n')
clear; cqn_repairmen; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>cqn_twoclass_hyperl</strong>\n');
fprintf('This example shows a model with a multiclass FCFS station.\n')
clear; cqn_twoclass_hyperl; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>cqn_threeclass_hyperl</strong>\n');
fprintf('This example shows the exact solution of a product-form queueing network.\n')
fprintf(1,'In this example we also calculate performance indexes by chain.\n')
clear; cqn_threeclass_hyperl; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>cqn_multiserver</strong>\n');
fprintf(1,'This example shows state space generation for a station.')
clear; cqn_multiserver;  fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>cqn_oneline</strong>\n');
fprintf(1,'This example shows a 1-line solution of a cyclic queueing network.\n');
clear; cqn_oneline; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>cqn_twoclass_erl</strong>\n');
fprintf(1,'This example shows a model with round-robin scheduling.\n');
clear; cqn_twoclass_erl; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>cqn_bcmp_theorem</strong>\n');
fprintf(1,'This example shows equivalnce of models due to product-form solution (BCMP theorem).\n');
clear; cqn_bcmp_theorem; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>cqn_repairmen_multi</strong>\n');
fprintf(1,'This example shows a multiclass variant of the repairmen model.\n');
clear; cqn_repairmen_multi; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>cqn_twoqueues</strong>\n');
fprintf(1,'This example shows a closed queueing network with two queues.\n');
clear; cqn_twoqueues; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>cqn_twoqueues_multi</strong>\n');
fprintf(1,'This example shows a multiclass closed queueing network with two queues.\n');
clear; cqn_twoqueues_multi; fprintf(1,'Pausing...'); pause(3.0);

%%
fprintf(1,'\n<strong>RUNNING: basic/openQN examples</strong>');
fprintf(1,'\n\nExample: <strong>oqn_basic</strong>\n');
clear; oqn_basic; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>oqn_oneline</strong>\n');
clear; oqn_oneline; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>oqn_cs_routing</strong>\n');
clear; oqn_cs_routing; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>oqn_trace_driven</strong>\n');
clear; oqn_trace_driven; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>oqn_vsinks</strong>\n');
fprintf(1,' This model examplifies how to specify models with multiple sinks (virtual sinks).\n');
clear; oqn_vsinks; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>oqn_fourqueues</strong>\n');
fprintf(1,' This model examplifies a large multiclass open model.\n');
clear; oqn_fourqueues; fprintf(1,'Pausing...'); pause(3.0);

%%
fprintf(1,'\n<strong>RUNNING: basic/mixedQN examples</strong>');
fprintf(1,'\n\nExample: <strong>mqn_basic</strong>\n');
clear; mqn_basic; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>mqn_multiserver_ps</strong>\n');
clear; mqn_multiserver_ps; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>mqn_multiserver_fcfs</strong>\n');
clear; mqn_multiserver_fcfs; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>mqn_singleserver_fcfs</strong>\n');
clear; mqn_singleserver_fcfs; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>mqn_singleserver_ps</strong>\n');
clear; mqn_singleserver_ps; fprintf(1,'Pausing...'); pause(3.0);

%%
try % LQNS must be available on the system path
fprintf(1,'\n<strong>RUNNING: basic/layeredModel examples</strong>');
fprintf(1,'\n\nExample: <strong>lqn_serial</strong>\n');
    clear; lqn_serial; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>lqn_multi_solvers</strong>\n');
    clear; lqn_multi_solvers; fprintf(1,'Pausing...'); pause(3.0);
catch
    warning('LQNS is not available on this computer. Skipping LQN tests.');
end

%%
fprintf(1,'\n<strong>RUNNING: basic/prioModel examples</strong>');
fprintf(1,' This model examplifies how to specify priorities.\n');
clear; prio_hol_open; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,' Another model that examplifies how to specify priorities.\n');
clear; prio_hol_closed; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>prio_identical</strong>\n');
fprintf(1,' This model examplifies identical priority classes.\n');
clear; prio_identical; fprintf(1,'Pausing...'); pause(3.0);

%%
fprintf(1,'\n<strong>RUNNING: basic/forkJoin examples</strong>');
fprintf(1,'\n\nExample: <strong>fj_basic_open</strong>\n');
fprintf(1,'This example shows the simulation of a fork-join open queueing network.\n');
clear; fj_basic_open; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>fj_twoclasses_forked</strong>\n');
fprintf(1,'This example shows the simulation of a multiclass fork-join open queueing network.\n');
clear; fj_twoclasses_forked; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>fj_basic_nesting</strong>\n');
fprintf(1,'This example shows the simulation of nested forks and joins.\n');
clear; fj_basic_nesting; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>fj_nojoin</strong>\n');
fprintf(1,'This example shows a model with a fork but without a join.\n');
clear; fj_nojoin; fprintf(1,'Pausing...'); pause(3.0);


%%
fprintf(1,'\n<strong>RUNNING: advanced/loadDependent examples</strong>');
fprintf(1,'\n\nExample: <strong>ld_multiserver_fcfs</strong>\n');
fprintf(1,'This example shows all solvers on a basic single-class closed model.\n')
clear; ld_multiserver_fcfs; fprintf(1,'Pausing...'); pause(3.0);

%%
fprintf(1,'\n<strong>RUNNING: advanced/initState examples</strong>');
fprintf(1,'\n\nExample: <strong>init_state_fcfs_exp</strong>\n');
fprintf(1,'This example shows the execution of the transient solver on a 2-class 2-node class-switching model.')
clear; init_state_fcfs_exp; fprintf(1,'Pausing...'); pause(3.0); try close(handleFig); end
fprintf(1,'\n\nExample: <strong>init_state_fcfs_nonexp</strong>\n');
fprintf(1,'This example shows the execution of the transient solver on a 2-class 2-node class-switching model.')
clear; init_state_fcfs_nonexp; fprintf(1,'Pausing...'); pause(3.0); try close(handleFig); end
fprintf(1,'\n\nExample: <strong>init_state_ps</strong>\n');
fprintf(1,'This example shows the execution with PS scheduling.')
clear; init_state_ps; fprintf(1,'Pausing...'); pause(3.0); try close(handleFig); end

%%
fprintf(1,'\n<strong>RUNNING: advanced/stateDepRouting examples</strong>');
fprintf(1,'This example analyzes round-robin scheduling.\n');
fprintf(1,'\n\nExample: <strong>sdroute_closed</strong>\n');
clear; sdroute_closed; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>sdroute_open</strong>\n');
fprintf(1,'This example shows state-dependent routing in an open model.\n');
clear; sdroute_open; fprintf(1,'Pausing...'); pause(3.0);

%%
fprintf(1,'\n<strong>RUNNING: advanced/stateProbabilities examples</strong>');
fprintf(1,'\n\nExample: <strong>statepr_aggr</strong>\n');
clear; statepr_aggr; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>statepr_aggr_large</strong>\n');
clear; statepr_aggr_large; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>statepr_sys_aggr</strong>\n');
clear; statepr_sys_aggr; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>statepr_sys_aggr_large</strong>\n');
clear; statepr_sys_aggr_large; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>statepr_allprobs_ps</strong>\n');
clear; statepr_allprobs_ps; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>statepr_allprobs_fcfs</strong>\n');
clear; statepr_allprobs_fcfs; fprintf(1,'Pausing...'); pause(3.0);

%%
fprintf(1,'\n<strong>RUNNING: advanced/cdfRespT examples</strong>');
fprintf(1,'\n\nExample: <strong>cdf_respt_closed</strong>\n');
clear; cdf_respt_closed; fprintf(1,'Pausing...'); pause(3.0); try close(handleFig); end
fprintf(1,'\n\nExample: <strong>cdf_respt_closed_threeclasses</strong>\n');
clear; cdf_respt_closed_threeclasses; fprintf(1,'Pausing...'); pause(3.0); try close(handleFig); end
fprintf(1,'\n\nExample: <strong>cdf_respt_open_twoclasses</strong>\n');
clear; cdf_respt_open_twoclasses; fprintf(1,'Pausing...'); pause(3.0); try close(handleFig); end
fprintf(1,'\n\nExample: <strong>cdf_respt_distrib</strong>\n');
clear; cdf_respt_distrib; fprintf(1,'Pausing...'); pause(3.0); try close(handleFig); end
fprintf(1,'\n\nExample: <strong>cdf_respt_populations</strong>\n');
clear; cdf_respt_populations; fprintf(1,'Pausing...'); pause(3.0); try close(handleFig); end

%%
fprintf(1,'\n<strong>RUNNING: advanced/randomEnv examples</strong>');
fprintf(1,'\n\nExample: <strong>renv_twostages_repairmen</strong>\n');
clear; renv_twostages_repairmen; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>renv_fourstages_repairmen</strong>\n');
clear; renv_fourstages_repairmen; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>renv_threestages_repairmen</strong>\n');
clear; renv_threestages_repairmen; fprintf(1,'Pausing...'); pause(3.0);

%%
fprintf(1,'\n<strong>RUNNING: advanced/misc examples</strong>');
fprintf(1,'\n\nExample: <strong>cqn_scheduling_dps</strong>\n');
fprintf(1,'This example illustrates the solution of DPS models.\n')
clear; cqn_scheduling_dps; fprintf(1,'Pausing...'); pause(3.0);
fprintf(1,'\n\nExample: <strong>cqn_mmpp2_service</strong>\n');
fprintf(1,'This example shows that LINE automatically checks if a solver is feasible for a given model.\n');
fprintf(1,'If not, an empty result set is returned.\n');
clear; cqn_mmpp2_service; fprintf(1,'Pausing...'); pause(3.0);

%%
%fprintf(1,'\n<strong>RUNNING: example_svcEstimation_*</strong>');
%fprintf(1,'\n\nExample: <strong>example_svcEstimation_1</strong>\n');
%fprintf('This example shows service demand estimation in a single class model using the utilization-based regression (UBR) method.\n')
%clear; example_svcEstimation_1; fprintf(1,'Pausing...'); pause(3.0);
%fprintf('This example shows service demand estimation in a multiclass model using the ERPS method.\n')
%clear; example_svcEstimation_2; fprintf(1,'Pausing...'); pause(3.0);
%fprintf('This example shows service demand estimation in a multiclass model using the utilization-based regression (UBR) method.\n')
%clear; example_svcEstimation_3; fprintf(1,'Pausing...'); pause(3.0);
%fprintf('This example shows service demand estimation in a multiclass model using the utilization-based optimization (UBO) method.\n')
%clear; example_svcEstimation_4; fprintf(1,'Pausing...'); pause(3.0);

%%
fprintf(1,'\nExamples completed.\n')
