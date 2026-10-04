cwd = fileparts(mfilename('fullpath'))
cd(cwd);
%%
clear; cache_replc_rr; clearvars -except solver AvgTable; save ./testsExamples/cache_replc_rr.mat
clear; cache_replc_fifo; clearvars -except solver AvgTable; save ./testsExamples/cache_replc_fifo.mat
clear; cache_replc_routing; clearvars -except solver AvgTable; save ./testsExamples/cache_replc_routing.mat
clear; cache_replc_lru; clearvars -except solver AvgTable; save ./testsExamples/cache_replc_lru.mat
clear; cache_compare_replc; clearvars -except solver AvgTable; save ./testsExamples/cache_compare_replc.mat
clear; lcq_singlehost; clearvars -except solver AvgTable; save ./testsExamples/lcq_singlehost.mat
clear; lcq_threehosts; clearvars -except solver AvgTable; save ./testsExamples/lcq_threehosts.mat
clear; cdf_respt_closed; save ./testsExamples/cdf_respt_closed.mat
clear; cdf_respt_closed_threeclasses; save ./testsExamples/cdf_respt_closed_threeclasses.mat
clear; cdf_respt_open_twoclasses; save ./testsExamples/cdf_respt_open_twoclasses.mat
clear; cdf_respt_distrib; save ./testsExamples/cdf_respt_distrib.mat
clear; cs_implicit; save ./testsExamples/cs_implicit.mat
clear; cs_multi_diamond; save ./testsExamples/cs_multi_diamond.mat
clear; cs_single_diamond; save ./testsExamples/cs_single_diamond.mat
clear; cs_transient_class; save ./testsExamples/cs_transient_class.mat
clear; cqn_repairmen; save ./testsExamples/cqn_repairmen.mat
clear; cqn_twoclass_hyperl; save ./testsExamples/cqn_twoclass_hyperl.mat
clear; cqn_threeclass_hyperl; save ./testsExamples/cqn_threeclass_hyperl.mat
clear; cqn_multiserver; save ./testsExamples/cqn_multiserver.mat
clear; cqn_oneline; save ./testsExamples/cqn_oneline.mat
clear; cqn_twoclass_erl; save ./testsExamples/cqn_twoclass_erl.mat
clear; cqn_bcmp_theorem; save ./testsExamples/cqn_bcmp_theorem.mat
clear; cqn_repairmen_multi; save ./testsExamples/cqn_repairmen_multi.mat
clear; cqn_twoqueues_multi; save ./testsExamples/cqn_twoqueues_multi.mat
clear; fj_basic_open; save ./testsExamples/fj_basic_open.mat
clear; fj_twoclasses_forked; save ./testsExamples/fj_twoclasses_forked.mat
clear; fj_basic_nesting; save ./testsExamples/fj_basic_nesting.mat
clear; fj_nojoin; save ./testsExamples/fj_nojoin.mat
clear; fj_basic_closed; save ./testsExamples/fj_basic_closed.mat
clear; fj_serialfjs_open; save ./testsExamples/fj_serialfjs_open.mat
clear; fj_cs_postfork; save ./testsExamples/fj_cs_postfork.mat % buggy
clear; fj_cs_multi_visits; save ./testsExamples/fj_cs_multi_visits.mat
clear; fj_route_overlap; save ./testsExamples/fj_route_overlap.mat
clear; fj_asymm; save ./testsExamples/fj_asymm.mat
clear; fj_delays; save ./testsExamples/fj_delays.mat
clear; fj_complex_serial; save ./testsExamples/fj_complex_serial.mat
clear; fj_threebranches; save ./testsExamples/fj_threebranches.mat
clear; fj_cs_prefork; save ./testsExamples/fj_cs_prefork.mat
clear; fj_deep_nesting; save ./testsExamples/fj_deep_nesting.mat
clear; fj_serialfjs_closed; save ./testsExamples/fj_serialfjs_closed.mat
clear; init_state_fcfs_exp; save ./testsExamples/init_state_fcfs_exp.mat
clear; init_state_fcfs_nonexp; save ./testsExamples/init_state_fcfs_nonexp.mat
clear; init_state_ps; save ./testsExamples/init_state_ps.mat
clear; lqn_serial; save ./testsExamples/lqn_serial.mat 
clear; lqn_multi_solvers; save ./testsExamples/lqn_multi_solvers.mat
clear; lqn_twotasks; save ./testsExamples/lqn_twotasks.mat
clear; lqn_bpmn; save ./testsExamples/lqn_bpmn.mat
clear; lqn_workflows; save ./testsExamples/lqn_workflows.mat
clear; lqn_setup; save ./testsExamples/lqn_setup.mat
clear; lqn_basic; save ./testsExamples/lqn_basic.mat
clear; ld_multiserver_fcfs; save ./testsExamples/ld_multiserver_fcfs.mat
clear; ld_multiserver_ps_twoclasses; save ./testsExamples/ld_multiserver_ps_twoclasses.mat
clear; ld_multiserver_ps; save ./testsExamples/ld_multiserver_ps.mat
clear; ld_class_dependence; save ./testsExamples/ld_class_dependence.mat
clear; ld_joint_dependence; save ./testsExamples/ld_joint_dependence.mat
clear; cqn_scheduling_dps; save ./testsExamples/cqn_scheduling_dps.mat
clear; cqn_mmpp2_service; save ./testsExamples/cqn_mmpp2_service.mat
clear; mqn_basic; save ./testsExamples/mqn_basic.mat
clear; mqn_multiserver_ps; save ./testsExamples/mqn_multiserver_ps.mat
clear; mqn_multiserver_fcfs; save ./testsExamples/mqn_multiserver_fcfs.mat
clear; mqn_singleserver_fcfs; save ./testsExamples/mqn_singleserver_fcfs.mat
clear; mqn_singleserver_ps; save ./testsExamples/mqn_singleserver_ps.mat
clear; oqn_basic; save ./testsExamples/oqn_basic.mat
clear; oqn_oneline; save ./testsExamples/oqn_oneline.mat
clear; oqn_cs_routing; save ./testsExamples/oqn_cs_routing.mat
clear; oqn_trace_driven; save ./testsExamples/oqn_trace_driven.mat
clear; oqn_vsinks; save ./testsExamples/oqn_vsinks.mat
clear; oqn_fourqueues; save ./testsExamples/oqn_fourqueues.mat
clear; prio_hol_open; save ./testsExamples/prio_hol_open.mat
clear; prio_hol_closed; save ./testsExamples/prio_hol_closed.mat
clear; prio_psprio; save ./testsExamples/prio_psprio.mat
clear; prio_identical; save ./testsExamples/prio_identical.mat
clear; renv_twostages_repairmen; save ./testsExamples/renv_twostages_repairmen.mat
clear; renv_fourstages_repairmen; save ./testsExamples/renv_fourstages_repairmen.mat
clear; renv_threestages_repairmen; save ./testsExamples/renv_threestages_repairmen.mat
clear; sdroute_closed; save ./testsExamples/sdroute_closed.mat
clear; sdroute_twoclasses_closed; save ./testsExamples/sdroute_twoclasses_closed.mat
clear; sdroute_open; save ./testsExamples/sdroute_open.mat
%clear; swt_basic; save ./testsExamples/swt_basic.mat % switchover times
clear; polling_exhaustive_exp; save ./testsExamples/polling_exhaustive_exp.mat
clear; polling_gated; save ./testsExamples/polling_gated.mat
clear; polling_klimited; save ./testsExamples/polling_klimited.mat
clear; polling_exhaustive_det; save ./testsExamples/polling_exhaustive_det.mat
clear; statepr_aggr; save ./testsExamples/statepr_aggr.mat
clear; statepr_aggr_large; save ./testsExamples/statepr_aggr_large.mat
clear; statepr_sys_aggr; save ./testsExamples/statepr_sys_aggr.mat
clear; statepr_sys_aggr_large; save ./testsExamples/statepr_sys_aggr_large.mat
clear; statepr_allprobs_ps; save ./testsExamples/statepr_allprobs_ps.mat
clear; statepr_allprobs_fcfs; save ./testsExamples/statepr_allprobs_fcfs.mat
clear; spn_basic_open; save ./testsExamples/spn_basic_open.mat
clear; spn_open_sevenplaces; save ./testsExamples/spn_open_sevenplaces.mat
clear; spn_twomodes; save ./testsExamples/spn_twomodes.mat
clear; spn_fourmodes; save ./testsExamples/spn_fourmodes.mat
clear; spn_inhibiting; save ./testsExamples/spn_inhibiting.mat
clear; spn_closed_fourplaces; save ./testsExamples/spn_closed_fourplaces.mat
clear; spn_closed_twoplaces; save ./testsExamples/spn_closed_twoplaces.mat
clear; spn_basic_closed; save ./testsExamples/spn_basic_closed.mat
clear; tut01_mm1_basics; save ./testsExamples/tut01_mm1_basics.mat
clear; tut02_mg1_multiclass_solvers; save ./testsExamples/tut02_mg1_multiclass_solvers.mat
clear; tut03_repairmen; save ./testsExamples/tut03_repairmen.mat
clear; tut04_lb_routing; save ./testsExamples/tut04_lb_routing.mat
clear; tut05_completes_flag; save ./testsExamples/tut05_completes_flag.mat
clear; tut06_cache_lru_zipf; save ./testsExamples/tut06_cache_lru_zipf.mat
clear; tut08_respt_cdf; save ./testsExamples/tut08_respt_cdf.mat
clear; tut09_opt_load_balancing; save ./testsExamples/tut09_opt_load_balancing.mat
clear; tut10_dep_process_analysis; save ./testsExamples/tut10_dep_process_analysis.mat
close all