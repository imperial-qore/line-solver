
%% tut01_mm1_basics
clear
% Runs under every lang: the tutorial calls LDES(), and SolverLDES pins
% options.lang='java' in its constructor because the engine is reached as a
% subprocess (common/ldes, falling back to common/ldes.jar) rather than through
% a per-language backend. LINEDefaultLang therefore does not reach it and there
% is nothing to skip. It used to call JMT(), which does reject lang='python'.
testName = 'test_tut01_mm1_basics';
run(testName);
%% tut02_mg1_multiclass_solvers
clear
% Runs under every lang: the tutorial's simulated row is LDES(), reached as a
% subprocess, so LINEDefaultLang does not select its backend. See tut01.
testName = 'test_tut02_mg1_multiclass_solvers';
run(testName);
%% tut03_repairmen
clear
testName = 'test_tut03_repairmen';
run(testName);
%% tut04_lb_routing
clear
% Runs under every lang: both routing solves are LDES(), reached as a
% subprocess, so LINEDefaultLang does not select its backend. See tut01.
testName = 'test_tut04_lb_routing';
run(testName);
%% tut05_completes_flag
clear
testName = 'test_tut05_completes_flag';
run(testName);
%% tut06_cache_lru_zipf
clear
testName = 'test_tut06_cache_lru_zipf';
run(testName);
%% tut08_respt_cdf
clear
% Runs under every lang: the simulated CDF now comes from LDES().getCdfRespT()
% rather than JMT, and LDES is a subprocess, so the lang='python' skip this
% case used to carry no longer applies.
testName = 'test_tut08_respt_cdf';
run(testName);
%% tut09_opt_load_balancing
clear
testName = 'test_tut09_opt_load_balancing';
run(testName);
%% tut10_dep_process_analysis
clear
testName = 'test_tut10_dep_process_analysis';
run(testName);
