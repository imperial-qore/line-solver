%{
%{
 % @brief MATLAB Coder script to generate MEX functions for sn_ module.
 %
 % This script generates MEX (MATLAB Executable) versions of sn_ functions
 % for improved performance. It configures the code generation settings
 % and specifies the expected input types (primarily the sn structure).
 %
 % See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.
%}
%}

%% Create configuration object of class 'coder.CodeConfig'.
cfg = coder.config('mex','ecoder',false);
cfg.GenerateReport = false;
cfg.ReportPotentialDifferences = false;
cfg.GenCodeOnly = false;

%% Define common types
L_type = coder.typeof(0,[Inf Inf],[1 1]);
N_type = coder.typeof(0,[1 Inf],[0 1]);
vec_type = coder.typeof(0,[1 Inf],[0 1]);
col_vec_type = coder.typeof(0,[Inf 1],[1 0]);
mat_type = coder.typeof(0,[Inf Inf],[1 1]);

%% Define sn structure type
% We define a comprehensive sn structure to cover most use cases.
% Fields must match the dynamic structure used in LINE.
S_sn = struct();
S_sn.nstations = 0;
S_sn.nclasses = 0;
S_sn.nchains = 0;
S_sn.njobs = coder.typeof(0,[1 Inf],[0 1]);
S_sn.nnodes = 0;
S_sn.nstateful = 0;
S_sn.nservers = coder.typeof(0,[Inf 1],[1 0]);
S_sn.rates = coder.typeof(0,[Inf Inf],[1 1]);
S_sn.scv = coder.typeof(0,[Inf Inf],[1 1]);
S_sn.sched = coder.typeof(0,[Inf 1],[1 0]);
S_sn.nodetype = coder.typeof(0,[Inf 1],[1 0]);
S_sn.refstat = coder.typeof(0,[Inf 1],[1 0]);
S_sn.chains = coder.typeof(0,[Inf Inf],[1 1]); 
S_sn.visits = coder.typeof({coder.typeof(0,[Inf Inf],[1 1])}, [Inf 1], [1 0]);
S_sn.nodevisits = coder.typeof({coder.typeof(0,[Inf Inf],[1 1])}, [Inf 1], [1 0]);
S_sn.inchain = coder.typeof({coder.typeof(0,[1 Inf],[0 1])}, [Inf 1], [1 0]);
S_sn.nodeToStation = coder.typeof(0,[1 Inf],[0 1]);
S_sn.nodeToStateful = coder.typeof(0,[1 Inf],[0 1]);
S_sn.stationToStateful = coder.typeof(0,[1 Inf],[0 1]);
S_sn.stationToNode = coder.typeof(0,[1 Inf],[0 1]);
S_sn.statefulToNode = coder.typeof(0,[1 Inf],[0 1]);
S_sn.refclass = coder.typeof(0,[1 Inf],[0 1]);

% Additional fields often used
S_sn.classprio = coder.typeof(0,[1 Inf],[0 1]);
S_sn.lldscaling = coder.typeof(0,[Inf Inf],[1 1]);
S_sn.fj = coder.typeof(0,[Inf Inf],[1 1]);
S_sn.routing = coder.typeof(0,[Inf Inf],[1 1]);
S_sn.procid = coder.typeof(0,[Inf Inf],[1 1]);
S_sn.phases = coder.typeof(0,[Inf Inf],[1 1]);
S_sn.phasessz = coder.typeof(0,[Inf Inf],[1 1]);
S_sn.phaseshift = coder.typeof(0,[Inf Inf],[1 1]);
S_sn.nvars = coder.typeof(0,[Inf Inf],[1 1]);
S_sn.rt = coder.typeof(0,[Inf Inf],[1 1]); 
S_sn.rtnodes = coder.typeof(0,[Inf Inf],[1 1]); 
S_sn.state = coder.typeof({coder.typeof(0,[Inf Inf],[1 1])}, [Inf 1], [1 0]); 
S_sn.space = coder.typeof({coder.typeof(0,[Inf Inf],[1 1])}, [Inf 1], [1 0]); 

sn_type = coder.typeof(S_sn);

%% Codegen commands

% Extraction and Parameter getters
codegen -config cfg sn_get_product_form_params -args {sn_type}
codegen -config cfg sn_get_demands_chain -args {sn_type}

% Property checkers
codegen -config cfg sn_has_product_form -args {sn_type}
codegen -config cfg sn_has_load_dependence -args {sn_type}
codegen -config cfg sn_has_multi_class_heter_fcfs -args {sn_type}
codegen -config cfg sn_has_priorities -args {sn_type}
codegen -config cfg sn_has_fork_join -args {sn_type}
codegen -config cfg sn_has_sd_routing -args {sn_type}
codegen -config cfg sn_is_population_model -args {sn_type}

%% ===== Additional boolean checkers (no SchedStrategy dependency) =====

% Model type checkers
codegen -config cfg sn_is_closed_model -args {sn_type}
codegen -config cfg sn_is_open_model -args {sn_type}
codegen -config cfg sn_is_mixed_model -args {sn_type}

% Class checkers
codegen -config cfg sn_has_closed_classes -args {sn_type}
codegen -config cfg sn_has_open_classes -args {sn_type}
codegen -config cfg sn_has_mixed_classes -args {sn_type}
codegen -config cfg sn_has_single_class -args {sn_type}
codegen -config cfg sn_has_multi_class -args {sn_type}
codegen -config cfg sn_has_multiple_closed_classes -args {sn_type}
codegen -config cfg sn_has_fractional_populations -args {sn_type}

% Chain checkers
codegen -config cfg sn_has_single_chain -args {sn_type}
codegen -config cfg sn_has_multi_chain -args {sn_type}

% Topology checkers
codegen -config cfg sn_has_class_switching -args {sn_type}
codegen -config cfg sn_has_multi_server -args {sn_type}
