%{
%{
 % @brief MATLAB Coder script to generate MEX functions for fj_ module.
 %
 % This script generates MEX (MATLAB Executable) versions of fj_ functions
 % for improved performance. It configures the code generation settings
 % and specifies the expected input types.
 %
 % See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.
%}
%}

%% Create configuration object of class 'coder.CodeConfig'.
cfg = coder.config('mex','ecoder',false);
cfg.GenerateReport = false;
cfg.ReportPotentialDifferences = false;
cfg.GenCodeOnly = false;

%% Define types
scal_type = coder.typeof(0);
vec_type = coder.typeof(0,[1 Inf],[0 1]);
str_type = coder.typeof('a', [1 Inf], [0 1]);

%% Codegen commands

% Basic Functions
codegen -config cfg fj_harmonic -args {scal_type}
codegen -config cfg fj_respt_2way -args {scal_type, scal_type}
codegen -config cfg fj_respt_vm -args {scal_type, scal_type, scal_type}
codegen -config cfg fj_respt_nt -args {scal_type, scal_type, scal_type}
codegen -config cfg fj_respt_varki -args {scal_type, scal_type, scal_type}
codegen -config cfg fj_bounds -args {scal_type, scal_type, scal_type}
codegen -config cfg fj_rmax -args {scal_type, scal_type, scal_type}
codegen -config cfg fj_synch_delay -args {scal_type, scal_type}
codegen -config cfg fj_sm_tput -args {scal_type, scal_type}

% Approximation Functions
codegen -config cfg fj_xmax_approx -args {scal_type, scal_type, scal_type, str_type}
codegen -config cfg fj_xmax_2 -args {scal_type, scal_type}
codegen -config cfg fj_xmax_exp -args {scal_type, scal_type}
codegen -config cfg fj_xmax_erlang -args {scal_type, scal_type, scal_type}
codegen -config cfg fj_xmax_hyperexp -args {scal_type, scal_type, scal_type, scal_type}
codegen -config cfg fj_xmax_pareto -args {scal_type, scal_type, scal_type}
codegen -config cfg fj_xmax_normal -args {scal_type, scal_type, scal_type, str_type}
codegen -config cfg fj_xmax_emma -args {scal_type, scal_type, str_type} % param assumed scalar for exp

% Complex/Variant Functions
codegen -config cfg fj_rmax_evd -args {scal_type, scal_type, scal_type, scal_type}
codegen -config cfg fj_rmax_erlang -args {scal_type, scal_type, scal_type, scal_type}
% fj_gk_bound - Skipped (returns struct or scalar depending on type, mixed return types)
% codegen -config cfg fj_gk_bound -args {scal_type, str_type}
codegen -config cfg fj_char_max -args {scal_type, vec_type, str_type}
