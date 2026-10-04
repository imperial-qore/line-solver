% @brief MATLAB Coder script to generate MEX functions for infer_ module.
%
% This script generates MEX (MATLAB Executable) versions of demand-
% inference functions for improved performance.
%
% Skipped functions (Coder-incompatible):
%   infer_get_qlen_arrival         - Cell-array `data` argument with ragged
%                                    per-class samples.
%   infer_gibbs                    - Cell-array `data` argument and use of
%                                    exist() for default-arg detection.
%   infer_rps                      - Calls lsqnonneg() (Optimization
%                                    Toolbox), not supported by Coder.
%   infer_mlps, infer_minps,
%   infer_minps_setup, infer_fmlps - Take a Network model OOP object and
%                                    call its method getNumberOfServers().
%   infer_quick_model              - Construct Network/Source/Queue/JobClass
%                                    OOP objects.
%   infer_fluid_ps_rt_likelihood   - Returns a function handle (ode_h) and
%                                    consumes the sn struct.
%   sn_set_service_coc             - Uses ProcessType enum, cell-of-cells
%                                    fields of the sn struct.
%
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.CodeConfig'.
cfg = coder.config('mex','ecoder',false);
cfg.GenerateReport = false;
cfg.ReportPotentialDifferences = false;
cfg.GenCodeOnly = false;

%% Define common types
scal_type    = coder.typeof(0);
vec_type     = coder.typeof(0,[1 Inf],[0 1]);
col_vec_type = coder.typeof(0,[Inf 1],[1 0]);
mat_type     = coder.typeof(0,[Inf Inf],[1 1]);

%% ===== Group: Closed-form demand estimators =====

% infer_qmle(Q, N, Z) -> D
codegen -config cfg infer_qmle -args {mat_type, vec_type, vec_type}

%% ===== Group: Queue-length reconstruction =====

% infer_compute_ql_at_arrival(at, at_jobid, rt, rt_jobid, class, R) -> ql
codegen -config cfg infer_compute_ql_at_arrival -args ...
    {col_vec_type, col_vec_type, col_vec_type, col_vec_type, col_vec_type, scal_type}
