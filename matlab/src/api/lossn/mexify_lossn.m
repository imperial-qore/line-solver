% @brief MATLAB Coder script to generate MEX functions for lossn_ module.
%
% This script generates MEX (MATLAB Executable) versions of loss-network
% functions for improved performance.
%
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.CodeConfig'.
cfg = coder.config('mex','ecoder',false);
cfg.GenerateReport = false;
cfg.ReportPotentialDifferences = false;
cfg.GenCodeOnly = false;

%% Define common types
vec_type = coder.typeof(0,[1 Inf],[0 1]);
col_vec_type = coder.typeof(0,[Inf 1],[1 0]);
mat_type = coder.typeof(0,[Inf Inf],[1 1]);

%% ===== Group: Erlang fixed-point approximation =====

% lossn_erlangfp(nu, A, C) -> [QLen, Loss, E, niter]
codegen -config cfg lossn_erlangfp -args {vec_type, mat_type, col_vec_type}
