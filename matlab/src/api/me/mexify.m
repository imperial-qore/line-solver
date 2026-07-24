% @brief MATLAB Coder script to generate MEX functions for me_ module.
%
% This script generates MEX (MATLAB Executable) versions of Maximum Entropy
% (ME) queueing-network functions for improved performance.
%
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.CodeConfig'.
cfg = coder.config('mex','ecoder',false);
cfg.GenerateReport = false;
cfg.ReportPotentialDifferences = false;
cfg.GenCodeOnly = false;

%% Define common types
scal_type = coder.typeof(0);
mat_type  = coder.typeof(0,[Inf Inf],[1 1]);
mat3_type = coder.typeof(0,[Inf Inf Inf],[1 1 1]);

% Options struct: must declare all fields actually accessed.
S_opt = struct();
S_opt.tol     = 0;
S_opt.maxiter = 0;
S_opt.verbose = false;
opt_type = coder.typeof(S_opt);

%% ===== Group: Maximum-Entropy Open Queueing Networks =====

% me_oqn(M, R, lambda0, Ca0, mu, Cs, P, c, insens, options)
%   -> [L, W, Ca, Cd, lambda, rho, iter]
vec_type = coder.typeof(0,[Inf 1],[1 0]);
lvec_type = coder.typeof(false,[Inf 1],[1 0]);
codegen -config cfg me_oqn -args ...
    {scal_type, scal_type, mat_type, mat_type, mat_type, mat_type, mat3_type, vec_type, lvec_type, opt_type}
