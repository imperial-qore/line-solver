% @brief MATLAB Coder script to generate MEX functions for mam_ module.
%
% This script generates MEX (MATLAB Executable) versions of MAM (Matrix
% Analytic Methods) functions for improved performance.
%
% Skipped functions (Coder-incompatible):
%   ldqbd, ldqbd_R, ldqbd_pi           - cell array arguments
%   qbd_rg, qbd_mapmap1, qbd_raprap1   - cell array arguments
%   qbd_bmapbmap1                       - cell array arguments
%   map_ccdf_derivative                 - cell/struct MAP argument
%   map_jointpdf_derivative             - cell/struct MAP argument
%   map_m1ps_h_recursive                - cell return values
%   map_m1ps_sojourn                    - function handle validators
%   map_m1ps_cdfrespt                   - function handles + nested functions
%   mmap_compress                       - struct/cell arguments
%
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.CodeConfig'.
cfg = coder.config('mex','ecoder',false);
cfg.GenerateReport = false;
cfg.ReportPotentialDifferences = false;
cfg.GenCodeOnly = false;

%% Define common types
scal_type = coder.typeof(0);
mat_type = coder.typeof(0,[Inf Inf],[1 1]);

%% ===== Group: QBD rate matrix computations =====

% qbd_R(B, L, F, iter_max) -> R
codegen -config cfg qbd_R -args {mat_type, mat_type, mat_type, scal_type}

% qbd_R_logred(B, L, F, iter_max) -> R
codegen -config cfg qbd_R_logred -args {mat_type, mat_type, mat_type, scal_type}

%% ===== Group: MAP computations =====

% map_compute_R(C, D, mu) -> R
codegen -config cfg map_compute_R -args {mat_type, mat_type, scal_type}

%% ===== Group: Delay-off setup =====

% qbd_setupdelayoff - Skipped (uses APH/Distribution OOP classes, QBD_CR, QBD_pi)
% codegen -config cfg qbd_setupdelayoff -args {scal_type, scal_type, scal_type, scal_type, scal_type, scal_type}

%% ===== Group: Feasibility check =====

% mmdp_isfeasible(Q, R) -> bool
codegen -config cfg mmdp_isfeasible -args {mat_type, mat_type}
