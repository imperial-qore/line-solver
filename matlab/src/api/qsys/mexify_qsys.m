% @brief MATLAB Coder script to generate MEX functions for qsys_ module.
%
% This script generates MEX (MATLAB Executable) versions of qsys_ functions
% for improved performance. All functions in this module are pure numerical
% computations (scalar or vector input/output).
%
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.CodeConfig'.
cfg = coder.config('mex','ecoder',false);
cfg.GenerateReport = false;
cfg.ReportPotentialDifferences = false;
cfg.GenCodeOnly = false;

%% Define common types
scal_type = coder.typeof(0);
vec_type = coder.typeof(0,[1 Inf],[0 1]);

%% ===== Group: Simple 2-scalar functions =====

% qsys_mm1(lambda, mu) -> [W, rho]
codegen -config cfg qsys_mm1 -args {scal_type, scal_type}

% qsys_gm1(sigma, mu) -> W
codegen -config cfg qsys_gm1 -args {scal_type, scal_type}

%% ===== Group: 3-scalar functions =====

% qsys_mmk(lambda, mu, k) -> [W, rho]
codegen -config cfg qsys_mmk -args {scal_type, scal_type, scal_type}

% qsys_mg1(lambda, mu, cs) -> [W, rhohat]
codegen -config cfg qsys_mg1 -args {scal_type, scal_type, scal_type}

% qsys_mm1k_loss(lambda, mu, K) -> [lossprob, rho]
codegen -config cfg qsys_mm1k_loss -args {scal_type, scal_type, scal_type}

%% ===== Group: 4-scalar functions =====

% qsys_gg1(lambda, mu, ca2, cs2) -> [W, rhohat]
codegen -config cfg qsys_gg1 -args {scal_type, scal_type, scal_type, scal_type}

% qsys_mg1k_loss_mgs(lambda, mu, mu_scv, K) -> [lossprob, rho]
codegen -config cfg qsys_mg1k_loss_mgs -args {scal_type, scal_type, scal_type, scal_type}

% qsys_mmcc_retrial_fp(lambda, mu, c, tol, maxiter) -> [blocProb, r, niter]
codegen -config cfg qsys_mmcc_retrial_fp -args {scal_type, scal_type, scal_type, scal_type, scal_type}

%% ===== Group: GI/G/1 approximations (4 scalars) =====
ARGS_GIG1 = {scal_type, scal_type, scal_type, scal_type};

codegen -config cfg qsys_gig1_lbnd -args ARGS_GIG1
codegen -config cfg qsys_gig1_ubnd_kingman -args ARGS_GIG1
codegen -config cfg qsys_gig1_approx_allencunneen -args ARGS_GIG1
codegen -config cfg qsys_gig1_approx_klb -args ARGS_GIG1
codegen -config cfg qsys_gig1_approx_kobayashi -args ARGS_GIG1
codegen -config cfg qsys_gig1_approx_marchal -args ARGS_GIG1
codegen -config cfg qsys_gig1_approx_heyman -args ARGS_GIG1
codegen -config cfg qsys_gig1_approx_kimura -args ARGS_GIG1
codegen -config cfg qsys_gig1_approx_gelenbe -args ARGS_GIG1

%% ===== Group: GI/G/1 Myskja (6 scalars) =====

% qsys_gig1_approx_myskja(lambda, mu, ca, cs, q0, qa) -> W
codegen -config cfg qsys_gig1_approx_myskja -args {scal_type, scal_type, scal_type, scal_type, scal_type, scal_type}

% qsys_gig1_approx_myskja2(lambda, mu, ca, cs, q0, qa) -> W
codegen -config cfg qsys_gig1_approx_myskja2 -args {scal_type, scal_type, scal_type, scal_type, scal_type, scal_type}

%% ===== Group: GI/G/k approximations (5 scalars) =====
ARGS_GIGK = {scal_type, scal_type, scal_type, scal_type, scal_type};

codegen -config cfg qsys_gigk_approx -args ARGS_GIGK
codegen -config cfg qsys_gigk_approx_kingman -args ARGS_GIGK
codegen -config cfg qsys_gigk_approx_whitt -args ARGS_GIGK
codegen -config cfg qsys_gigk_approx_cosmetatos -args ARGS_GIGK

%% ===== Group: M/G/1 scheduling disciplines (vector args) =====
ARGS_MG1_SCHED = {vec_type, vec_type, vec_type};

% qsys_mg1_prio(lambda, mu, cs) -> [W, rho]
codegen -config cfg qsys_mg1_prio -args ARGS_MG1_SCHED

% qsys_mg1_fb(lambda, mu, cs) -> [W, rho_total]
codegen -config cfg qsys_mg1_fb -args ARGS_MG1_SCHED

% qsys_mg1_lrpt(lambda, mu, cs) -> [W, rho_total]
codegen -config cfg qsys_mg1_lrpt -args ARGS_MG1_SCHED

% qsys_mg1_psjf(lambda, mu, cs) -> [W, rho_total]
codegen -config cfg qsys_mg1_psjf -args ARGS_MG1_SCHED

% qsys_mg1_setf(lambda, mu, cs) -> [W, rho_total]
codegen -config cfg qsys_mg1_setf -args ARGS_MG1_SCHED

% qsys_mg1_srpt(lambda, mu, cs) -> [W, rho_total]
codegen -config cfg qsys_mg1_srpt -args ARGS_MG1_SCHED
