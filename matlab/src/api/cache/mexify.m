% @brief MATLAB Coder script to generate MEX functions for cache_ module.
%
% This script generates MEX (MATLAB Executable) versions of cache_ functions
% for improved performance. Cache functions compute hit/miss probabilities
% and normalizing constants for cache replacement models.
%
% Skipped functions (Coder-incompatible):
%   cache_mva        - State.cartesian() OOP dependency
%   cache_t_hlru     - nested bisection helper (function handles)
%   cache_ttl_lrua   - fsolve() + rng() + function handles
%   cache_rrm_meanfield - ode23s() script (not a function)
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
mat_type = coder.typeof(0,[Inf Inf],[1 1]);

%% ===== Group: Normalizing constant computations =====

% cache_erec(gamma, m) -> E
codegen -config cfg cache_erec -args {vec_type, scal_type}

% cache_is(gamma, m, samples) -> [E, lE]
codegen -config cfg cache_is -args {vec_type, scal_type, scal_type}

% cache_spm(gamma, m, xi0) -> [Z, lZ, xi]
codegen -config cfg cache_spm -args {vec_type, scal_type, vec_type}

%% ===== Group: Hit probability computations =====

% cache_prob_erec - Skipped (cache_erec compile-time recursion exceeds codegen limit)
% codegen -config cfg cache_prob_erec -args {vec_type, scal_type}

% cache_prob_fpi(gamma, m) -> prob
codegen -config cfg cache_prob_fpi -args {vec_type, scal_type}

% cache_prob_spm - Skipped (calls cache_spm with reduced gamma, triggers cache_erec recursion limit)
% codegen -config cfg cache_prob_spm -args {vec_type, scal_type, vec_type}

% cache_prob_is - Skipped (calls cache_is which uses randperm, not codegen-compatible)
% codegen -config cfg cache_prob_is -args {vec_type, scal_type, scal_type}

%% ===== Group: Miss rate computations =====

% cache_miss - Skipped (calls cache_erec directly, recursion limit)
% codegen -config cfg cache_miss -args {vec_type, scal_type, vec_type}

% cache_miss_fpi - Skipped (calls cache_miss -> cache_erec)
% codegen -config cfg cache_miss_fpi -args {vec_type, scal_type, vec_type}

% cache_miss_spm - Skipped (calls cache_miss -> cache_erec)
% codegen -config cfg cache_miss_spm -args {vec_type, scal_type, vec_type}

% cache_miss_is - Skipped (calls cache_miss -> cache_erec)
% codegen -config cfg cache_miss_is -args {vec_type, scal_type, vec_type, scal_type}

%% ===== Group: MVA miss computation =====

% cache_mva_miss(p, m, R) -> [M, Mk]
codegen -config cfg cache_mva_miss -args {vec_type, scal_type, mat_type}

%% ===== Group: Parameter computation =====

% cache_gamma_lp(lambda, R) -> [gamma, u, n, h]
codegen -config cfg cache_gamma_lp -args {vec_type, mat_type}

%% ===== Group: Fixed-point iterations =====

% cache_xi_fp(gamma, m, xi) -> [xi, pi0, pij, it]
codegen -config cfg cache_xi_fp -args {vec_type, scal_type, vec_type}

% cache_xi_iter(gamma, m, tmax) -> z
codegen -config cfg cache_xi_iter -args {vec_type, scal_type, scal_type}

%% ===== Group: TTL approximations =====

% cache_ttl_hlru - Skipped (nested bisection helper, not codegen-compatible)
% codegen -config cfg cache_ttl_hlru -args {vec_type, scal_type}

%% ===== Group: ODE helper =====

% cache_rrm_meanfield_ode(t, x, lambda, m, n, h) -> dxdt
codegen -config cfg cache_rrm_meanfield_ode -args {scal_type, vec_type, vec_type, scal_type, scal_type, scal_type}
