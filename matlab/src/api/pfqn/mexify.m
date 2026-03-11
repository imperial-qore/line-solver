%{
%{
 % @brief MATLAB Coder script to generate MEX function for pfqn_bs.
 %
 % This script generates a MEX (MATLAB Executable) version of pfqn_bs
 % for improved performance. It configures the code generation settings
 % and specifies the expected input types and dimensions.
 %
 % See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.
%}
%}
% UNTITLED   Generate static library pfqn_bs
%
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.CodeConfig'.
cfg = coder.config('mex','ecoder',false);
cfg.GenerateReport = false;
cfg.ReportPotentialDifferences = false;
cfg.GenCodeOnly = false;

%% Define common types
L_type = coder.typeof(0,[Inf Inf],[1 1]);
N_type = coder.typeof(0,[1 Inf],[0 1]);
Z_type = coder.typeof(0,[1 Inf],[0 1]);
scal_type = coder.typeof(0);
vec_type = coder.typeof(0,[1 Inf],[0 1]);
col_vec_type = coder.typeof(0,[Inf 1],[1 0]);
mat_type = coder.typeof(0,[Inf Inf],[1 1]);
vec_flex_type = coder.typeof(0,[Inf Inf],[1 1]); % Flexible 1D/2D
str_type = coder.typeof('a', [1 Inf], [0 1]); % Variable string

%% Group: (L,N,Z)
ARGS_LNZ = {L_type, N_type, Z_type};

codegen -config cfg pfqn_ca -args ARGS_LNZ
codegen -config cfg pfqn_le -args ARGS_LNZ
codegen -config cfg pfqn_mmint2 -args ARGS_LNZ
codegen -config cfg pfqn_panacea -args ARGS_LNZ
% pfqn_kt - Skipped (calls pfqn_aql which uses variable-size cell arrays)
% codegen -config cfg pfqn_kt -args ARGS_LNZ
codegen -config cfg pfqn_xzabalow -args ARGS_LNZ
codegen -config cfg pfqn_xzabaup -args ARGS_LNZ
% pfqn_xzgsblow - Skipped (sqrt may produce complex values, codegen requires real)
% codegen -config cfg pfqn_xzgsblow -args ARGS_LNZ
% pfqn_xzgsbup - Skipped (same complex value issue as pfqn_xzgsblow)
% codegen -config cfg pfqn_xzgsbup -args ARGS_LNZ
% pfqn_propfair - Skipped (fmincon 'interior-point' not supported in codegen, only 'sqp')
% codegen -config cfg pfqn_propfair -args ARGS_LNZ
% pfqn_lap - Skipped (fzero internal nonscalar logical not supported in codegen)
% codegen -config cfg pfqn_lap -args ARGS_LNZ

%% Other simple signatures
% pfqn_grnmol - Skipped (sprod helper has uninitialized variable D in nargin>1 path)
% codegen -config cfg pfqn_grnmol -args {L_type, N_type}
% pfqn_sqni - Skipped (element count mismatch at line 54)
% codegen -config cfg pfqn_sqni -args {N_type, L_type, Z_type}
codegen -config cfg pfqn_gldsingle -args {L_type, N_type, mat_type}
codegen -config cfg pfqn_mvaldmx_ec -args {vec_type, L_type, mat_type}

%% Group: Complex Signatures

% pfqn_bs (L,N,Z,tol,maxiter,QN0,type)
codegen -config cfg pfqn_bs -args {L_type, N_type, Z_type, scal_type, scal_type, mat_type, scal_type}

% pfqn_bsfcfs (L,N,Z,tol,maxiter,QN,weight)
codegen -config cfg pfqn_bsfcfs -args {L_type, N_type, Z_type, scal_type, scal_type, mat_type, mat_type}

% pfqn_mva - Skipped (find() returns variable-size result used as scalar index)
% codegen -config cfg pfqn_mva -args {L_type, N_type, Z_type, vec_flex_type}

% pfqn_mvams - Skipped (depends on pfqn_mvamx which depends on pfqn_mva)
% codegen -config cfg pfqn_mvams -args {vec_type, L_type, N_type, Z_type, vec_flex_type, vec_flex_type}

% pfqn_mvald - Skipped (find() returns variable-size + end+1 indexing)
% codegen -config cfg pfqn_mvald -args {L_type, N_type, Z_type, mat_type, scal_type}

% pfqn_mvamx - Skipped (calls pfqn_mva which uses find() with variable-size result)
% codegen -config cfg pfqn_mvamx -args {vec_type, L_type, N_type, Z_type, vec_flex_type}

% pfqn_mvaldms (lambda,L,N,Z,S)
codegen -config cfg pfqn_mvaldms -args {vec_type, L_type, N_type, Z_type, vec_flex_type}

% pfqn_mvaldmx (lambda,L,N,Z,mu,S)
codegen -config cfg pfqn_mvaldmx -args {vec_type, L_type, N_type, Z_type, mat_type, vec_flex_type}

% pfqn_cub - Skipped (anonymous functions not supported in codegen)
% codegen -config cfg pfqn_cub -args {L_type, N_type, Z_type, scal_type, scal_type}

% pfqn_ls - Skipped (anonymous functions, mvnrnd, mvnpdf not supported in codegen)
% codegen -config cfg pfqn_ls -args {L_type, N_type, Z_type, vec_type}

% pfqn_mci - Skipped (persistent in nested function not supported in codegen)
% codegen -config cfg pfqn_mci -args {L_type, N_type, Z_type, vec_type, scal_type}

% pfqn_linearizer (L,N,Z,type,tol,maxiter)
codegen -config cfg pfqn_linearizer -args {L_type, N_type, Z_type, scal_type, scal_type, scal_type}

% pfqn_linearizerms (L,N,Z,nservers,type,tol,maxiter)
codegen -config cfg pfqn_linearizerms -args {L_type, N_type, Z_type, vec_flex_type, scal_type, scal_type, scal_type}

% pfqn_conwayms (L,N,Z,nservers,type,tol,maxiter)
codegen -config cfg pfqn_conwayms -args {L_type, N_type, Z_type, vec_flex_type, scal_type, scal_type, scal_type}

% pfqn_aql - Skipped (cell array elements cannot be verified as fully defined before use)
% codegen -config cfg pfqn_aql -args {L_type, N_type, Z_type, scal_type, scal_type, mat_type}

% pfqn_gld - Skipped (uses SolverNC.defaultOptions OOP + pfqn_nc which uses Solver.parseOptions)
% codegen -config cfg pfqn_gld -args {L_type, N_type, mat_type}

% pfqn_comom - Skipped (unique with non-fixed-size input dimension)
% codegen -config cfg pfqn_comom -args {L_type, N_type, Z_type, scal_type}

% pfqn_recal (L,N,Z,m0)
codegen -config cfg pfqn_recal -args {L_type, N_type, Z_type, vec_flex_type}

% pfqn_rd - Skipped (depends on pfqn_nc which uses OOP Solver.parseOptions)
% codegen -config cfg pfqn_rd -args {L_type, N_type, Z_type, mat_type}

% pfqn_mmint2_gausslaguerre - Skipped (gengausslegquadrule size mismatch)
% codegen -config cfg pfqn_mmint2_gausslaguerre -args {L_type, N_type, Z_type, scal_type}

% pfqn_mmint2_gausslegendre - Skipped (which() not supported in codegen)
% codegen -config cfg pfqn_mmint2_gausslegendre -args {L_type, N_type, Z_type, scal_type}

% pfqn_egflinearizer (L,N,Z,type,tol,maxiter,alpha) - alpha is vector
codegen -config cfg pfqn_egflinearizer -args {L_type, N_type, Z_type, vec_flex_type, scal_type, scal_type, vec_flex_type}

% pfqn_gflinearizer (L,N,Z,type,tol,maxiter,alpha) - alpha is scalar
codegen -config cfg pfqn_gflinearizer -args {L_type, N_type, Z_type, vec_flex_type, scal_type, scal_type, scal_type}

% pfqn_ab_amva - Skipped (containers.Map not supported in codegen)
% codegen -config cfg pfqn_ab_amva -args {L_type, N_type, L_type, vec_flex_type, vec_flex_type, scal_type, str_type}

% pfqn_linearizermx (lambda,L,N,Z,nservers,type,tol,maxiter,method)
codegen -config cfg pfqn_linearizermx -args {vec_type, L_type, N_type, Z_type, vec_flex_type, vec_flex_type, scal_type, scal_type, str_type}

%% LCFS Group (alpha,beta,N)
ARGS_LCFS = {vec_type, vec_type, vec_type};
% pfqn_lcfsqn_nc - Skipped (issym + cell array element verification)
% codegen -config cfg pfqn_lcfsqn_nc -args ARGS_LCFS
codegen -config cfg pfqn_lcfsqn_ca -args ARGS_LCFS
codegen -config cfg pfqn_lcfsqn_mva -args ARGS_LCFS

%% ===== Additional pfqn functions =====

%% Group: Load-dependent linearization
% ljd_linearize(nvec, cutoffs) -> idx
codegen -config cfg ljd_linearize -args {vec_type, vec_type}

% ljd_delinearize(idx, cutoffs) -> nvec
codegen -config cfg ljd_delinearize -args {scal_type, vec_type}

% ljcd_interpolate(nvec, cutoffs, table, K) -> Xval
codegen -config cfg ljcd_interpolate -args {vec_type, vec_type, mat_type, scal_type}

%% Group: Normalizing constant variants

% pfqn_comomrm - Skipped (end+1 + find() in for-loop)
% codegen -config cfg pfqn_comomrm -args {L_type, N_type, Z_type, scal_type, scal_type}

% pfqn_comomrm_orig - Skipped (undefined variable h_1 + complex codegen issues)
% codegen -config cfg pfqn_comomrm_orig -args {L_type, N_type, Z_type, scal_type}

% pfqn_nc_sanitize(lambda, L, N, Z, atol) -> [lambda, L, N, Z, lGremaind]
codegen -config cfg pfqn_nc_sanitize -args {vec_type, L_type, N_type, Z_type, scal_type}

% pfqn_fnc - Skipped (sparse mrdivide not supported in codegen)
% codegen -config cfg pfqn_fnc -args {vec_type, scal_type}

%% Group: Queue-dependent and probability methods

% pfqn_qd - Skipped (anonymous functions in default args not supported in codegen)
% codegen -config cfg pfqn_qd -args {L_type, N_type, vec_type, vec_type, mat_type}

% pfqn_procomom - Skipped (nonscalar logical operations not supported in codegen)
% codegen -config cfg pfqn_procomom -args {L_type, N_type, Z_type, scal_type}

% pfqn_procomom2 - Skipped (sparse + tic not supported in codegen)
% codegen -config cfg pfqn_procomom2 -args {L_type, N_type, Z_type, mat_type, vec_flex_type}

%% Group: Schmidt approximation

% pfqn_schmidt - Skipped (table() + repmat dimension change not supported in codegen)
% codegen -config cfg pfqn_schmidt -args {L_type, N_type, vec_type, vec_type, mat_type}

% pfqn_schmidt_ext - Skipped (same as pfqn_schmidt)
% codegen -config cfg pfqn_schmidt_ext -args {L_type, N_type, vec_type, vec_type}

%% Group: Result expansion/compression

% pfqn_expand(QN, UN, CN, mapping, M_original) -> [QN_full, UN_full, CN_full]
codegen -config cfg pfqn_expand -args {mat_type, mat_type, mat_type, vec_type, scal_type}

% pfqn_unique(L, mu, gamma) -> [L_unique, mu_unique, gamma_unique, mi, mapping]
codegen -config cfg pfqn_unique -args {L_type, mat_type, mat_type}
