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
codegen -config cfg pfqn_kt -args ARGS_LNZ
codegen -config cfg pfqn_xzabalow -args ARGS_LNZ
codegen -config cfg pfqn_xzabaup -args ARGS_LNZ
codegen -config cfg pfqn_xzgsblow -args ARGS_LNZ
codegen -config cfg pfqn_xzgsbup -args ARGS_LNZ
codegen -config cfg pfqn_propfair -args ARGS_LNZ
codegen -config cfg pfqn_lap -args ARGS_LNZ

%% Other simple signatures
codegen -config cfg pfqn_grnmol -args {L_type, N_type}
codegen -config cfg pfqn_sqni -args {N_type, L_type, Z_type}
codegen -config cfg pfqn_gldsingle -args {L_type, N_type, mat_type}
codegen -config cfg pfqn_mvaldmx_ec -args {vec_type, L_type, mat_type}

%% Group: Complex Signatures

% pfqn_bs (L,N,Z,tol,maxiter,QN0,type)
codegen -config cfg pfqn_bs -args {L_type, N_type, Z_type, scal_type, scal_type, mat_type, scal_type}

% pfqn_bsfcfs (L,N,Z,tol,maxiter,QN,weight)
codegen -config cfg pfqn_bsfcfs -args {L_type, N_type, Z_type, scal_type, scal_type, mat_type, mat_type}

% pfqn_mva (L,N,Z,mi)
codegen -config cfg pfqn_mva -args {L_type, N_type, Z_type, vec_flex_type}

% pfqn_mvams (lambda,L,N,Z,mi,S)
codegen -config cfg pfqn_mvams -args {vec_type, L_type, N_type, Z_type, vec_flex_type, vec_flex_type}

% pfqn_mvald (L,N,Z,mu,stabilize)
codegen -config cfg pfqn_mvald -args {L_type, N_type, Z_type, mat_type, scal_type}

% pfqn_mvamx (lambda,L,N,Z,mi)
codegen -config cfg pfqn_mvamx -args {vec_type, L_type, N_type, Z_type, vec_flex_type}

% pfqn_mvaldms (lambda,L,N,Z,S)
codegen -config cfg pfqn_mvaldms -args {vec_type, L_type, N_type, Z_type, vec_flex_type}

% pfqn_mvaldmx (lambda,L,N,Z,mu,S)
codegen -config cfg pfqn_mvaldmx -args {vec_type, L_type, N_type, Z_type, mat_type, vec_flex_type}

% pfqn_cub (L,N,Z,order,atol)
codegen -config cfg pfqn_cub -args {L_type, N_type, Z_type, scal_type, scal_type}

% pfqn_ls (L,N,Z,I)
codegen -config cfg pfqn_ls -args {L_type, N_type, Z_type, vec_type}

% pfqn_mci (D,N,Z,I,variant)
codegen -config cfg pfqn_mci -args {L_type, N_type, Z_type, vec_type, scal_type}

% pfqn_linearizer (L,N,Z,type,tol,maxiter)
codegen -config cfg pfqn_linearizer -args {L_type, N_type, Z_type, scal_type, scal_type, scal_type}

% pfqn_linearizerms (L,N,Z,nservers,type,tol,maxiter)
codegen -config cfg pfqn_linearizerms -args {L_type, N_type, Z_type, vec_flex_type, scal_type, scal_type, scal_type}

% pfqn_conwayms (L,N,Z,nservers,type,tol,maxiter)
codegen -config cfg pfqn_conwayms -args {L_type, N_type, Z_type, vec_flex_type, scal_type, scal_type, scal_type}

% pfqn_aql (L,N,Z,TOL,MAXITER,QN0)
codegen -config cfg pfqn_aql -args {L_type, N_type, Z_type, scal_type, scal_type, mat_type}

% pfqn_gld (L,N,mu) - Options skipped
codegen -config cfg pfqn_gld -args {L_type, N_type, mat_type}

% pfqn_comom (L,N,Z,atol)
codegen -config cfg pfqn_comom -args {L_type, N_type, Z_type, scal_type}

% pfqn_recal (L,N,Z,m0)
codegen -config cfg pfqn_recal -args {L_type, N_type, Z_type, vec_flex_type}

% pfqn_rd (L,N,Z,mu)
codegen -config cfg pfqn_rd -args {L_type, N_type, Z_type, mat_type}

% pfqn_mmint2_gausslaguerre (L,N,Z,m)
codegen -config cfg pfqn_mmint2_gausslaguerre -args {L_type, N_type, Z_type, scal_type}

% pfqn_mmint2_gausslegendre (L,N,Z,m)
codegen -config cfg pfqn_mmint2_gausslegendre -args {L_type, N_type, Z_type, scal_type}

% pfqn_egflinearizer (L,N,Z,type,tol,maxiter,alpha) - alpha is vector
codegen -config cfg pfqn_egflinearizer -args {L_type, N_type, Z_type, vec_flex_type, scal_type, scal_type, vec_flex_type}

% pfqn_gflinearizer (L,N,Z,type,tol,maxiter,alpha) - alpha is scalar
codegen -config cfg pfqn_gflinearizer -args {L_type, N_type, Z_type, vec_flex_type, scal_type, scal_type, scal_type}

% pfqn_ab_amva (D,N,V,nservers,sched,fcfsSchmidt,marginalProbMethod)
codegen -config cfg pfqn_ab_amva -args {L_type, N_type, L_type, vec_flex_type, vec_flex_type, scal_type, str_type}

% pfqn_linearizermx (lambda,L,N,Z,nservers,type,tol,maxiter,method)
codegen -config cfg pfqn_linearizermx -args {vec_type, L_type, N_type, Z_type, vec_flex_type, vec_flex_type, scal_type, scal_type, str_type}

%% LCFS Group (alpha,beta,N)
ARGS_LCFS = {vec_type, vec_type, vec_type};
codegen -config cfg pfqn_lcfsqn_nc -args ARGS_LCFS
codegen -config cfg pfqn_lcfsqn_ca -args ARGS_LCFS
codegen -config cfg pfqn_lcfsqn_mva -args ARGS_LCFS