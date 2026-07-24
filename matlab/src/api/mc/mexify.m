% @brief MATLAB Coder script to generate MEX functions for mc_ module.
%
% This script generates MEX (MATLAB Executable) versions of Markov chain
% functions for improved performance.
%
% Skipped functions (Coder-incompatible):
%   ctmc_courtois, ctmc_kms         - cell array args (MS), dtmc_solve dependency
%   ctmc_takahashi                  - cell array args, dtmc_solve, ctmc_randomization
%   ctmc_multi                      - cell array args, ctmc_solve_reducible
%   ctmc_simulate                   - exprnd() not supported
%   ctmc_uniformization             - speye(), sparse() not supported
%   ctmc_transient                  - function handles (@ctmc_transientode)
%   ctmc_ssg, ctmc_ssg_reachability - State class OOP methods
%   dtmc_solve_reducible            - stronglyconncomp() not supported
%   ctmc_solve_reducible            - depends on dtmc_solve_reducible
%   ctmc_solve_reducible_blkdecomp  - stronglyconncomp() not supported
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

%% ===== Infinitesimal generator =====

% ctmc_makeinfgen(Q) -> Q
codegen -config cfg ctmc_makeinfgen -args {mat_type}

% ctmc_stochcomp(Q, I) -> [S, Q11, Q12, Q21, Q22, T]
codegen -config cfg ctmc_stochcomp -args {mat_type, coder.typeof(0,[1 Inf],[0 1])}

%% ===== Random generator =====

% ctmc_rand(n) -> Q
codegen -config cfg ctmc_rand -args {scal_type}
