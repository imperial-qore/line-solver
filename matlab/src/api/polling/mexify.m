% @brief MATLAB Coder script to generate MEX functions for polling_ module.
%
% This script generates MEX (MATLAB Executable) versions of polling-system
% functions for improved performance.
%
% Skipped functions (Coder-incompatible):
%   polling_qsys_1limited    - Cell-array-of-MAPs arguments (each MAP is a
%                              cell of two matrices), unsupported by Coder.
%   polling_qsys_exhaustive  - Same: cell-array-of-MAPs arguments.
%   polling_qsys_gated       - Same: cell-array-of-MAPs arguments.
%
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.CodeConfig'.
cfg = coder.config('mex','ecoder',false);
cfg.GenerateReport = false;
cfg.ReportPotentialDifferences = false;
cfg.GenCodeOnly = false;

%% No codegen-compatible functions in this domain.
