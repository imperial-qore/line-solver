% @brief MATLAB Coder script to generate MEX functions for fes_ module.
%
% This script generates MEX (MATLAB Executable) versions of FES (Flow
% Equivalent Server) functions for improved performance.
%
% Skipped functions (Coder-incompatible):
%   fes_build_isolated       - Uses sn struct + NodeType / GlobalConstants
%                              OOP enums; relies on dynamic typing.
%   fes_compute_throughputs  - Uses a dynamic cell-array BFS queue,
%                              try/catch around pfqn_mva and a
%                              cell-of-vector output.
%   fes_validate             - Uses sn struct + NodeType OOP enum and
%                              dynamic isfield()/sprintf() error reporting.
%
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.CodeConfig'.
cfg = coder.config('mex','ecoder',false);
cfg.GenerateReport = false;
cfg.ReportPotentialDifferences = false;
cfg.GenCodeOnly = false;

%% No codegen-compatible functions in this domain.
