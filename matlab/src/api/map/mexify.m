% @brief MATLAB Coder script to generate MEX functions for map_ module.
%
% This script generates MEX (MATLAB Executable) versions of MAP (Markovian
% Arrival Process) functions for improved performance.
%
% Skipped functions (Coder-incompatible):
%   map_m1ps_h_recursive  - Returns a 2-D cell array of variable-size
%                           vectors. Coder cannot infer cell-element types
%                           constructed inside nested loops with @cell().
%   map_m1ps_sojourn      - Uses inputParser/varargin, function handles in
%                           validators, dynamic cell sizing, and the cell
%                           output of map_m1ps_h_recursive.
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

%% ===== Group: MAP rate-matrix computation =====

% map_compute_R(C, D, mu) -> R
codegen -config cfg map_compute_R -args {mat_type, mat_type, scal_type}
