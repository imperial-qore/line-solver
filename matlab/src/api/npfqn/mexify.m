% @brief MATLAB Coder script to generate MEX functions for npfqn_ module.
%
% This script generates MEX (MATLAB Executable) versions of npfqn_
% (Non-Product-Form Queueing Network) functions for improved performance.
%
% Skipped functions (Coder-incompatible):
%   npfqn_nonexp_approx       - Uses sn struct + SchedStrategy /
%                               GlobalConstants OOP enums and dispatches on
%                               string method, none of which is supported
%                               by Coder.
%   npfqn_traffic_merge       - Cell-array-of-MMAPs argument and dynamic
%                               struct dispatching on config.merge /
%                               config.compress strings.
%   npfqn_traffic_merge_cs    - Cell-array-of-MMAPs argument with
%                               dynamic-typed config struct.
%   npfqn_traffic_split_cs    - Builds variable-shape cell-of-cell output
%                               (varargout) which Coder cannot infer.
%
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.CodeConfig'.
cfg = coder.config('mex','ecoder',false);
cfg.GenerateReport = false;
cfg.ReportPotentialDifferences = false;
cfg.GenCodeOnly = false;

%% No codegen-compatible functions in this domain.
