% @brief MATLAB Coder script to generate MEX functions for lsn_ module.
%
% This script generates MEX (MATLAB Executable) versions of LSN (Layered
% Software Network) functions for improved performance.
%
% Skipped functions (Coder-incompatible):
%   lsn_max_multiplicity  - Uses LayeredNetworkElement OOP enum and a
%                           heterogeneous lsn struct with cell-array
%                           field 'arrival', neither of which is
%                           supported by MATLAB Coder.
%
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.CodeConfig'.
cfg = coder.config('mex','ecoder',false);
cfg.GenerateReport = false;
cfg.ReportPotentialDifferences = false;
cfg.GenCodeOnly = false;

%% No codegen-compatible functions in this domain.
