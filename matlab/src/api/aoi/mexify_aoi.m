% @brief MATLAB Coder script to generate MEX functions for aoi_ module.
%
% This script generates MEX (MATLAB Executable) versions of Age of
% Information (AoI) metric functions for improved performance.
% Only pure scalar-in/scalar-out functions are included.
%
% Skipped functions (Coder-incompatible):
%   aoi_lst_det, aoi_lst_exp, aoi_lst_erlang, aoi_lst_ph
%       - Return function handles
%   aoi_fcfs_mgi1, aoi_fcfs_gim1, aoi_lcfspr_mgi1, aoi_lcfspr_gim1
%       - Take function handle (LST) arguments
%   aoi_lcfss_mgi1, aoi_lcfss_gim1, aoi_lcfsd_mgi1, aoi_lcfsd_gim1
%       - Take function handle (LST) arguments
%   aoi_is_aoi        - SchedStrategy/NodeType enum dependencies
%   aoi_extract_params - SchedStrategy, sn.proc cell array
%   aoi_dist2ph        - Cell array unpacking
%
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.CodeConfig'.
cfg = coder.config('mex','ecoder',false);
cfg.GenerateReport = false;
cfg.ReportPotentialDifferences = false;
cfg.GenCodeOnly = false;

%% Define common types
scal_type = coder.typeof(0);

%% ===== Group: FCFS AoI metrics (scalar) =====

% aoi_fcfs_mm1(lambda, mu) -> [meanAoI, varAoI, peakAoI]
codegen -config cfg aoi_fcfs_mm1 -args {scal_type, scal_type}

% aoi_fcfs_md1(lambda, d) -> [meanAoI, varAoI, peakAoI]
codegen -config cfg aoi_fcfs_md1 -args {scal_type, scal_type}

% aoi_fcfs_dm1(tau, mu) -> [meanAoI, varAoI, peakAoI]
codegen -config cfg aoi_fcfs_dm1 -args {scal_type, scal_type}

%% ===== Group: LCFS-PR AoI metrics (scalar) =====

% aoi_lcfspr_mm1(lambda, mu) -> [meanAoI, varAoI, peakAoI]
codegen -config cfg aoi_lcfspr_mm1 -args {scal_type, scal_type}

% aoi_lcfspr_md1(lambda, d) -> [meanAoI, varAoI, peakAoI]
codegen -config cfg aoi_lcfspr_md1 -args {scal_type, scal_type}

% aoi_lcfspr_dm1(tau, mu) -> [meanAoI, varAoI, peakAoI]
codegen -config cfg aoi_lcfspr_dm1 -args {scal_type, scal_type}
