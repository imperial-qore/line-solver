function bool = sn_has_joint_dependence(sn)
% BOOL = SN_HAS_JOINT_DEPENDENCE(SN)
%
% Checks if the network has joint-dependent service rates
% This includes both LJD (Limited Joint Dependence) and LJCD
% (Limited Joint Class Dependence) scaling.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

hasLjd = ~isempty(sn.ljdscaling) && any(~cellfun(@isempty, sn.ljdscaling));
hasLjcd = ~isempty(sn.ljcdscaling) && any(~cellfun(@isempty, sn.ljcdscaling));
bool = hasLjd || hasLjcd;
end
