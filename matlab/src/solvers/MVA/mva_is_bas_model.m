function tf = mva_is_bas_model(sn)
% MVA_IS_BAS_MODEL Detect a closed single-chain network with Blocking-After-Service
% (BAS) finite-buffer blocking, which solver_sqd handles but exact/AMVA MVA does not.
%
% Copyright (c) 2012-2026, QORE Lab, Imperial College London
% All rights reserved.

tf = false;
if sn.nchains ~= 1 || sn.nclosedjobs <= 0 || isempty(sn.droprule)
    return;
end
if any(isinf(sn.njobs))
    return; % open class present
end
tf = any(sn.droprule(:) == DropStrategy.BAS);
end
