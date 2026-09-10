function [hasRef, refstat_k, refclass_c] = ln_layer_refcell(layerSn, classidx)
% [HASREF, REFSTAT_K, REFCLASS_C] = LN_LAYER_REFCELL(LAYERSN, CLASSIDX)
%
% Locate the (station,class) cell of the chain reference throughput used to
% normalise a residence time in a layer. HASREF is false when the chain of
% CLASSIDX carries no reference class, as on an open chain, where
% LAYERSN.REFCLASS is 0 and the cell does not exist; the caller must then
% fall back on WN rather than index TN.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

hasRef = false;
refstat_k = 0;
refclass_c = 0;

c = find(layerSn.chains(:, classidx), 1);
if isempty(c)
    return
end
refclass_c = layerSn.refclass(c);
refstat_k = layerSn.refstat(classidx);
hasRef = refclass_c > 0 && refstat_k > 0;
end
