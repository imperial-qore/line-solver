function [hasRef, refstat_k, refclass_c] = ln_layer_refcell(layerSn, classidx, normclass)
% [HASREF, REFSTAT_K, REFCLASS_C] = LN_LAYER_REFCELL(LAYERSN, CLASSIDX, NORMCLASS)
%
% Locate the (station,class) cell of the chain reference throughput used to
% normalise a residence time in a layer. HASREF is false when the chain of
% CLASSIDX carries no reference class, as on an open chain, where
% LAYERSN.REFCLASS is 0 and the cell does not exist; the caller must then
% fall back on WN rather than index TN.
%
% NORMCLASS is optional and names, per class, the TASK class whose invocation
% rate that class is measured against. It exists because a chain may hold more
% than one caller task: the refpath interlock method merges the client chains of
% every caller descending from one reference task, and refreshStruct allows at
% most ONE reference class per chain, so the other callers cannot be found
% through LAYERSN.REFCLASS any more. Their residence times are still "per
% invocation of the owning task", which is what NORMCLASS records. Without it,
% and with the chain's single reference class standing in, every merged caller's
% residence would be normalised by the head's rate and come out scaled by the
% product of the call multiplicities along the path.
%
% Omitted, or 0 for this class, it falls back to the chain's reference class,
% which is the class LAYERSN.REFCLASS would have named anyway when no merging
% took place. The argument is therefore inert on every other method.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

hasRef = false;
refstat_k = 0;
refclass_c = 0;

refstat_k = layerSn.refstat(classidx);
if nargin >= 3 && ~isempty(normclass) && classidx <= numel(normclass) && normclass(classidx) > 0
    refclass_c = normclass(classidx);
    hasRef = refclass_c > 0 && refstat_k > 0;
    return
end

c = find(layerSn.chains(:, classidx), 1);
if isempty(c)
    refstat_k = 0;
    return
end
refclass_c = layerSn.refclass(c);
hasRef = refclass_c > 0 && refstat_k > 0;
end
