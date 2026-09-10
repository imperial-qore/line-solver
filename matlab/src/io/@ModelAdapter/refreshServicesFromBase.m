function ok = refreshServicesFromBase(nonfjmodel, prov)
% OK = REFRESHSERVICESFROMBASE(NONFJMODEL, PROV)
%
% Re-feed a transformed model produced by ModelAdapter.mmt from the current
% service parameters of the base model it was derived from, so that the
% transformation can be reused across the iterations of an outer fixed point
% instead of rebuilt. SolverLN re-solves each layer once per iteration and only
% the rates change between iterations; the fork topology the transformation
% encodes does not. Rebuilding costs a full model.copy() every time.
%
% PROV is the provenance struct returned by ModelAdapter.mmt. Three kinds of
% slot, and all three must be handled differently:
%   prov.serviceSrc     base-derived  -> re-read from the base model
%   prov.immediateSlots transformation-owned -> RESET to Immediate. The fork
%                       loop in SolverMVA.runAnalyzer overwrites every join with
%                       its current synchronisation delay on each pass, so a
%                       reused model still holds the previous outer iteration's
%                       converged value. Leaving these alone silently warm-starts
%                       the fork loop and shifts the results.
%   prov.auxArrival     owned by the forkLambda fixed point -> reset to what a
%                       cold mmt call would have set.
% Reading a transformation-owned slot from the base model is the converse error
% and corrupts the transform outright.
%
% Returns false when the provenance cannot be applied, in which case the caller
% must fall back to a cold ModelAdapter.mmt.

ok = false;
if isempty(prov) || ~isstruct(prov) || ~isfield(prov,'baseModel') || isempty(prov.baseModel)
    return
end
base = prov.baseModel;
nclasses = length(base.classes);

% Base-derived service slots.
for k = 1:size(prov.serviceSrc,1)
    i = prov.serviceSrc(k,1); c = prov.serviceSrc(k,2); r = prov.serviceSrc(k,3);
    if i > length(nonfjmodel.nodes) || c > length(nonfjmodel.classes) || r > nclasses
        return
    end
    svc = base.nodes{i}.getService(base.classes{r});
    if isempty(svc)
        return
    end
    % copy() mirrors the cold path, which hands every slot its own distribution
    % object rather than aliasing one across classes.
    nonfjmodel.nodes{i}.setService(nonfjmodel.classes{c}, svc.copy());
end

% Transformation-owned initial conditions.
for k = 1:size(prov.immediateSlots,1)
    i = prov.immediateSlots(k,1); c = prov.immediateSlots(k,2);
    if i > length(nonfjmodel.nodes) || c > length(nonfjmodel.classes)
        return
    end
    nonfjmodel.nodes{i}.setService(nonfjmodel.classes{c}, Immediate());
end

% Auxiliary-class arrivals, as a cold call would leave them.
if ~isempty(prov.auxArrival)
    source = nonfjmodel.getSource;
    if isempty(source)
        return
    end
    for k = 1:size(prov.auxArrival,1)
        c = prov.auxArrival(k,1); r = prov.auxArrival(k,2); disableAux = prov.auxArrival(k,3);
        if c > length(nonfjmodel.classes)
            return
        end
        if disableAux
            source.setArrival(nonfjmodel.classes{c}, Disabled.getInstance);
        elseif r <= length(prov.forkLambdaInit)
            source.setArrival(nonfjmodel.classes{c}, Exp(prov.forkLambdaInit(r)));
        end
    end
end
ok = true;
end
