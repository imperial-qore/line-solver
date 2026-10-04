function isoi = nc_is_oi_model(sn)
% ISOI = NC_IS_OI_NC_MODEL(SN)
%
% True when the model is a closed queueing network that contains at least one
% order-independent (OI) station (SchedStrategy.OI, or PAS with an empty swap
% graph) and every other station is a BCMP product-form station: infinite
% server (delay), PS, LCFS-PR, SIRO, or class-independent-rate FCFS. Such
% models are solved exactly by SOLVER_NC_OI_NC_ANALYZER. The OI requirement
% keeps pure-BCMP networks on the standard (faster) normalizing-constant path.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
isoi = false;
if any(isinf(sn.njobs))
    return
end
hasOI = false;
for ist = 1:sn.nstations
    if sn.sched(ist) == SchedStrategy.INF
        continue
    elseif sn.sched(ist) == SchedStrategy.PAS || sn.sched(ist) == SchedStrategy.OI
        ind = sn.stationToNode(ist);
        if ind < 1 || ind > numel(sn.nodeparam) || ~isstruct(sn.nodeparam{ind}) ...
                || ~isfield(sn.nodeparam{ind}, 'swapGraph')
            return
        end
        sg = sn.nodeparam{ind}.swapGraph;
        if isempty(sg) || any(sg(:) ~= 0)
            return   % genuine pass-and-swap: not order-independent
        end
        hasOI = true;
    elseif any(sn.sched(ist) == [SchedStrategy.PS, SchedStrategy.LCFSPR, SchedStrategy.SIRO, SchedStrategy.FCFS])
        if any(sn.sched(ist) == [SchedStrategy.FCFS, SchedStrategy.SIRO])
            % Product form requires a class-independent FCFS/SIRO rate.
            rr = sn.rates(ist, :);
            rr = rr(isfinite(rr) & sn.njobs > 0);
            if ~isempty(rr) && (max(rr) - min(rr)) > 1e-9 * max(rr)
                return   % class-dependent FCFS/SIRO: not product form
            end
        end
    else
        return       % an unsupported (non-product-form) station
    end
end
isoi = hasOI;
end
