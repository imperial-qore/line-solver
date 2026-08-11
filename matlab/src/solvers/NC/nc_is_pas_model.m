function ispas = nc_is_pas_model(sn)
% ISPAS = NC_IS_PAS_MODEL(SN)
%
% True when the model is a closed two-station pass-and-swap (P&S) tandem: both
% stations are OI/PAS (SchedStrategy.OI or SchedStrategy.PAS) and there is no
% other station. The swap graph may be empty, since an order-independent (OI)
% queue is exactly the P&S specialization with an empty/zero swap graph: the
% importance-sampling analyzer SOLVER_NC_PAS_IS_ANALYZER covers both, as
% PFQN_PAS_IS with H=0 samples all microstates and reduces to PFQN_OI_IS.
%
% Used to bind the importance-sampling selectors ('is', and 'sampling' which
% maps to 'is' in the presence of OI/PAS stations). Note that a P&S tandem with
% a NON-EMPTY swap graph is reducible (Comte & Dorsman, 2021) and lies outside
% the exact path of NC_IS_OI_MODEL / SOLVER_NC_OI_ANALYZER, so it also takes
% this analyzer on 'default'; a pure-OI tandem on 'default'/'exact' is caught
% earlier by the exact OI analyzer.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
ispas = false;
if any(isinf(sn.njobs))
    return
end
if sn.nstations ~= 2
    return
end
for ist = 1:sn.nstations
    if sn.sched(ist) ~= SchedStrategy.PAS && sn.sched(ist) ~= SchedStrategy.OI
        return
    end
    ind = sn.stationToNode(ist);
    if ind < 1 || ind > numel(sn.nodeparam) || ~isstruct(sn.nodeparam{ind})
        return
    end
    if isempty(sn.nodeparam{ind}.svcRateFun)
        return   % no OI rank-rate function: cannot evaluate the balance function
    end
end
ispas = true;
end
