function reason = nc_sdr_refusal(sn)
% REASON = NC_SDR_REFUSAL(SN)
%
% May SOLVER_NC_SDR_ANALYZER run on this model? '' when it may, otherwise the
% reason it may not, in the words the analyzer refuses with. A model that
% declares no state-dependent routing is not that route's business and gets ''.
%
% The product form of Krzesinski (1987), "Multiclass Queueing Networks with
% State-Dependent Routing", Performance Evaluation 7:125-143, eq. (16), is
% stated over the queue lengths of a CLOSED network whose customers keep their
% class, at BCMP centres whose f_i(n_i) admits chain-dependent rates only at the
% symmetric disciplines. Each clause below is one of those premises.
%
% ONE PREDICATE, TWO CALLERS. The analyzer asks it first and turns a non-empty
% answer into an error; NC_METHOD_REFUSAL asks it so that MODEL.HELP never
% offers a pair the run refuses. A new rule goes here, not at a call site.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

reason = '';
if ~isfield(sn,'sdr') || isempty(sn.sdr)
    return
end
M = sn.nstations;
R = sn.nclasses;
if any(isinf(sn.njobs))
    reason = 'State-dependent routing is defined for closed networks only; the model has an open class.';
    return
end
if sn.nchains ~= R
    reason = ['State-dependent routing does not support class switching: the product form of Krzesinski (1987) ', ...
        'is stated over closed chains whose customers keep their class. Merge the switching classes into one class.'];
    return
end
if sn.nstateful ~= M
    reason = ['State-dependent routing requires every stateful node to be a station: the product form is over ', ...
        'queue lengths, and a stateless node holds none.'];
    return
end
% A BCMP center served FCFS must hold one rate for every chain; the paper's
% f_i(n_i) admits chain-dependent rates only at the symmetric disciplines.
S = zeros(M,R);
for i = 1:M
    for r = 1:R
        if sn.rates(i,r) > 0 && isfinite(sn.rates(i,r))
            S(i,r) = 1/sn.rates(i,r);
        end
    end
end
for i = 1:M
    if sn.sched(i) == SchedStrategy.FCFS
        act = find(sn.njobs > 0 & any(S(i,:) > 0, 1));
        if numel(unique(S(i,act))) > 1
            reason = sprintf(['Station %s is FCFS with chain-dependent service times, which has no BCMP product ', ...
                'form. Use PS, LCFSPR or INF, or equalize the service times.'], sn.nodenames{sn.stationToNode(i)});
            return
        end
    end
end
end
