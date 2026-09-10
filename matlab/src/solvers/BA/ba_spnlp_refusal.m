function reason = ba_spnlp_refusal(sn, method)
% REASON = BA_SPNLP_REFUSAL(SN, METHOD)
%
% The structural premises of the 'spnlp' family, in one place: the reason the
% resolved METHOD cannot bound the net SN, or '' when it can. Empty for every
% name outside the family.
%
% ONE PREDICATE, TWO CALLERS. solver_ba_spnlp_analyzer asks it before building
% the polytope and raises what it returns; BA_METHOD_REFUSAL asks it so that
% supportsModelMethod, findSolver and listValidMethods report the same
% sentence. Until it existed the mode rules lived in SPN_LPBND alone, so on a
% net with an Erlang mode model.help called 'spnlp.upper' runnable and the run
% raised; and with no rule at all a queueing network answered
% supportsModelMethod('spnlp.upper') with yes.
%
% THE RULES ARE SPN_LPBND'S OWN, restated over sn.nodeparam so that a
% predicate can answer them without assembling an LP: the model must hold
% Transition and Place nodes; a Place must not carry an embedded queue (the
% relaxation has one variable per (place, class) marking and no notion of a
% local queue); and every firing mode must be timed with a finite constant
% rate -- an IMMEDIATE mode, a marking-dependent rate (setFiringRateDependence)
% and a multi-server mode are each refused by name, as is a firing law that is
% not phase-type at all. A PHASE-TYPE law with more than one phase is where the
% two variants part: the Markovian pair 'spnlp.upper'/'spnlp.lower' needs the
% marking alone to be the state and refuses it, the operational pair
% 'spnlp.op.*' reads only its mean and admits it. SPN_LPBND keeps its own
% copies of these tests as the library's guard; the sentences here are the
% same ones it raises.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

reason = '';
switch method
    case {'spnlp.upper','spnlp.lower'}
        markovian = true;
    case {'spnlp.op.upper','spnlp.op.lower'}
        markovian = false;
    otherwise
        return
end

if ~any(sn.nodetype == NodeType.Transition)
    reason = sprintf(['Method ''%s'' bounds a stochastic Petri net; this model has no ' ...
        'Transition node. Use the queueing-network bound families, or SolverMVA/SolverNC.'], ...
        method);
    return
end
places = find(sn.nodetype == NodeType.Place);
if isempty(places)
    reason = sprintf('Method ''%s'' cannot bound this net: the model holds no Place node.', method);
    return
end
for pp = places(:)'
    ist = sn.nodeToStation(pp);
    if ist >= 1 && sn.sched(ist) ~= SchedStrategy.INF
        reason = sprintf(['Method ''%s'' does not support queueing places: place %s serves under ' ...
            '%s, and the relaxation carries one variable per (place, class) marking with no ' ...
            'notion of an embedded queue.'], method, sn.nodenames{pp}, ...
            SchedStrategy.toText(sn.sched(ist)));
        return
    end
end

nmodes = 0;
transitions = find(sn.nodetype == NodeType.Transition);
for ind = transitions(:)'
    param = sn.nodeparam{ind};
    if ~isstruct(param) || ~isfield(param,'nmodes')
        continue
    end
    for m = 1:param.nmodes
        nmodes = nmodes + 1;
        if isfield(param,'timing') && numel(param.timing) >= m && ...
                param.timing(m) == TimingStrategy.IMMEDIATE
            reason = sprintf(['Method ''%s'' cannot bound this net: mode %d of node %d is ' ...
                'IMMEDIATE; the moment relaxation is written for a net whose transitions all ' ...
                'have finite rates, so vanishing states must be eliminated first.'], method, m, ind);
            return
        end
        if isfield(param,'firingdep') && numel(param.firingdep) >= m && ~isempty(param.firingdep{m})
            reason = sprintf(['Method ''%s'' cannot bound this net: mode %d of node %d has a ' ...
                'marking-dependent firing rate; the uniformization step needs one rate per mode.'], ...
                method, m, ind);
            return
        end
        if isfield(param,'nmodeservers') && numel(param.nmodeservers) >= m && ...
                param.nmodeservers(m) ~= 1
            reason = sprintf(['Method ''%s'' cannot bound this net: mode %d of node %d has %g ' ...
                'servers; the relaxation is derived under single-server semantics, where the ' ...
                'firing rate is mu*q, and its infinite-server form is not implemented.'], ...
                method, m, ind, param.nmodeservers(m));
            return
        end
        proc = [];
        if isfield(param,'firingproc') && numel(param.firingproc) >= m
            proc = param.firingproc{m};
        end
        if ~iscell(proc) || numel(proc) < 2 || ~isnumeric(proc{1}) || ~ismatrix(proc{1}) || ...
                size(proc{1},1) ~= size(proc{1},2)
            reason = sprintf(['Method ''%s'' cannot bound this net: mode %d of node %d has a ' ...
                'firing law that is not phase-type, so neither the Markovian nor the operational ' ...
                'bound can read a mean firing rate from it. Use a phase-type law, or an exact ' ...
                'solver.'], method, m, ind);
            return
        end
        D0 = full(proc{1}); D1 = full(proc{2});
        nph = size(D0,1);
        if markovian && nph > 1
            reason = sprintf(['Method ''%s'' cannot bound this net: mode %d of node %d has a ' ...
                'phase-type firing time, and the Markovian relaxation is written over the ' ...
                'marking alone. Use ''spnlp.op.upper''/''spnlp.op.lower'', which need only the ' ...
                'mean.'], method, m, ind);
            return
        end
        if nph == 1
            rate = D1(1);
        else
            pv = [];
            if isfield(param,'firingpie') && numel(param.firingpie) >= m
                pv = param.firingpie{m};
            end
            if isempty(pv), pv = [1, zeros(1, nph-1)]; end
            pv = pv(:)' / sum(pv);
            rate = 1 / (pv * ((-D0) \ ones(nph,1)));
        end
        if ~(rate > 0) || ~isfinite(rate)
            reason = sprintf(['Method ''%s'' cannot bound this net: mode %d of node %d has mean ' ...
                'firing rate %g; a bound needs a finite positive one.'], method, m, ind, rate);
            return
        end
    end
end
if nmodes == 0
    reason = sprintf('Method ''%s'' cannot bound this net: the net has no firing mode.', method);
end
end
