function [ok, reason] = fluid_petri_applicable(sn, options)
% [OK, REASON] = FLUID_PETRI_APPLICABLE(SN, OPTIONS)
%
% Whether the fluid Petri route can answer this model, and why not when it
% cannot.
%
% The route covers the class of nets the exact engines cover -- Place,
% Transition, enabling, inhibiting and firing arcs, timed and immediate modes,
% single/k/infinite server firing, phase-type and MAP firing times,
% marking-dependent firing rates, bounded places, and a Source/Sink pair for an
% open net. What it does NOT cover is named here rather than discovered inside
% the solve, so a refused model says which declaration it was refused for.
%
% A QUEUEING STATION IS THE ONE STRUCTURAL EXCLUSION. A net whose tokens also
% visit a Queue or a Delay is two formalisms at once, and LINE has no reference
% semantics for the hand-off: the Petri arcs and the routing matrix would each
% describe part of the movement and nothing pins how a token becomes a job. The
% same goes for a QUEUEING PLACE, whose embedded queue only SolverLDES
% simulates.
%
% Parameters:
%   sn      - NetworkStruct
%   options - solver options
%
% Returns:
%   ok     - true when the model can be solved by SOLVER_FLUID_PETRI
%   reason - the refusal, empty when OK
%
% See also SOLVER_FLUID_PETRI, FLUID_DAE_APPLICABLE.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

ok = false;
reason = '';

if ~any(sn.nodetype == NodeType.Transition)
    reason = 'the model has no Transition node, so it is not a Petri net';
    return
end

allowed = [NodeType.Place, NodeType.Transition, NodeType.Source, NodeType.Sink];
for ind = 1:sn.nnodes
    if ~any(sn.nodetype(ind) == allowed)
        reason = sprintf(['node %s is a %s. The fluid Petri route solves the marking of a Petri net, and a ' ...
            'model that also holds queueing stations is two formalisms at once with no reference semantics ' ...
            'for the hand-off; use SolverCTMC, SolverJMT, SolverSSA or SolverLDES'], ...
            sn.nodenames{ind}, NodeType.toText(sn.nodetype(ind)));
        return
    end
end

% A queueing place declares a service process, which is what turns it into a
% station with an embedded queue and a depository.
for ind = find(sn.nodetype == NodeType.Place)'
    ist = sn.nodeToStation(ind);
    for k = 1:sn.nclasses
        if ~isnan(sn.rates(ist,k)) && sn.rates(ist,k) > 0
            reason = sprintf(['place %s is a QUEUEING place (it declares a service process), whose embedded ' ...
                'queue this drift does not carry; use SolverLDES'], sn.nodenames{ind});
            return
        end
    end
end

% options.config.hide_immediate is NOT consulted here. It is on by default now,
% and this route never builds the event set ODE_ELIMINATE_IMMEDIATE acts on: it
% carries immediate firings as algebraic unknowns instead, which is the stronger
% treatment. FLUID_HIDE_IMMEDIATE answers false for a Petri net for that reason,
% so the flag reaches nothing here and refusing on it would refuse every SPN.

ok = true;
end
