function reason = qns_immfeed_refusal(sn)
% REASON = QNS_IMMFEED_REFUSAL(SN)
% Why SolverQNS cannot serve a model with immediate feedback, or '' when the
% model has none.
%
% Immediate feedback (sn.immfeed) keeps a self-looping job on its server
% instead of re-queueing it, and neither path of SolverQNS can state that: the
% JMVA document qnsolver reads carries a mean demand and a visit count per
% chain, and the LQN QN2LQN writes turns the routing into OR-fork precedences
% of pseudo-activities on the reference task, where a repeated visit is a new
% call. Either would answer for re-queueing under this solver's name.
%
% ONE PREDICATE, TWO CALLERS: SolverQNS.supportsModelMethod (the gate, hence
% model.help and SolverAUTO) and SolverQNS.runAnalyzer (the run, for a caller
% with enableChecks off). SolverJMT keeps its own wording in jmtMethodRefusal.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

reason = '';
if isfield(sn,'immfeed') && ~isempty(sn.immfeed) && any(sn.immfeed(:))
    reason = ['SolverQNS does not support immediate feedback (sn.immfeed): neither the JMVA ' ...
        'document qnsolver reads nor the LQN QN2LQN writes can keep a self-looping job on ' ...
        'its server. Use SolverCTMC or SolverSSA, whose state space carries the self-loop.'];
end
end
