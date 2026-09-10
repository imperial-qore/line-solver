function traits = solverTraits(self)
% TRAITS = SOLVERTRAITS()
%
% Structural traits of the model that drive solver ranking. Computed once
% per call from the NetworkStruct so that chooseSolver* stays a table of
% rankings rather than a second place where model inspection is written.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

traits = struct('hasCache', false, 'hasFCR', false, 'hasFork', false, ...
    'hasMAP', false, 'prio', 'none', 'isClosed', false, 'isOpen', false, ...
    'isMixed', false, 'hasMultiServer', false, 'isProductForm', false, ...
    'popPerChain', 0, 'totalJobs', 0, 'singleChain', false);

if ~isa(self.model, 'Network')
    return
end

this_model = self.model;
sn = this_model.getStruct(false);

traits.hasCache = any(sn.nodetype == NodeType.Cache);
traits.hasFCR = ~isempty(sn.nregions) && sn.nregions > 0;
traits.hasFork = this_model.hasFork;
traits.hasMultiServer = this_model.hasMultiServer;
traits.isProductForm = this_model.hasProductFormSolution();
traits.singleChain = this_model.hasSingleChain;

hasOpen = this_model.hasOpenClasses;
hasClosed = this_model.hasClosedClasses;
traits.isOpen = hasOpen && ~hasClosed;
traits.isClosed = hasClosed && ~hasOpen;
traits.isMixed = hasOpen && hasClosed;

njobs = this_model.getNumberOfJobs();
traits.totalJobs = sum(njobs(isfinite(njobs)));
nchains = sum(this_model.getNumberOfChains);
if nchains > 0
    traits.popPerChain = traits.totalJobs / nchains;
end

% Autocorrelated arrival or service: only MAM keeps the correlation, every
% other analytical solver sees the marginal only.
traits.hasMAP = any(sn.procid(:) == ProcessType.MAP | sn.procid(:) == ProcessType.MMPP2);

% Preemptive priority is a strictly narrower capability than HOL, so the two
% rank differently and must be distinguished here.
sched = sn.sched;
if any(sched == SchedStrategy.FCFSPRPRIO) || any(sched == SchedStrategy.FCFSPIPRIO) || ...
        any(sched == SchedStrategy.LCFSPRPRIO) || any(sched == SchedStrategy.LCFSPIPRIO)
    traits.prio = 'preempt';
elseif any(sched == SchedStrategy.PSPRIO) || any(sched == SchedStrategy.DPSPRIO) || ...
        any(sched == SchedStrategy.GPSPRIO)
    traits.prio = 'ps';
elseif any(sched == SchedStrategy.HOL) || any(sched == SchedStrategy.LCFSPRIO) || ...
        any(sched == SchedStrategy.SRPTPRIO)
    traits.prio = 'hol';
end
end
