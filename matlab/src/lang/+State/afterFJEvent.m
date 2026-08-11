function [outGlobalStates, outrate, outprob] = afterFJEvent(sn, fjentry, glspace, isSimulation)
% [OUTGLOBALSTATES, OUTRATE, OUTPROB] = AFTERFJEVENT(SN, FJENTRY, GLSPACE, ISSIMULATION)
%
% Fire a fork synchronization (one entry of sn.fjsync, built by
% ModelAdapter.fjtag) on the global state GLSPACE, a cell array with the
% local state of every stateful node. The firing atomically consumes one
% parent job of class r held at the (stateful) Fork node and emits one
% sibling per branch, in the auxiliary classes of the entry's tag, at the
% branch head nodes.
%
% Enabling condition (evaluated on the global state):
%  1. the Fork holds at least one class-r parent job;
%  2. the entry's tag is the LOWEST free tag for this (fork, class):
%     a tag is free iff its auxiliary classes have zero occupancy
%     network-wide. Canonical lowest-free-tag allocation ensures exactly
%     one fjsync entry per (fork, class) is enabled in any state.
%
% Outputs: OUTGLOBALSTATES is a cell array of successor global states
% (full glspace cells); OUTRATE the firing rate (immediate) and OUTPROB
% the probability of each outcome (phase-entry mixtures of the sibling
% service processes at the branch heads). When ISSIMULATION is true a
% single outcome is returned, sampled internally.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

outGlobalStates = {};
outrate = [];
outprob = [];

R = sn.nclasses;
f = fjentry.fork;
r = fjentry.class;
isf_f = sn.nodeToStateful(f);

% 1. parent job held at the fork
forkstate = glspace{isf_f};
if forkstate(end-R+r) < 1
    return
end

% 2. lowest-free-tag test: occupancy of each tag's auxiliary classes,
% scanned across all stateful nodes
auxall = fjentry.auxall; % B x T auxiliary class indices
T = size(auxall,2);
t = fjentry.tag;
nglobal = zeros(1,R);
for isf=1:sn.nstateful
    [~, nir] = State.toMarginalAggr(sn, sn.statefulToNode(isf), glspace{isf});
    nglobal = nglobal + nir(1,1:R);
end
occ = sum(reshape(nglobal(auxall), size(auxall,1), T), 1);
if occ(t) > 0
    return % tag in use
end
if any(occ(1:t-1) == 0)
    return % a lower tag is free: that entry fires instead
end

% consume the parent job at the fork
newgl = glspace;
newgl{isf_f}(end-R+r) = newgl{isf_f}(end-R+r) - 1;

% emit fjentry.weight (tasksPerLink) siblings per branch: sequential
% application over the partial outcome list handles branches sharing the
% same head node, repeated emissions on the same branch, and expands
% phase-entry mixtures of non-exponential sibling services
partials = {newgl};
partprob = 1;
B = length(fjentry.branchheads);
emissions = repmat(1:B, 1, fjentry.weight);
for b=emissions
    bh = fjentry.branchheads(b);
    isf_b = sn.nodeToStateful(bh);
    a = fjentry.auxclasses(b);
    newpartials = {};
    newpartprob = [];
    for pp=1:length(partials)
        curgl = partials{pp};
        [arvspace, ~, arvprob] = State.afterEvent(sn, bh, curgl{isf_b}, EventType.ARV, a, isSimulation, []);
        if isempty(arvspace)
            % sibling arrival blocked (cannot occur under the per-tag
            % auxiliary class capacity invariant); disable the firing
            outGlobalStates = {};
            outrate = [];
            outprob = [];
            return
        end
        for io=1:size(arvspace,1)
            nextgl = curgl;
            nextgl{isf_b} = arvspace(io,:);
            newpartials{end+1,1} = nextgl; %#ok<AGROW>
            if length(arvprob) >= io
                newpartprob(end+1,1) = partprob(pp)*arvprob(io); %#ok<AGROW>
            else
                newpartprob(end+1,1) = partprob(pp); %#ok<AGROW>
            end
        end
    end
    partials = newpartials;
    partprob = newpartprob;
end

outGlobalStates = partials;
outprob = partprob(:) * fjentry.prob;
outrate = GlobalConstants.Immediate * ones(length(partials),1);

if isSimulation
    if length(partials) > 1
        % sample one outcome (afterEvent already samples per-node mixtures
        % when isSimulation is true, so this is defensive)
        cum_prob = cumsum(outprob) / sum(outprob);
        firing_ctr = 1 + max([0,find( rand > cum_prob' )]);
        outGlobalStates = outGlobalStates(firing_ctr);
        outrate = outrate(firing_ctr);
    end
    % the phase-entry choice has already been sampled inside afterEvent,
    % so the firing competes at the full immediate rate
    outprob = 1;
end

end
