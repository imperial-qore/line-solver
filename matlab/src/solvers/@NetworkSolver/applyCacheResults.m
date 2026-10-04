function applied = applyCacheResults(self, cacheRec)
% APPLIED = APPLYCACHERESULTS(CACHEREC)
%
% Copy a cache surface produced on ANOTHER model onto this solver's own model,
% matched by node NAME, and refresh the struct so the hit/miss routing follows.
%
% Why this exists. A Cache does not report in Q/U/R/T/A/W/C/X; it reports on the
% Cache NODE, and getAvgCacheTable reads it back off self.model.nodes{ind}. A
% driver that solves a TRANSFORMED image of the model (mapEnvApprox, which hands
% SolverENV a map2renv environment whose stage 1 is a throwaway deep copy) leaves
% this model's nodes untouched, so the cache surface is lost and, worse, the
% hit/miss class throughputs silently fall back to the nodevisits split -- the
% 0.5/0.5 guess link() laid down before any cache was analysed.
%
% MATCHED BY NAME, never by index: the image and the caller agree on node order
% today, but that is buildStage's invariant rather than a law, and a mismatch
% must be an error rather than a cache result written onto a Sink.
%
% ABSENT MUST CLEAR, NOT KEEP. A quantity the coupling did not produce is
% written as empty, so a previous solver's value cannot survive into this one's
% table and leave the probabilities failing to sum to 1.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

applied = false;
if isempty(cacheRec)
    return
end

sn = self.model.getStruct();
byName = containers.Map('KeyType','char','ValueType','any');
for ind = 1:sn.nnodes
    if sn.nodetype(ind) == NodeType.Cache
        byName(self.model.nodes{ind}.getName()) = ind;
    end
end

for j = 1:numel(cacheRec)
    rec = cacheRec(j);
    if ~isKey(byName, rec.name)
        line_error(mfilename, sprintf(['The environment image reports a cache ' ...
            'named ''%s'', which this model does not have. The image and the ' ...
            'model must carry the same Cache nodes.'], rec.name));
    end
    node = self.model.nodes{byName(rec.name)};
    node.setResultHitProb(asRow(rec.hitprob));
    node.setResultMissProb(asRow(rec.missprob));
    node.setResultDelayedHitProb(asRow(rec.delayedhitprob));
    if isfield(rec, 'hitproblist')
        node.setResultHitProbList(rec.hitproblist);
    end
    applied = applied || ~isempty(rec.hitprob);
end

if applied
    % The visits that carry the hit and miss classes are DERIVED from the
    % split, so the struct must be rebuilt before any node-level throughput
    % is reconstructed. This is what a native cache solve does.
    self.model.refreshStruct(true);
end
end

function v = asRow(v)
if isempty(v)
    v = [];
    return
end
v = full(v(:)).';
end
