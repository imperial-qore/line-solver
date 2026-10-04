function resetStruct(self)
% The linked routing table and the reward definitions are model inputs parked in sn, not derived from it, so they
% survive the discard as in resetModel; losing rtorig made SolverLN's relink() of a load-dependent layer index []
if isstruct(self.sn) && isfield(self.sn, 'rtorig') && ~isempty(self.sn.rtorig)
    kept.rtorig = self.sn.rtorig;
    if isfield(self.sn, 'reward')
        kept.reward = self.sn.reward;
    end
    self.sn = kept;
else
    self.sn = [];
end
self.hasStruct = false;
% Discarding the struct is a structural change like any other: bump the
% version so that caches keyed on it (see MNetwork.structVersion) are
% invalidated even if the next compilation reproduces the same contents.
self.structVersion = self.structVersion + 1;
end
