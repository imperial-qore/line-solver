function resetStruct(self)
self.sn = [];
self.hasStruct = false;
% Discarding the struct is a structural change like any other: bump the
% version so that caches keyed on it (see MNetwork.structVersion) are
% invalidated even if the next compilation reproduces the same contents.
self.structVersion = self.structVersion + 1;
end
