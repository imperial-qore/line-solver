function sn = getStruct(self, wantInitialState)
% QN = GETSTRUCT(WANTINITSTATE)

if ~self.hasStruct
    refreshStruct(self);
end

if nargin == 1 || wantInitialState
    [self.sn.state, self.sn.stateprior] = self.getState;
end

sn = self.sn;
end
