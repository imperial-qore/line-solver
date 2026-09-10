function sn = getStruct(self, wantInitialState)
% QN = GETSTRUCT(WANTINITSTATE)

if ~self.hasStruct
    refreshStruct(self);
    % REFUSE A CLASS ROUTED TO A STATION THAT CANNOT SERVE IT. Runs HERE, not at
    % the end of refreshStruct, because only here is self.sn final: measured
    % inside refreshStruct the guard saw nodevisits that did not yet carry the
    % flow it must test, and passed a model it refuses one frame later. Gated on
    % ~hasStruct so it costs one pass per BUILD, not one per getStruct call.
    checkServiceReachable(self);
end

if nargin == 1 || wantInitialState
    [self.sn.state, self.sn.stateprior, self.sn.space] = self.getState;
end

sn = self.sn;

end
