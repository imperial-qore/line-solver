function bool = snHasClassSwitching(sn)
% BOOL = HASCLASSWITCHING()

bool = sn.nclasses ~= sn.nchains;
end