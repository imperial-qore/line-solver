function bool = snHasClosedClasses(sn)
% BOOL = HASCLOSEDCLASSES()

bool = any(isfinite(sn.njobs));
end