function bool = snHasOpenClasses(sn)
% BOOL = HASOPENCLASSES()

bool = any(isinf(sn.njobs));
end