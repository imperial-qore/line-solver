function bool = snHasHomogeneousScheduling(sn, strategy)
% BOOL = HASHOMOGENEOUSSCHEDULING(STRATEGY)

bool = length(findstring(sn.sched,strategy)) == sn.nstations;
end