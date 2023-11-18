function bool = snHasLCFSPR(sn)
% BOOL = HASLCFSPR()

bool = any(sn.schedid==SchedStrategy.ID_LCFSPR);
end