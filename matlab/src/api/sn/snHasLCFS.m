function bool = snHasLCFS(sn)
% BOOL = HASLCFS()

bool = any(sn.schedid==SchedStrategy.ID_LCFS);
end