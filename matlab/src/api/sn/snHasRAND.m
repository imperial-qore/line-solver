function bool = snHasRAND(sn)
% BOOL = HASRAND()

bool = any(sn.schedid==SchedStrategy.ID_SIRO);
end