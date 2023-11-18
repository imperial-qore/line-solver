function bool = snHasPS(sn)
% BOOL = HASPS()

bool = any(sn.schedid==SchedStrategy.ID_PS);
end