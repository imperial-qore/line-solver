function bool = snHasDPS(sn)
% BOOL = HASDPS()

bool = any(sn.schedid==SchedStrategy.ID_DPS);
end