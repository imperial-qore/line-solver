function bool = snHasGPS(sn)
% BOOL = HASGPS()

bool = any(sn.schedid==SchedStrategy.ID_GPS);
end