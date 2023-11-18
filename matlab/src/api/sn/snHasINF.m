function bool = snHasINF(sn)
% BOOL = HASINF()

bool = any(sn.schedid==SchedStrategy.ID_INF);
end