function bool = snHasHOL(sn)
% BOOL = HASHOL()

bool = any(sn.schedid==SchedStrategy.ID_HOL);
end