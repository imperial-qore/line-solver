function bool = snHasSJF(sn)
% BOOL = HASSJF()

bool = any(sn.schedid==SchedStrategy.ID_SJF);
end