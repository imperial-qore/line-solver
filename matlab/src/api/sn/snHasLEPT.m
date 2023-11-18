function bool = snHasLEPT(sn)
% BOOL = HASLEPT()

bool = any(sn.schedid==SchedStrategy.ID_LEPT);
end