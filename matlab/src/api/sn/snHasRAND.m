function bool = snHasSIRO(sn)
% BOOL = SNHASSIRO()

bool = any(sn.schedid==SchedStrategy.ID_SIRO);
end
