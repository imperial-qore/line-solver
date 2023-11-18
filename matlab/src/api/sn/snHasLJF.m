function bool = snHasLJF(sn)
% BOOL = HASLJF()

bool = any(sn.schedid==SchedStrategy.ID_LJF);
end