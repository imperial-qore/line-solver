function bool = snHasFCFS(sn)
% BOOL = HASFCFS()

bool = any(sn.schedid==SchedStrategy.ID_FCFS);
end