function bool = snHasMultiClassFCFS(sn)
% BOOL = SNHASMULTICLASSFCFS()

iset = find(sn.schedid == SchedStrategy.ID_FCFS);
bool = false;
for i=iset(:)'
    bool = bool | sum(sn.rates(i,:)>0)>1;
end
end