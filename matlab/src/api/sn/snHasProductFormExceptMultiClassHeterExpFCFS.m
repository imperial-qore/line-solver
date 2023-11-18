function bool = snHasProductFormExceptMultiClassHeterExpFCFS(sn)
% BOOL = HASPRODUCTFORMEXCEPTMULTICLASSHETEREXPFCFS()

bool = all(sn.schedid==SchedStrategy.ID_INF | sn.schedid==SchedStrategy.ID_PS | sn.schedid==SchedStrategy.ID_FCFS | sn.schedid==SchedStrategy.ID_LCFSPR | sn.schedid==SchedStrategy.ID_EXT);
bool = bool & ~snHasPriorities(sn);
bool = bool & ~snHasForkJoin(sn);
iset = find(sn.schedid == SchedStrategy.ID_FCFS);
for i=iset(:)'
    icset = isfinite(sn.scv(i,:)) & sn.scv(i,:)>0;
    bool = bool & all(sn.scv(i,icset) > 1-GlobalConstants.FineTol) & all(sn.scv(i,icset) < 1+GlobalConstants.FineTol);
end
end