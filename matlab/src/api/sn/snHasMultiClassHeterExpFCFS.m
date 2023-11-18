function bool = snHasMultiClassHeterExpFCFS(sn)
% BOOL = SNHASMULTICLASSHETEREXPFCFS()

iset = find(sn.schedid == SchedStrategy.ID_FCFS);
if isempty(iset)
    bool = false;
else
    bool = true;
    for i=iset(:)'
        bool = bool & range([sn.rates(i,:)])>0 & all(sn.scv(i,:) > 1-GlobalConstants.FineTol) & all(sn.scv(i,:) < 1+GlobalConstants.FineTol);
    end
end
end