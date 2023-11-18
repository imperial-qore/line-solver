function bool = snHasMultiClassHeterFCFS(sn)
% BOOL = SNHASMULTICLASSHETERFCFS()

iset = find(sn.schedid == SchedStrategy.ID_FCFS);
if isempty(iset)
    bool = false;
else
    bool = false;
    for i=iset(:)'
        bool = bool | range([sn.rates(i,:)])>0;
    end
end
end
