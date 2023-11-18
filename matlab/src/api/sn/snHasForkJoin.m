function bool = snHasForkJoin(sn)
% BOOL = SNHASFORKJOIN()
bool = any(sn.fj(:) > 0);
end