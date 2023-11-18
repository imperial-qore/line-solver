function bool = snIsClosedModel(sn)
bool = all(isfinite(sn.njobs));
end