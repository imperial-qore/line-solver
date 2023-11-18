function bool = snHasMultiServer(sn)
% BOOL = SNHASMULTISERVER()

bool = any(sn.nservers(isfinite(sn.nservers)) > 1);
end