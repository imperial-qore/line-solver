function bool = snHasPriorities(sn)
% BOOL = SNHASPRIORITIES()
bool = any(sn.classprio(:) > 0);
end