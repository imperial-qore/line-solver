function bool = snHasMixedClasses(sn)
% BOOL = HASMIXEDCLASSES()

bool = snHasClosedClasses(sn) && snHasOpenClasses(sn);
end
    
