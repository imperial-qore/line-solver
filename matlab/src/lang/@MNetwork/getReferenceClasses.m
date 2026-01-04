function refclass = getReferenceClasses(self)
% REFSTAT = GETREFERENCECLASSES()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

K = getNumberOfClasses(self);
refclass = false(K,1);
for k=1:K
    refclass(k,1) =  self.classes{k}.isrefclass;
end
end
