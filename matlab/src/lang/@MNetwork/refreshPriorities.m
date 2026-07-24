function classprio = refreshPriorities(self)
% CLASSPRIO = REFRESHPRIORITIES()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

K = getNumberOfClasses(self);
classprio = zeros(1,K);
for r=1:K
    classprio(r) = self.getClassByIndex(r).priority;
end
if ~isempty(self.sn)
    self.sn.classprio = classprio;
end
end
