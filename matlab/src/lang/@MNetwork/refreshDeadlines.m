function classdeadline = refreshDeadlines(self)
% CLASSDEADLINE = REFRESHDEADLINES()
% Extract deadline values from all job classes

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

K = getNumberOfClasses(self);
classdeadline = zeros(1,K);
for r=1:K
    classdeadline(r) = self.getClassByIndex(r).deadline;
end
if ~isempty(self.sn)
    self.sn.classdeadline = classdeadline;
end
end
