function N = getNumberOfJobs(self)
% N = GETNUMBEROFJOBS()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

K = getNumberOfClasses(self);
N = zeros(K,1); % changed later
classes = self.classes;
for k=1:K
    switch classes{k}.type
        case JobClassType.CLOSED
            N(k) = classes{k}.population;
        case JobClassType.OPEN
            N(k) = Inf;
        case JobClassType.DISABLED
            N(k) = 0;
    end
end
end
