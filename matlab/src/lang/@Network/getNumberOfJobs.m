function N = getNumberOfJobs(self)
% N = GETNUMBEROFJOBS()

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

K = getNumberOfClasses(self);
N = zeros(K,1); % changed later
classes = self.classes;
for k=1:K
    switch classes{k}.type
        case 'closed'
            N(k) = classes{k}.population;
        case 'open'
            N(k) = Inf;
    end
end
end
