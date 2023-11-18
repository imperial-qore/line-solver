function self = enableMetric(self, Y)
% SELF = ENABLEMETRIC(Y)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

if iscell(Y)
    Y={Y{:}};
    for i=1:length(Y)
        Y{i}.enabled = true;
    end
else
    Y.enabled = true;
end
end
