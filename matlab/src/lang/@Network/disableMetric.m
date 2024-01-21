function self = disableMetric(self, Y)
% SELF = DISABLEMETRIC(Y)

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

if iscell(Y)
    Y={Y{:}};
    for i=1:length(Y)
        Y{i}.disable=true;
    end
else
    Y.disable=true;
end
end
