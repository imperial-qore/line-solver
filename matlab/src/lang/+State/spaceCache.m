function SS = spaceCache(n, m)
% SS = SPACECACHE(N, M)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

SS = [];
try
    HS = nchoosek(1:sum(n),sum(m)); % hotset
    for i=1:size(HS,1)
        SS = [SS; perms(HS(i,:))];
    end
catch ME
    line_error(mfilename,'The model is too large, LINE cannot generate the state space explicitly.');
end

end
