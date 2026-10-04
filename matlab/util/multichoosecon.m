function v = multichoosecon(n,S)
% v = MULTICHOOSECON(n,S)
% Pick vectors of S elements from the available units in vector n
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
v = [];
coder.varsize('v');
if S == 1
    rowCount = 0;
    for i = 1:size(n,2)
        if n(i) > 0
            rowCount = rowCount + 1;
            v(rowCount,1:size(n,2)) = 0; %#ok<AGROW>
            v(rowCount,i) = 1;
        end
    end
    return
end

for i = 1:size(n,2)
    if n(i) > 0
        n_1 = n;
        n_1(i) = n_1(i) - 1;
        T = multichoosecon(n_1,S-1);
        y = zeros(size(T,1),size(n,2));
        y(:,i) = 1;
        v = [v; y+T]; %#ok<AGROW>
    end
end
%v= sortrows(unique(v,'rows'));
end
