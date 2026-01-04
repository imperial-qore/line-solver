function v = multichoosecon(n,S)
% v = MULTICHOOSECON(n,S)
% Pick vectors of S elements from the available units in vector n
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if S == 1
    v = [];
    for i=find(n)
        v(end+1,1:size(n,2)) = 0;
        v(end,i) = 1;
    end
    return
end

v =[];
for i = find(n) % for all
    n_1 = n; n_1(i) = n_1(i) - 1; 
    T = multichoosecon(n_1,S-1);
    y = zeros(size(T,1),size(n,2)); y(:,i) = 1;
    v = [v; y+T];
end
%v= sortrows(unique(v,'rows'));
end