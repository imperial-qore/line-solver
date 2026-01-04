function inspace1 = cartesian(inspace1, inspace2)
% INSPACE1 = CARTESIAN(INSPACE1, INSPACE2)
%
% Cartesian product of two matrices. It replicates elements of the first 
% input matrix and pairs them with each row of the second input matrix.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin<2 && iscell(inspace1)
    C = inspace1;
    inspace1 = [];
    for c=1:length(C)
        inspace1 = State.cartesian(inspace1, C{c});
    end
    return
end

if isempty(inspace1)
    inspace1 = inspace2;
    return
end
if isempty(inspace2)
    return
end

n1 = size(inspace1,1); m1 = size(inspace1,2);
n2 = size(inspace2,1); m2 = size(inspace2,2);
inspace1 = repmat(inspace1, n2, 1);

curStates = 1:n1;
for s=1:n2
    inspace1(curStates,(m1+1):(m1+m2)) = repmat(inspace2(s,:),length(curStates),1);
    curStates = (curStates) + n1;
end
end
