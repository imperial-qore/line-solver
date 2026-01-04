function v = multichoose(n,k)
% v = MULTICHOOSE(n,k)
% Chooses k elements out of n with repetition
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

v=[];
if n==1
    v=k;
    return
elseif k==0
    v=zeros(1,n);
else
    last=0;
    for i=0:k
        w=multichoose(n-1,k-i);
        for j=1:size(w,1)
            v(end+1,:)=[i w(j,:)];
        end %for
    end %for
end %if
end
