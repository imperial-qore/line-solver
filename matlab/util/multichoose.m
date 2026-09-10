function v = multichoose(n,k)
% v = MULTICHOOSE(n,k)
% Chooses k elements out of n with repetition
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

v=[];
coder.varsize('v');
% Cooperative wall-clock budget checkpoint (session-level deadline set by the
% per-solver runAnalyzer): combinatorial enumeration can explode in a single
% top-level call, so the budget must be enforced inside the recursion.
persistent tochk
if isempty(tochk)
    tochk = 0;
end
tochk = tochk + 1;
if tochk >= 4096
    tochk = 0;
    if lineTimeoutExceeded()
        line_error(mfilename,'Enumeration exceeded the wall-clock time budget (options.timeout).');
    end
end
if n==1
    v=k;
    return
elseif k==0
    v=zeros(1,n);
else
    last=0; %#ok<NASGU>
    for i=0:k
        w=multichoose(n-1,k-i);
        for j=1:size(w,1)
            v=[v; [i w(j,:)]]; %#ok<AGROW>
        end %for
    end %for
end %if
end
