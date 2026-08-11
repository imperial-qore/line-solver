function SS = spaceClosedSingle(M, N)
% SS = SPACECLOSEDSINGLE(M, N)

% Copyright (c) 2012-2014, Imperial College London
% All rights reserved.
if M==0
    SS = [];
else
    SS = multichoose(M, N);
end
end

function [v] = multichoose(n,k)
% [V] = MULTICHOOSE(N,K)

v=[];
% Cooperative wall-clock budget checkpoint (session-level deadline set by the
% per-solver runAnalyzer); see matlab/util/multichoose.m for rationale.
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
    last=0;
    for i=0:k
        w=multichoose(n-1,k-i);
        for j=1:size(w,1)
            v(end+1,:)=[i w(j,:)];
        end %for
    end %for
end %if
end
