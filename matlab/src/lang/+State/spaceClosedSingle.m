function SS = spaceClosedSingle(M, N, caps)
% SS = SPACECLOSEDSINGLE(M, N)
% SS = SPACECLOSEDSINGLE(M, N, CAPS)
%
% All the ways N identical jobs of ONE class can sit at M nodes. CAPS, when
% given, is a 1xM per-node bound on the class: the recursion is then pruned AT
% THE BRANCH rather than after it, so a class that can only occupy one node out
% of M enumerates M rows instead of nchoosek(N+M-1,M-1).

% Copyright (c) 2012-2014, Imperial College London
% All rights reserved.
if M==0
    SS = [];
elseif nargin < 3 || isempty(caps)
    SS = multichoose(M, N);
else
    SS = multichoose_capped(M, N, caps(:)');
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

function [v] = multichoose_capped(n,k,caps)
% [V] = MULTICHOOSE_CAPPED(N,K,CAPS) - MULTICHOOSE with a per-bin bound.
%
% The bound is applied BEFORE recursing, which is the whole point: a bin of
% capacity 0 takes only the zero item, so a class-switching chain whose classes
% each visit one station enumerates its own station rather than every station.
% The unbounded form generated the full lattice and left the caller to reject
% the impossible rows one at a time, at one State.fromMarginal call each.

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
v=[];
if k > sum(caps(1:n))
    return
end
if n==1
    if k <= caps(1)
        v=k;
    end
    return
elseif k==0
    v=zeros(1,n);
else
    hi = min(k, caps(1));
    for i=0:hi
        w=multichoose_capped(n-1,k-i,caps(2:end));
        for j=1:size(w,1)
            v(end+1,:)=[i w(j,:)]; %#ok<AGROW>
        end %for
    end %for
end %if
end
