function [pi,kmax]=ctmc_uniformization(pi0,Q,t,tol,maxiter)
% [PI,KMAX]=CTMC_UNIFORMIZATION(PI0,Q,T,TOL,MAXITER)
%
% MAXITER caps the Poisson series truncation depth; pass a nonpositive
% value (or omit it) to size it adaptively as max(100, q*t+10*sqrt(q*t)+20).

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if nargin<4%~exist('tol','var')
    tol = 1e-12;
end
if nargin<5%~exist('maxiter','var')
    maxiter = -1;
end
q=1.1*max(abs(diag(Q)));
% Split the horizon so exp(-q*t) never underflows within a segment
% (exp(-745)==0 in double precision): exp(Q*t)=(exp(Q*t/nSeg))^nSeg
MAXQT = 500;
if q*t > MAXQT
    nSeg = ceil(q*t/MAXQT);
    tSeg = t/nSeg;
    pi = pi0;
    for seg=1:nSeg
        [pi,kmax] = ctmc_uniformization(pi,Q,tSeg,tol,maxiter);
    end
    return
end
if maxiter<=0
    % The Poisson(q*t) mass concentrates around q*t with spread O(sqrt(q*t));
    % a fixed cap silently truncates the series for large horizons
    maxiter = max(100, ceil(q*t+10*sqrt(q*t)+20));
end
Qs=speye(size(Q))+sparse(Q)/q;
k=0;
s=1;
r=1;
iter=0;
kmax=1;
while iter<maxiter
    iter=iter+1;
    k=k+1;
    r=r*(q*t)/k;
    s=s+r;
    if (1-exp(-q*t)*s)<=tol
        kmax=k;
        break;
    end
    % Best-effort truncation depth if the loop exhausts maxiter: summing k
    % terms is always more accurate than the single term a stale kmax=1 gives
    kmax=k;
end

pi=pi0*(exp(-q*t));
P=pi0;
ri=exp(-q*t);
for j=1:kmax
    P=P*Qs;
    ri=ri*(q*t/j);
    pi=pi+ri*P;
end
end
