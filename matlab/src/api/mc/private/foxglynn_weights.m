function [left,right,w] = foxglynn_weights(lambda,tol,maxiter,normalize)
% [LEFT,RIGHT,W] = FOXGLYNN_WEIGHTS(LAMBDA,TOL,MAXITER)
% [LEFT,RIGHT,W] = FOXGLYNN_WEIGHTS(LAMBDA,TOL,MAXITER,NORMALIZE)
%
% Fox-Glynn truncation window and Poisson weights for a Poisson(LAMBDA)
% mixing distribution at tail-mass tolerance TOL. MAXITER caps the right
% truncation point; pass a nonpositive value to leave it uncapped.
%
% NORMALIZE defaults to true, which rescales the window to sum to one, as
% Fox and Glynn prescribe: the truncated tails are then redistributed over
% the window and the mixture is a proper average. Pass false to obtain the
% true Poisson probabilities instead, so that 1-sum(W) is the discarded tail
% rather than being absorbed; CTMC_FAU needs that, its error being reported
% as missing mass rather than as a bound.
%
% Private to matlab/src/api/mc: shared by CTMC_FOXGLYNN, which mixes over the
% Poisson process subordinating a uniformized chain, and by CTMC_FAU, which
% mixes over the same process subordinating the bidiagonal birth chain that
% carries its adaptive uniformization rates.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if nargin<2
    tol = 1e-12;
end
if nargin<3
    maxiter = -1;
end
if nargin<4
    normalize = true;
end
if tol<=0
    tol = 1e-12;
end
left = foxglynn_left(lambda,tol);
right = foxglynn_right(lambda,tol);
if maxiter>0 && right>maxiter
    right = maxiter;
    left = min(left,right);
end
w = foxglynn_poisson(lambda,left,right,normalize);
end

function e = foxglynn_chernoff(lambda,k)
% Chernoff exponent lambda*h(k/lambda) with h(u)=u*log(u)-u+1, so that exp(-e)
% dominates P{X>=k} for k>lambda and P{X<=k} for k<lambda under
% X~Poisson(lambda). It certifies the Fox-Glynn estimates below, which are
% asymptotic and valid only for lambda>=25.
if k<=0
    e = lambda;
else
    e = lambda - k + k*log(k/lambda);
end
end

function r = foxglynn_right(lambda,tol)
% Right truncation point R with P{X>R}<=tol/2. The starting guess is the
% Fox-Glynn (1988) right tail estimate: with m=floor(lambda) and shift
% s(k)=k*sqrt(2*lambda)+3/2, the tail P{X>=m+ceil(s(k))} is bounded by
% a*d*exp(-k^2/2)/(k*sqrt(2*pi)) with a=(1+1/lambda)*exp(1/16)*sqrt(2).
% Jansen's correction (2011), applied here as the factor
% d=1/(1-exp(-(2/9)*s(k))), repairs the original statement, which drops this
% factor for the k-dependent shift and is therefore optimistic at moderate
% lambda; d tends to one as lambda grows, so the corrected bound agrees with
% Fox and Glynn's asymptotically. The guess is then tightened and, if needed,
% grown until the Chernoff bound holds, so R is certified in any regime.
target = log(2/tol);
m = floor(lambda);
r = m;
if lambda >= 25
    a = (1+1/lambda)*exp(1/16)*sqrt(2);
    spread = sqrt(2*lambda);
    for k=1:64
        shift = k*spread + 1.5;
        d = 1/(1-exp(-(2/9)*shift));
        bound = a*d*exp(-0.5*k^2)/(k*sqrt(2*pi));
        if bound <= 0.5*tol
            r = m + ceil(shift);
            break
        end
    end
end
while r>m && foxglynn_chernoff(lambda,r)>=target
    r = r-1;
end
while foxglynn_chernoff(lambda,r+1)<target
    r = r+1;
end
end

function l = foxglynn_left(lambda,tol)
% Left truncation point L with P{X<L}<=tol/2, zero when no truncation is
% admissible. The starting guess is the Fox-Glynn (1988) left tail estimate
% with b=(1+1/lambda)*exp(1/(8*lambda)) and shift k*sqrt(lambda)+3/2 below the
% mode. Unlike the right tail this one needs no Jansen factor, the left tail of
% a Poisson being lighter than its normal approximation.
target = log(2/tol);
m = floor(lambda);
if foxglynn_chernoff(lambda,0) < target
    l = 0;
    return
end
l = 0;
if lambda >= 25
    b = (1+1/lambda)*exp(1/(8*lambda));
    spread = sqrt(lambda);
    for k=1:64
        bound = b*exp(-0.5*k^2)/(k*sqrt(2*pi));
        if bound <= 0.5*tol
            l = m - floor(k*spread + 1.5);
            break
        end
    end
    l = max(l,0);
end
while l>0 && foxglynn_chernoff(lambda,l-1) < target
    l = l-1;
end
while l<m && foxglynn_chernoff(lambda,l) >= target
    l = l+1;
end
end

function w = foxglynn_poisson(lambda,left,right,normalize)
% Poisson(lambda) weights on [left,right]. Following
% Fox-Glynn they are built by the two-sided recursion w(k-1)=w(k)*k/lambda and
% w(k+1)=w(k)*lambda/(k+1) anchored at the mode, so neither exp(-lambda) nor
% lambda^k/k! is ever evaluated and the overflow and underflow that limit the
% direct series cannot occur. Anchoring at w(mode)=1 keeps the extreme weights
% near tol, far above the denormal threshold, making Fox and Glynn's rescaling
% of the mode weight unnecessary. The normalizing sum is accumulated in
% increasing order of magnitude. Without normalization the anchor is instead
% scaled by the true mode probability, evaluated once through a log-gamma so
% that the factorial is never formed, which turns the same recursion into the
% Poisson probabilities themselves.
len = right-left+1;
w = zeros(1,len);
m = min(max(floor(lambda),left),right);
w(m-left+1) = 1;
for k=m:-1:left+1
    w(k-1-left+1) = w(k-left+1)*k/lambda;
end
for k=m:right-1
    w(k+1-left+1) = w(k-left+1)*lambda/(k+1);
end
if normalize
    w = w/sum(sort(w));
else
    w = w*exp(-lambda + m*log(lambda) - gammaln(m+1));
end
end
