function value = softmin(x,y,alpha)
% VALUE = SOFTMIN(X,Y,ALPHA)
%
% Smooth approximation of min(x,y). The literal weighted-average form
%   (x*exp(-alpha*x) + y*exp(-alpha*y)) / (exp(-alpha*x) + exp(-alpha*y))
% overflows once alpha*|x| leaves the double range and then evaluates
% Inf*0 = NaN. Writing lo = min(x,y), gap = |x-y| and w = exp(-alpha*gap)
% gives the algebraically identical
%   (lo + hi*w) / (1 + w) = lo + gap*w/(1 + w),
% whose exponent argument is never positive, so w stays in (0,1] and the
% limit w -> 0 returns min(x,y) exactly.
lo = min(x,y);
hi = max(x,y);
gap = hi - lo;
% exp(-t) underflows to exactly 0 for t > 745.13; beyond that softmin is min.
if ~(gap < 745 / alpha)
    value = lo;
    return
end
w = exp(-alpha*gap);
value = lo + gap*w/(1 + w);
end
