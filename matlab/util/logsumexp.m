function lse = logsumexp(x)
% L = LOGSUMEXP(X)
% Approximate the logarithm of a sum of exponentials
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%xstar = max(x);
%lse = xstar + log(sum(exp(x-xstar)));

% Implementation described in:
% P. Blanchard et al., ACCURATE COMPUTATION OF THE LOG-SUM-EXP AND SOFTMAX
% FUNCTIONS, IMA Journal of Numerical Analysis, draa038, https://doi.org/10.1093/imanum/draa038
n = length(x);
[a, k] = max(x);
w = zeros(1,n);
s = 0;
for i=1:n
    w(i) =  exp(x(i)-a);
    if i~=k
        s = s + w(i);
    end
end
lse = a + log1p(s);
end
