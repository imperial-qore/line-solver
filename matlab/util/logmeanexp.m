function lme = logmeanexp(x)
% L = LOGMEANEXP(X)
% Approximate the logarithm of a mean of exponentials
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

    lme = logsumexp(x) - log(length(x));
end