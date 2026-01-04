function r=nchoosekln(n,m)
% r=NCHOOSEKLN(n,m)
% Logarithm of the binomial coefficient {n \choose m}
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
r=gammaln(1+n)-gammaln(1+(n-m))-gammaln(1+m);
end
