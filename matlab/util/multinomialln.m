function r = multinomialln(m)
% r = MULTINOMIALLN(M)
% Logarithm of multinomial coefficient sum(M)!/(M(1)!* M(2)! * ... * M(n)!)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
r = gammaln(1+sum(m));
    for i=1:length(m)
        r = r - gammaln(1+m(i));
    end
end