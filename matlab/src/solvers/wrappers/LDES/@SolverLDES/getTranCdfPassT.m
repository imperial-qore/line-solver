function RD = getTranCdfPassT(self, R)
% RD = GETTRANCDFPASST(R)
%
% Empirical passage time CDF: for LDES, passage times are approximated by
% response times, so this delegates to getTranCdfRespT -- the same measured
% ecdf the engine's per-job samples define, as in the JAR and the C++ CLI.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    RD = self.getTranCdfRespT();
else
    RD = self.getTranCdfRespT(R);
end
end
