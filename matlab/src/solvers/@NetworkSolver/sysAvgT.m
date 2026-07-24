function [AvgSysChainTable, CT,XT] = sysAvgT(self,R,T)
% [AVGSYSCHAINTABLE, CT,XT] = SYSAVGT(SELF,R,T)
% Alias for getAvgSysTable
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin == 1
    [AvgSysChainTable, CT, XT] = self.getAvgSysTable();
elseif nargin == 2
    [AvgSysChainTable, CT, XT] = self.getAvgSysTable(R);
else
    [AvgSysChainTable, CT, XT] = self.getAvgSysTable(R,T);
end
end
