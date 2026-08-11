function [AvgChain,QTc,UTc,RTc,WTc,ATc,TTc] = getNodeChainAvgT(self,Q,U,R,T)
% [AVGCHAIN,QTC,UTC,RTC,WTC,ATC,TTC] = GETNODECHAINAVGT(SELF,Q,U,R,T)
% Alias for getAvgNodeChainTable
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin == 1
    [AvgChain,QTc,UTc,RTc,WTc,ATc,TTc] = self.getAvgNodeChainTable();
elseif nargin == 2
    [AvgChain,QTc,UTc,RTc,WTc,ATc,TTc] = self.getAvgNodeChainTable(Q);
else
    [AvgChain,QTc,UTc,RTc,WTc,ATc,TTc] = self.getAvgNodeChainTable(Q,U,R,T);
end
end
