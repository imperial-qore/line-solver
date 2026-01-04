function [AvgChain,QTc,UTc,RTc,WTc,ATc,TTc] = chainAvgT(self,Q,U,R,T)
% [AVGCHAIN,QTC,UTC,RTC,WTC,ATC,TTC] = CHAINAVGT(SELF,Q,U,R,T)
% Alias for getAvgChainTable
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin == 1
    [AvgChain,QTc,UTc,RTc,WTc,ATc,TTc] = self.getAvgChainTable();
elseif nargin == 2
    [AvgChain,QTc,UTc,RTc,WTc,ATc,TTc] = self.getAvgChainTable(Q);
else
    [AvgChain,QTc,UTc,RTc,WTc,ATc,TTc] = self.getAvgChainTable(Q,U,R,T);
end
end
