function [AvgTable,QT,UT,RT,WT,AT,TT] = getNodeAvgT(self,Q,U,R,T,A,W,keepDisabled)
% [AVGTABLE,QT,UT,RT,WT,AT,TT] = GETNODEAVGT(SELF,Q,U,R,T,A,W,KEEPDISABLED)
% Alias for getAvgNodeTable
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin<8
    keepDisabled = false;
end

if nargin == 1
    [AvgTable,QT,UT,RT,WT,AT,TT] = self.getAvgNodeTable();
elseif nargin == 2
    [AvgTable,QT,UT,RT,WT,AT,TT] = self.getAvgNodeTable(Q);
elseif nargin == 7
    [AvgTable,QT,UT,RT,WT,AT,TT] = self.getAvgNodeTable(Q,U,R,T,A,W);
else
    [AvgTable,QT,UT,RT,WT,AT,TT] = self.getAvgNodeTable(Q,U,R,T,A,W,keepDisabled);
end
end
