function [AvgTable,QT,UT,RT,WT,AT,TT] = aT(self,Q,U,R,T,A,W,keepDisabled)
% [AVGTABLE,QT,UT,RT,WT,AT,TT] = AT(SELF,Q,U,R,T,A,W,KEEPDISABLED)
% Alias for getAvgTable
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin<8
    keepDisabled = false;
end

if nargin == 1
    [AvgTable,QT,UT,RT,WT,AT,TT] = self.getAvgTable();
elseif nargin == 2
    [AvgTable,QT,UT,RT,WT,AT,TT] = self.getAvgTable(Q);
elseif nargin == 7
    [AvgTable,QT,UT,RT,WT,AT,TT] = self.getAvgTable(Q,U,R,T,A,W);
else
    [AvgTable,QT,UT,RT,WT,AT,TT] = self.getAvgTable(Q,U,R,T,A,W,keepDisabled);
end
end
