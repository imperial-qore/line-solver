function [QNclass_t, UNclass_t, TNclass_t] = getTranAvg(self,Qt,Ut,Tt)
% [QNCLASS_T, UNCLASS_T, TNCLASS_T] = GETTRANAVG(SELF,QT,UT,TT)

% Return transient average station metrics
%
% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

% temporarily switch to closing method
if nargin == 1
    [Qt,Ut,Tt] = self.getTranHandles;
end

options = self.options;
switch options.method 
    case 'default'
        self.options.method = 'closing';
    otherwise
        line_warning(mfilename,'getTranAvg is not offered by the specified method. Setting the solution method to ''''closing''''.\n');        
        self.options.method = 'closing';
        self.resetResults;
end

[QNclass_t, UNclass_t, TNclass_t] = getTranAvg@NetworkSolver(self,Qt,Ut,Tt);
self.options = options;
end
