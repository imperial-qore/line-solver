function [QNclass_t, UNclass_t, TNclass_t] = getTranAvg(self, Qt, Ut, Tt)
% [QNCLASS_T, UNCLASS_T, TNCLASS_T] = GETTRANAVG(SELF, QT, UT, TT)
%
% Returns transient mean performance metrics over time using
% QBD-based Taylor series analysis (libQBD adaptive).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin == 1
    [Qt, Ut, Tt] = self.getTranHandles;
end

options = self.options;
% Force MATLAB path for transient analysis (Java path does not support transient)
self.options.lang = 'matlab';
self.options.method = 'ldqbd';

[QNclass_t, UNclass_t, TNclass_t] = getTranAvg@NetworkSolver(self, Qt, Ut, Tt);
self.options = options;
end
