function [QNclass_t, UNclass_t, TNclass_t] = getTranAvg(self, Qt, Ut, Tt)
% [QNCLASS_T, UNCLASS_T, TNCLASS_T] = GETTRANAVG(SELF, QT, UT, TT)
%
% Returns transient mean performance metrics over time using
% QBD-based Taylor series analysis (libQBD adaptive).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.


% lang='cpp' takes the transient means from line-cli (-s mam -a tran), which
% runs the same ldqbd Taylor-series analysis this getter forces below, and is
% held to the same single-class precondition.
%
% The tables are stored where the native path stores its own and the ANSWER IS
% BUILT BY @NetworkSolver/getTranAvg, exactly as in the line below: that is what
% wraps a [value, t] table into a metricVal struct and puts NaN on a station the
% ldqbd analysis carries no curve for, which here is the Source.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    [Qcpp, Ucpp, Tcpp] = CPPLINE.tranAvg(self.name, self.model, self.options);
    if nargin == 1
        [Qt, Ut, Tt] = self.getTranHandles;
    end
    empt = cell(size(Qcpp));
    self.setTranAvgResults(Qcpp, Ucpp, empt, Tcpp, empt, empt, 0);
    [QNclass_t, UNclass_t, TNclass_t] = getTranAvg@NetworkSolver(self, Qt, Ut, Tt);
    return
end

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
