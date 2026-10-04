function finalizeAvgResults(self, QN, UN, RN, TN, CN, XN, runtime, method, iter, actualmethod, lG, WN)
% FINALIZEAVGRESULTS(QN,UN,RN,TN,CN,XN,RUNTIME,METHOD,ITER,ACTUALMETHOD,LG,WN)
%
% Common closing of every runAnalyzer: derive the arrival rates from the
% throughputs, store the averages under the reported method name, record the
% normalizing constant and flag a run that exhausted the wall-clock budget.
%
% ACTUALMETHOD is the analyzer-selected method appended to 'default/' when the
% caller asked for 'default'; pass '' when the analyzer did not report one. LG
% and WN are optional: pass [] when the solver has no normalizing constant or
% does not report residence times (only MVA does today).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 11, actualmethod = ''; end
if nargin < 12, lG = []; end
if nargin < 13, WN = []; end

sn = self.model.getStruct();
AN = sn_get_arvr_from_tput(sn, TN, self.getAvgTputHandles());

if strcmp(method,'default') && ~isempty(actualmethod) && ~strcmp(actualmethod,'default')
    reported = ['default/' actualmethod];
else
    reported = method;
end
self.setAvgResults(QN,UN,RN,TN,AN,WN,CN,XN,runtime,reported,iter);

if ~isempty(lG)
    self.result.Prob.logNormConstAggr = real(lG);
end

options = self.getOptions;
if lineTimeoutExceeded(options)
    self.result.Avg.timedOut = true;
    line_warning(mfilename,'Solver exceeded the wall-clock time budget (options.timeout=%gs); returning the interim solution.\n', options.timeout);
end
end
