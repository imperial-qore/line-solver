function [QN,UN,RN,TN,CN,XN,runtime,actualMethod] = solver_qns_analyzer(sn, options)
% [QN,UN,RN,TN,CN,XN,RUNTIME,ACTUALMETHOD] = SOLVER_QNS_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Tstart = tic;

line_debug('LQNS analyzer starting: method=%s', options.method);

[QN,UN,RN,TN,CN,XN,actualMethod] = solver_qns(sn, options);

runtime = toc(Tstart);
end