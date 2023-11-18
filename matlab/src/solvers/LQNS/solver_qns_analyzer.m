function [QN,UN,RN,TN,CN,XN,runtime] = solver_qns_analyzer(sn, options)
% [QN,UN,RN,TN,CN,XN,RUNTIME] = SOLVER_QNS_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

Tstart = tic;

[QN,UN,RN,TN,CN,XN] = solver_qns(sn, options);

runtime = toc(Tstart);
end