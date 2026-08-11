function [QN,UN,RN,TN,CN,XN,runtime] = solver_custom_analyzer(sn, options)
% [QN,UN,RN,TN,CN,XN,RUNTIME] = SOLVER_CUSTOM_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Tstart = tic;

line_printf('Any activity prior or after launching the solution algorithm needs to be implemented in %s.m\n.',mfilename);
[QN,UN,RN,TN,CN,XN] = solver_custom(sn, options);

QN(isnan(QN))=0;
CN(isnan(CN))=0;
RN(isnan(RN))=0;
UN(isnan(UN))=0;
XN(isnan(XN))=0;
TN(isnan(TN))=0;

runtime = toc(Tstart);
end