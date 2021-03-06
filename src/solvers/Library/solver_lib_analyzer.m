function [QN,UN,RN,TN,CN,XN,runtime] = solver_lib_analyzer(sn, options)
% [QN,UN,RN,TN,CN,XN,RUNTIME] = SOLVER_LIB_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2021, Imperial College London
% All rights reserved.

Tstart = tic;

[QN,UN,RN,TN,CN,XN] = solver_lib(sn, options);

QN(isnan(QN))=0;
CN(isnan(CN))=0;
RN(isnan(RN))=0;
UN(isnan(UN))=0;
XN(isnan(XN))=0;
TN(isnan(TN))=0;

runtime = toc(Tstart);

if options.verbose > 0
    line_printf('\nSolver LIB analysis completed. Runtime: %f seconds.\n',runtime);
end
end
