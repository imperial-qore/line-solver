function out = ncDispatch(self, sn, options)
% OUT = NCDISPATCH(SN, OPTIONS)
%
% One inner solve of the normalizing-constant analyzer, in the contract that
% @NetworkSolver/fjFixedPoint.m expects (see @SolverMVA/mvaDispatch.m). It is
% used on the fork-join path only: the MMT transformation returns a plain
% mixed queueing network, so the specialised NC routes (order-independent,
% cache, loss network, finite capacity regions) cannot apply to it, and the
% load-dependent route is selected here exactly as runAnalyzer does.
%
% The auxiliary open classes that carry the parallelism start at an arrival
% rate of GlobalConstants.FineTol. On an open transformed model the analyzer
% resolves them with the exact open-network formulas, so no method override is
% needed; on a closed or mixed one the normalizing constant has to separate a
% chain whose throughput is ~1e-8, which is where the accuracy of this path
% has to be judged (see the fork-join tests).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

subopts = options;
if ~isempty(sn.lldscaling) || ~isempty(sn.cdscaling) || ~isempty(sn.jdscaling)
    [QN,UN,RN,TN,CN,XN,lG,runtime,lastiter,actualmethod] = solver_ncld_analyzer(sn, subopts);
else
    [QN,UN,RN,TN,CN,XN,lG,runtime,lastiter,actualmethod] = solver_nc_analyzer(sn, subopts);
end

out = struct('QN', QN, 'UN', UN, 'RN', RN, 'TN', TN, 'CN', CN, 'XN', XN, ...
    'lG', lG, 'runtime', runtime, 'lastiter', lastiter, 'method', options.method, ...
    'actualmethod', '');
if exist('actualmethod','var') && ~isempty(actualmethod)
    out.actualmethod = actualmethod;
end
end
