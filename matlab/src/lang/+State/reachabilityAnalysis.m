function [SS] = reachabilityAnalysis(sn, cutoff, options)
% [SS] = RECHABILITYANALSYSIS(QN, CUTOFF)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% LINE state space generator limited to states reachable from the initial
% state.
%
% SS: state space
% SSh: hashed state space
% sn: updated sn

%%
if ~exist('cutoff','var') && any(isinf(Np)) % if the model has open classes
    line_error(mfilename,'Unspecified cutoff for open classes in state space generator.');
end

if numel(cutoff)==1
    cutoff = cutoff * ones(sn.nstations, sn.nclasses);
end
T0=tic;
if nargin<2 %~exist('options','var')
    options = self.getOptions;
end

if self.enableChecks && ~SolverSSA.supports(self.model)
    line_error(mfilename,'This model contains features not supported by the solver.');
end

self.runAnalyzerChecks(options);

switch options.lang
    case 'java'
    case 'matlab'
        sn = getStruct(self);
        [StateSpaceAggr,arvRates,depRates] = solver_ssa(sn, options);       
        SS = sn.space;
end
end
