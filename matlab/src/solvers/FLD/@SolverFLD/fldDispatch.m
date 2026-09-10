function out = fldDispatch(self, sn, options)
% OUT = FLDDISPATCH(SN, OPTIONS)
%
% One inner solve of the fluid analyzer, in the contract that
% @NetworkSolver/fjFixedPoint.m expects (see @SolverMVA/mvaDispatch.m and
% @SolverNC/ncDispatch.m). It is used on the fork-join path only.
%
% The MMT transformation returns a plain mixed network built from Source,
% Delay, Queue, Router and ClassSwitch, every one of which the fluid drift
% already carries, so nothing here is fork-join specific: the callback solves
% the transformed struct exactly as runAnalyzer solves an ordinary one. The
% ODE has no product-form premise, which is why the transformed model, whose
% auxiliary open classes start at an arrival rate of GlobalConstants.FineTol,
% needs no approximate-mixed override of the kind ncDispatch documents.
%
% Only the steady-state means are returned. The transient tables (Qt, Ut, Tt)
% are indexed by the ORIGINAL stations and classes, whereas each pass of the
% fixed point integrates a different transformed network, so a trajectory read
% off the last pass would not be the trajectory of the model the caller built.
% getTranAvg therefore stays refused on a fork-join model.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0 = tic;
subopts = options;

% The transformed struct is rebuilt by the driver on every pass, so the
% non-Markovian conversion runAnalyzer applies to the original struct has to
% be re-applied here: the fluid drift reads mu*phi as a flow and a matrix
% exponential has neither.
subopts.config.phfit = 'ph';
sn = sn_nonmarkov_toph(sn, subopts);

[QN, UN, RN, TN, CN, XN, ~, ~, ~, ~, lastSol, lastiter] = solver_fluid_analyzer(sn, subopts);

if isempty(lastSol)
    M = sn.nstations; K = sn.nclasses;
    QN = NaN*ones(M,K); UN = NaN*ones(M,K); RN = NaN*ones(M,K); TN = NaN*ones(M,K);
    CN = NaN*ones(1,K); XN = NaN*ones(1,K);
end

out = struct('QN', QN, 'UN', UN, 'RN', RN, 'TN', TN, 'CN', CN, 'XN', XN, ...
    'lG', NaN, 'runtime', toc(T0), 'lastiter', lastiter, ...
    'method', subopts.method, 'actualmethod', subopts.method);
end
