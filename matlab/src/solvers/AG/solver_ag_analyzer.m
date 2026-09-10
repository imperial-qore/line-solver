function [QN,UN,RN,TN,CN,XN,runtime,method,totiter,percResults] = solver_ag_analyzer(sn, options)
% [QN,UN,RN,TN,CN,XN,RUNTIME,METHOD,TOTITER,PERCRESULTS] = SOLVER_AG_ANALYZER(SN, OPTIONS)
%
% Analyzer for the agent-based (RCAT) solver. Every (station, class) pair is an
% isolated component and the components are coupled only through the reversed
% rates of the synchronizing actions, so the whole analysis is the fixed point
% over that scalar vector; see solver_ag.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Tstart = tic;
percResults = [];

if ~isfield(options, 'config')
    options.config = struct();
end

% RCAT builds a CTMC per component out of (D0,D1), so a preserved Det would
% reach it with no matrix at all and be read back as its mean rate. Unlike the
% MAM analyzer, which decides this per method, every AG method needs the Det
% fitted, so the value is not negotiable and is only left alone when the caller
% has pinned it deliberately.
if ~isfield(options.config, 'preserveDet')
    options.config.preserveDet = false;
end
if ~isfield(options.config, 'phfit')
    % A concentrated matrix exponential is not a generator, so a component
    % assembled from it is a rational generator whose stationary solution is a
    % signed vector. RCAT needs a genuine phase-type, as SSA, Fluid and JMT do.
    options.config.phfit = 'ph';
end

% The conversion RETAGS procid to APH/ME/MAP, so afterwards nothing names the law
% the user declared; only DET survives the retagging. Keep the declared tags for
% the gate in SolverAG.supportsModelMethod, which reads them by name.
procidDeclared = sn.procid;
sn = sn_nonmarkov_toph(sn, options);
sn.procidDeclared = procidDeclared;

method = options.method;
if isempty(method)
    method = 'default';
end

line_debug('AG analyzer starting: method=%s', method);

switch method
    case {'default', 'inap', 'inapplus', 'inapinf', 'exact'}
        % RCAT methods per Marin, Rota Bulo, Balsamo (INAP, MASCOTS 2012);
        % see _kb/06-solver-catalog.md for rationale
        [QN,UN,RN,TN,CN,XN,totiter] = solver_ag(sn, options);
    otherwise
        line_error(mfilename,'Unknown method: %s', method);
end

% 'default' resolves to inap, and 'exact' falls back to it with a warning
% (autocat is unreachable). Report what actually ran: 'exact' is registered
% globally as an exact method, so leaving the name in place would banner an
% iterative approximation as exact.
if any(strcmpi(method, {'default','exact'}))
    method = 'inap';
end

for i=1:sn.nstations
    switch sn.sched(i)
        case SchedStrategy.EXT
            TN(i,:) = sn.rates(i,:);
    end
end

QN(isnan(QN))=0;
CN(isnan(CN))=0;
RN(isnan(RN))=0;
UN(isnan(UN))=0;
XN(isnan(XN))=0;
TN(isnan(TN))=0;

runtime = toc(Tstart);
end
