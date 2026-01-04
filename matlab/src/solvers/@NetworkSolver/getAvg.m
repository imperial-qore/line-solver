function [QNclass,UNclass,RNclass,TNclass,ANclass,WNclass] = getAvg(self,Q,U,R,T,A,W)
% [QNCLASS,UNCLASS,RNCLASS,TNCLASS,ANCLASS,WNCLASS] = GETAVG(SELF,Q,U,R,T,A,W)
%
% Compute steady-state average metrics (queue length, utilization, response time,
% throughput, arrival rate, residence time) for all stations and job classes.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.model.getStruct();

if strcmp(self.options.lang,'java') && ~strcmp(self.name,'SolverDES')
    T0=tic;
    M = sn.nstations;
    R = sn.nclasses;
    % object already created by NetworkSolver.setLang()
    self.obj.model.reset();
    self.obj.getOptions.verbose = jline.VerboseLevel.STD;
    SolverResult = self.obj.getAvg();
    QN = JLINE.from_jline_matrix(SolverResult.QN);
    UN = JLINE.from_jline_matrix(SolverResult.UN);
    RN = JLINE.from_jline_matrix(SolverResult.RN);
    TN = JLINE.from_jline_matrix(SolverResult.TN);
    AN = JLINE.from_jline_matrix(SolverResult.AN);
    WN = JLINE.from_jline_matrix(SolverResult.WN);
    runtime=SolverResult.runtime;
    method=SolverResult.method;
    self.setAvgResults(QN,UN,RN,TN,AN,WN,[],[],runtime,method,1);
    QNclass = reshape(QN,M,R);
    UNclass = reshape(UN,M,R);
    RNclass = reshape(RN,M,R);
    TNclass = reshape(TN,M,R);
    ANclass = reshape(AN,M,R);
    WNclass = reshape(WN,M,R);
    return
end

%%
if nargin == 1 % no parameter
    if   isempty(self.model.handles) || ~isfield(self.model.handles,'Q') || ...
            ~isfield(self.model.handles,'U') || ~isfield(self.model.handles,'R') || ...
            ~isfield(self.model.handles,'T') || ~isfield(self.model.handles,'A') || ...
            ~isfield(self.model.handles,'W')
        reset(self); % reset results in case there are partial results saved
    end
    [Q,U,R,T,A,W] = self.getAvgHandles;
elseif nargin == 2
    handlers = Q;
    [Q,U,R,T,A,W] = deal(handlers{:}); % set Q=handlers{1}, U=handlers{2}, ...
end

if isfield(self.options,'timespan')
    if isfinite(self.options.timespan(2))
        line_error(mfilename,'The getAvg method does not support the timespan option, use the getTranAvg method instead.');
    end
else
    self.options.timespan = [0,Inf];
end

if ~self.hasAvgResults() || ~self.options.cache
    runAnalyzer(self);
    % the next line is required because getAvg can alter the chain
    % structure in the presence of caches so we need to reload sn
    sn = self.model.getStruct;
    if ~self.hasAvgResults
        line_error(mfilename,'Unable to return results for this model.');
    end
end % else return cached value


M = sn.nstations;
K = sn.nclasses;

% Check if this is an SPN model (Places don't have response times)
hasSPN = any(sn.nodetype == NodeType.Place) || any(sn.nodetype == NodeType.Transition);

if ~isempty(R)
    RNclass = filterMetric(R, self.result.Avg.R, [], sn, K, M);
else
    RNclass = [];
end

if ~isempty(Q)
    % For SPNs, don't zero Q based on R because Places don't have response times
    if hasSPN
        zeroMaskQ = [];
    else
        zeroMaskQ = RNclass < 10 * GlobalConstants.FineTol;
    end
    QNclass = filterMetric(Q, self.result.Avg.Q, zeroMaskQ, sn, K, M);
else
    QNclass = [];
end

if ~isempty(U)
    % For SPNs, don't zero U based on R because Places don't have response times
    if hasSPN
        zeroMaskU = [];
    else
        zeroMaskU = RNclass < 10 * GlobalConstants.FineTol;
    end
    UNclass = filterMetric(U, self.result.Avg.U, zeroMaskU, sn, K, M);
else
    UNclass = [];
end

if ~isempty(T)
    TNclass = filterMetric(T, self.result.Avg.T, [], sn, K, M);
else
    TNclass = [];
end

if ~isempty(A)
    zeroMask = false(size(RNclass));
    zeroMask(sn.nodeToStation(sn.nodetype==NodeType.Source),:) = true;
    ANclass = filterMetric(A, self.result.Avg.A, zeroMask, sn, K, M);
else
    ANclass = [];
end

if ~isempty(W)
    WNclass = sn_get_residt_from_respt(sn, RNclass, W);
else
    WNclass = [];
end

if ~isempty(UNclass)
    unstableQueues = find(sum(UNclass,2)>0.99 * sn.nservers);
    if any(unstableQueues) && any(isinf(sn.njobs))
        line_warning(mfilename,'The model has unstable queues, performance metrics may grow unbounded.\n')
    end
end
end

function outData = filterMetric(handle, metric, zeroMask, sn, K, M)
% post-process Avg measure
outData = zeros(M, K);
for k = 1:K
    for i = 1:M
        if ~handle{i,k}.disabled && ~isempty(metric)
            outData(i,k) = metric(i,k);
        else
            outData(i,k) = NaN;
        end
    end
end

% NaN values indicate that a metric is disabled
outData(isnan(outData)) = 0;
% set to zero entries associated to immediate transitions
outData(zeroMask) = 0;
% round to zero numerical perturbations
outData(outData < GlobalConstants.FineTol) = 0;

% set to zero metrics for classes that are unreachable
% but skip this check for fork-join models where the visits calculation
% doesn't correctly capture the parent class visiting Join and downstream nodes
% also skip when no chains are defined (routing not specified)
% also skip for SPN models where Places don't have traditional visits
hasForkJoin = any(sn.nodetype == NodeType.Fork) && any(sn.nodetype == NodeType.Join);
hasSPN = any(sn.nodetype == NodeType.Place) || any(sn.nodetype == NodeType.Transition);

if sn.nchains > 0 && ~hasSPN  % Only check reachability if chains are defined and not SPN
    for k = 1:K
        c = sn.chains(:, k)>0;
        if any(c)  % Only check if class k belongs to a chain
            for i = 1:M
                if sn.visits{c}(i,k) == 0
                    % For fork-join models, don't zero out if the metric has a non-zero value
                    % from simulation - the visits calculation doesn't capture fork-join semantics
                    % where Join outputs the parent class
                    if hasForkJoin && metric(i,k) > GlobalConstants.FineTol
                        continue; % Trust the simulation result
                    end
                    outData(i,k) = 0;
                end
            end
        end
    end
end
end
