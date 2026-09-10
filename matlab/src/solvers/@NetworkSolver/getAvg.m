function [QNclass,UNclass,RNclass,TNclass,ANclass,WNclass] = getAvg(self,Q,U,R,T,A,W)
% [QNCLASS,UNCLASS,RNCLASS,TNCLASS,ANCLASS,WNCLASS] = GETAVG(SELF,Q,U,R,T,A,W)
%
% Compute steady-state average metrics (queue length, utilization, response time,
% throughput, arrival rate, residence time) for all stations and job classes.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

self.assertNotChainModel('getAvg');

if isa(self.model, 'LayeredNetwork')
    % Java-backed LQN simulation (SolverLDES): per-LQN-element vectors,
    % same convention as SolverLN/SolverLQNS getAvg
    self.obj.getAvg(); % runs the LN LDES analyzer
    avgTable = self.obj.getLNAvgTable();
    [QNclass,UNclass,RNclass,WNclass,ANclass,TNclass] = JLINE.arrayListToResults(avgTable);
    return
end

sn = self.model.getStruct();

if strcmp(self.options.lang,'java') && ~strcmp(self.name,'SolverLDES')
    T0=tic;
    M = sn.nstations;
    R = sn.nclasses;
    % Force fresh Java model conversion each time to ensure current state
    % Only re-create the Java model if the model is not already Java-native
    if ~self.model.isJavaNative()
        self.model.obj = [];
        self.setLang();
    end
    % Carry this solver's verbosity into the JAR object rather than forcing STD:
    % a layer solve asked to stay quiet (SolverLN, options.verbose=false) would
    % otherwise print one progress line per layer per iteration.
    switch self.options.verbose
        case VerboseLevel.SILENT
            self.obj.getOptions.verbose = jline.VerboseLevel.SILENT;
        case VerboseLevel.DEBUG
            self.obj.getOptions.verbose = jline.VerboseLevel.DEBUG;
        otherwise
            self.obj.getOptions.verbose = jline.VerboseLevel.STD;
    end
    % Carry the JLINE-side fork-join (MMT) iterate across the model rebuild so
    % layers warm-start as in lang='matlab'. see _kb/05-solvers-overview.md
    isMVA = isa(self,'SolverMVA');
    if isMVA && self.options.config.fj_warmstart && ~isempty(self.fjForkLambda)
        self.obj.setForkWarmStart(JLINE.from_line_matrix(self.fjForkLambda));
    end
    SolverResult = self.obj.getAvg();
    if isMVA && self.model.hasFork
        self.fjForkLambda = JLINE.from_jline_matrix(self.obj.getForkWarmStart());
    end
    QN = JLINE.from_jline_matrix(SolverResult.QN);
    UN = JLINE.from_jline_matrix(SolverResult.UN);
    RN = JLINE.from_jline_matrix(SolverResult.RN);
    TN = JLINE.from_jline_matrix(SolverResult.TN);
    AN = JLINE.from_jline_matrix(SolverResult.AN);
    WN = JLINE.from_jline_matrix(SolverResult.WN);
    runtime=SolverResult.runtime;
    method=SolverResult.method;
    % Extract cache hit/miss probabilities from Java model
    for ind = 1:sn.nnodes
        if sn.nodetype(ind) == NodeType.Cache
            jnode = self.model.obj.getNodeByIndex(ind-1);
            hitRatioVec = JLINE.from_jline_matrix(jnode.getHitRatio());
            missRatioVec = JLINE.from_jline_matrix(jnode.getMissRatio());
            % Store per-class hit/miss probabilities matching MATLAB format:
            % only parent classes with hitClass>0 have non-zero entries
            hitClass = self.model.nodes{ind}.getHitClass;
            nk = length(hitClass);
            hitprob = zeros(1, nk);
            missprob = zeros(1, nk);
            for k = 1:nk
                if hitClass(k) > 0 && k <= length(hitRatioVec)
                    hitprob(k) = hitRatioVec(k);
                    missprob(k) = missRatioVec(k);
                end
            end
            self.model.nodes{ind}.setResultHitProb(hitprob);
            self.model.nodes{ind}.setResultMissProb(missprob);
            % Retrieval-cache extras (delayed-hit ratio, per-list hit ratio and
            % expected latency) copied raw from the Java cache so getAvgCacheTable
            % matches the native path. Empty for plain caches.
            self.model.nodes{ind}.setResultDelayedHitProb(JLINE.from_jline_matrix(jnode.getDelayedHitRatio()));
            self.model.nodes{ind}.setResultHitProbList(JLINE.from_jline_matrix(jnode.getHitRatioByList()));
            self.model.nodes{ind}.setResultItemProb(JLINE.from_jline_matrix(jnode.getItemProb()));
            self.model.nodes{ind}.setResultResidT(JLINE.from_jline_matrix(jnode.getResidT()));
        end
    end
    if any(sn.nodetype == NodeType.Cache)
        self.model.refreshStruct(true);
    end
    self.setAvgResults(QN,UN,RN,TN,AN,WN,[],[],runtime,method,1);
    % Finite Capacity Region (FCR) metrics appended as pseudo-stations
    % M+1..M+F (JMT only). see _kb/06-solver-catalog.md for rationale
    if sn.nregions > 0 && strcmp(self.name,'SolverJMT')
        jnode = self.obj.getAvgNode();
        Qn = JLINE.from_jline_matrix(jnode.QN);
        Un = JLINE.from_jline_matrix(jnode.UN);
        Rn = JLINE.from_jline_matrix(jnode.RN);
        Wn = JLINE.from_jline_matrix(jnode.WN);
        Tn = JLINE.from_jline_matrix(jnode.TN);
        An = JLINE.from_jline_matrix(jnode.AN);
        for f = 1:sn.nregions
            self.result.Avg.Q(M+f,:) = Qn(sn.nnodes+f,:);
            self.result.Avg.U(M+f,:) = Un(sn.nnodes+f,:);
            self.result.Avg.R(M+f,:) = Rn(sn.nnodes+f,:);
            self.result.Avg.W(M+f,:) = Wn(sn.nnodes+f,:);
            self.result.Avg.T(M+f,:) = Tn(sn.nnodes+f,:);
            self.result.Avg.A(M+f,:) = An(sn.nnodes+f,:);
        end
    end
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
    % A non-renewal (MAP/MMPP/MMAP) process that the resolved method cannot
    % consume is solved through its random-environment image rather than
    % rejected; every other unsupported feature keeps its rejection in
    % runAnalyzer. see _kb/05-solvers-overview.md for rationale
    if self.needsMapEnv(self.options)
        self.mapEnvApprox(self.options);
    else
        runAnalyzer(self);
    end
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
    maskprobe('Q', QNclass, filterMetric(Q, self.result.Avg.Q, [], sn, K, M), class(self), RNclass);
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
    maskprobe('U', UNclass, filterMetric(U, self.result.Avg.U, [], sn, K, M), class(self), RNclass);
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
    % Finite-server open station unstable when offered load rho >= 1: Util capped
    % at 1.0, QLen/RespT set Inf. see _kb/06-solver-catalog.md for rationale
    anyUnstable = false;
    if any(isinf(sn.njobs)) && ~isempty(TNclass)
        srcStations = sn.nodeToStation(sn.nodetype==NodeType.Source);
        openCls = isinf(sn.njobs(:).');
        for i = 1:size(UNclass,1)
            c = sn.nservers(i);
            if ~isfinite(c) || c <= 0 || any(srcStations==i)
                continue % infinite-server / delay / source station
            end
            % A STATION WITH A BINDING BUFFER CANNOT BE UNSTABLE, however heavily
            % it is offered: the buffer bounds the queue and the excess is blocked
            % or lost. rho = T/(c*rate) reaches exactly 1 at such a station -- that
            % is what a saturated server WITH a finite buffer looks like -- so
            % without this the queue length of an M/M/1/K in overload is reported
            % as Inf while the solver's own answer is the cap. An exact solver
            % escapes it only because its carried throughput lands a hair below the
            % service rate; the fluid answer sits exactly on it.
            bounded = false;
            if isfield(sn,'cap') && numel(sn.cap) >= i && isfinite(sn.cap(i))
                bounded = true;
            end
            if ~bounded && isfield(sn,'classcap') && size(sn.classcap,1) >= i
                bounded = any(isfinite(sn.classcap(i,:)));
            end
            if ~bounded && isfield(sn,'regionmembers') && ~isempty(sn.regionmembers)
                for f = 1:numel(sn.regionmembers)
                    if ~isempty(sn.regionmembers{f}) && numel(sn.regionmembers{f}) >= i ...
                            && sn.regionmembers{f}(i)
                        bounded = true;
                        break
                    end
                end
            end
            if bounded
                continue
            end
            % A STATION CUSTOMERS ABANDON CANNOT BE UNSTABLE EITHER, and for the
            % same reason: reneging bounds the queue however heavily it is
            % offered, the excess leaving instead of accumulating. rho =
            % T/(c*rate) reaches exactly 1 there -- that is what a saturated
            % server WITH abandonment looks like -- so without this the Erlang A
            % queue length that qsys_erlanga and qsys_ggisgi_fluid compute
            % exactly would be overwritten with Inf.
            if isfield(sn,'impatienceClass') && ~isempty(sn.impatienceClass) ...
                    && size(sn.impatienceClass,1) >= i ...
                    && any(sn.impatienceClass(i,:) == ImpatienceType.RENEGING)
                continue
            end
            % A load-dependent station serves faster than its nominal rate;
            % without the peak scaling the test reads the rate at population
            % one and calls a stable station saturated (e.g. the discrete-time
            % p(n)=p*min(n,s) server of Daduna's example 2.10).
            lldpeak = 1;
            if isfield(sn,'lldscaling') && ~isempty(sn.lldscaling) && size(sn.lldscaling,1) >= i
                row = sn.lldscaling(i,:);
                row = row(isfinite(row) & row > 0);
                if ~isempty(row)
                    lldpeak = max(row);
                end
            end
            rho = zeros(1,size(UNclass,2));
            for r = 1:size(UNclass,2)
                if sn.rates(i,r) > 0 && TNclass(i,r) > 0
                    rho(r) = TNclass(i,r) / (c * lldpeak * sn.rates(i,r));
                end
            end
            rhoOpen = sum(rho(openCls));
            rhoTot = sum(rho);
            if rhoOpen >= 1 && rhoTot > 0
                anyUnstable = true;
                UNclass(i,:) = rho / rhoTot; % station total capped to 1.0
                % Saturated station: QLen/RespT Inf for open classes.
                % see _kb/06-solver-catalog.md for rationale
                if ~isempty(RNclass)
                    RNclass(i,openCls) = Inf;
                end
                if ~isempty(QNclass)
                    QNclass(i,openCls) = Inf;
                end
            end
        end
    end
    if anyUnstable
        line_warning(mfilename,'The model has unstable queues (utilization >= 1); station utilization is reported capped at 1.0, queue length and response time as Inf.\n')
    end
end
end

function maskprobe(tag, withMask, noMask, solverName, RN)
% TEMPORARY MEASUREMENT ONLY (task #32 blast radius). Returns nothing and
% changes no result: it records where dropping the R-based zeroMask at
% getAvg.m:181/:193 would change the reported metric. Appends to a file rather
% than a global so an example calling `clear all` cannot wipe the record.
f = getenv('LINE_MASKPROBE_FILE');
if isempty(f), return; end
d = abs(noMask - withMask);
[ii, rr] = find(d > 0);
if isempty(ii), return; end
nm = getenv('LINE_MASKPROBE_NAME');
fid = fopen(f, 'a');
if fid < 0, return; end
for z = 1:numel(ii)
    fprintf(fid, '%s\t%s\t%s\t%d\t%d\t%.17g\t%.17g\t%.17g\n', nm, solverName, tag, ...
        ii(z), rr(z), withMask(ii(z), rr(z)), noMask(ii(z), rr(z)), RN(ii(z), rr(z)));
end
fclose(fid);
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

% Zero metrics for unreachable classes, but trust the simulation for
% fork-join, cache, spawn-fed and SPN models where the routing-based visit
% equations do not apply. see _kb/06-solver-catalog.md for rationale
hasForkJoin = any(sn.nodetype == NodeType.Fork) && any(sn.nodetype == NodeType.Join);
hasSPN = any(sn.nodetype == NodeType.Place) || any(sn.nodetype == NodeType.Transition);
hasCache = any(sn.nodetype == NodeType.Cache);

spawnFedChain = false(sn.nchains, 1);
if isfield(sn, 'classspawn') && ~isempty(sn.classspawn)
    for k2 = 1:min(K, length(sn.classspawn))
        sc2 = sn.classspawn(k2);
        if sc2 >= 0 && sc2 < K
            spawnFedChain(sn.chains(:, sc2+1) > 0) = true;
        end
    end
end

if sn.nchains > 0 && ~hasSPN  % Only check reachability if chains are defined and not SPN
    for k = 1:K
        c = sn.chains(:, k)>0;
        if any(c)  % Only check if class k belongs to a chain
            for i = 1:M
                if sn.visits{c}(i,k) == 0
                    % For fork-join models, don't zero out if the metric has a non-zero value
                    % from simulation - the visits calculation doesn't capture fork-join semantics
                    % where Join outputs the parent class
                    if (hasForkJoin || hasCache || any(spawnFedChain(c))) && ~isempty(metric) && metric(i,k) > GlobalConstants.FineTol
                        continue; % Trust the simulation result
                    end
                    outData(i,k) = 0;
                end
            end
        end
    end
end
end
