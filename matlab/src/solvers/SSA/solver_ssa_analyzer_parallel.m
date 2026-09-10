function [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,snc,StartN,PreemptN]=solver_ssa_analyzer_parallel(sn, init_state, laboptions)
% [XN,UN,QN,RN,TN,CN]=SOLVER_SSA_ANALYZER_PARALLEL(SN, INIT_STATE, LABOPTIONS)
%
% Worker-count-invariant parallel SSA.
%
% The simulation budget is split into a FIXED number of independent
% replications R (laboptions.config.nreplicas). Replication r simulates
% ceil(samples/R) events and is seeded deterministically with (seed+r-1):
% inside a parfor body spmdIndex resolves to 1 in solver_ssa, so the
% effective seed is exactly laboptions.seed+r-1, independent of which
% worker executes the iteration. The R replications are distributed over
% whatever workers the parallel pool provides (and run on the client if no
% pool exists), and the per-replication estimates are averaged.
%
% Because replication r always uses the same seed and sample budget no
% matter how many workers are available, the returned averages depend only
% on (seed, samples, R) and are invariant to the worker count. This is the
% key difference from the previous spmd implementation, whose result varied
% with the pool size (it both divided the budget by, and seeded labs from,
% the runtime number of labs).

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

snc = sn;

% Fixed number of independent replications (worker-count invariant). The
% default lives in SolverOptions ('SSA' case); guard here so the analyzer
% remains usable if called with a hand-built options struct.
if isfield(laboptions,'config') && isfield(laboptions.config,'nreplicas') ...
        && ~isempty(laboptions.config.nreplicas)
    R = max(1, round(laboptions.config.nreplicas));
else
    R = 8;
end

baseSeed = laboptions.seed;
perRepSamples = ceil(laboptions.samples / R);
line_debug('SSA parallel: %d replications x %d events each (requested budget %d, effective %d)', ...
    R, perRepSamples, laboptions.samples, perRepSamples*R);

% per-replication outputs (sliced so they can be populated inside parfor)
XNr = cell(1,R);
UNr = cell(1,R);
QNr = cell(1,R);
RNr = cell(1,R);
TNr = cell(1,R);
CNr = cell(1,R);
StartNr = cell(1,R);
PreemptNr = cell(1,R);
sncr = cell(1,R);

parfor r = 1:R
    repoptions = laboptions;
    repoptions.samples = perRepSamples;
    repoptions.verbose = VerboseLevel.SILENT;
    % Deterministic per-replication seed (see header). solver_ssa adds
    % (lab_idx-1) which is 0 inside parfor, so this seed is used verbatim.
    repoptions.seed = baseSeed + r - 1;
    [XNr{r},UNr{r},QNr{r},RNr{r},TNr{r},CNr{r},sncr{r},StartNr{r},PreemptNr{r}] = ...
        run_replica(snc, init_state, repoptions);
end

% average the per-replication estimates
QN = cellsum(QNr)/R;
UN = cellsum(UNr)/R;
RN = cellsum(RNr)/R;
TN = cellsum(TNr)/R;
CN = cellsum(CNr)/R;
XN = cellsum(XNr)/R;
% the derived rates average across replications like every other estimate
StartN = cellsum(StartNr)/R;
PreemptN = cellsum(PreemptNr)/R;

% average cache actual hit/miss probabilities across replications
for k=1:snc.nclasses
    for isf=1:sn.nstateful
        if sn.nodetype(isf) == NodeType.Cache
            ind = sn.statefulToNode(isf);
            sn.nodeparam{ind}.actualhitprob(k) = 0;
            sn.nodeparam{ind}.actualmissprob(k) = 0;
            for l=1:R
                qntmp = sncr{l};
                if length(qntmp.nodeparam{ind}.hitclass)>=k
                    sn.nodeparam{ind}.actualhitprob(k) = sn.nodeparam{ind}.actualhitprob(k) + (1/R) * qntmp.nodeparam{ind}.actualhitprob(k);
                    sn.nodeparam{ind}.actualmissprob(k) = sn.nodeparam{ind}.actualmissprob(k) + (1/R) * qntmp.nodeparam{ind}.actualmissprob(k);
                end
            end
        end
    end
end
snc = sn;
tranSysState=[];
tranSync=[];
end

function [XN,UN,QN,RN,TN,CN,sncl,StartN,PreemptN]=run_replica(sn, init_state, repoptions)
% Run a single SSA replication and reduce it to per-station/class averages.
% Kept as a local function so the (parfor-unfriendly) nested indexing runs
% in an ordinary function workspace.

M = sn.nstations;
K = sn.nclasses;
PH = sn.proc;
S = sn.nservers;
NK = sn.njobs';
rates = sn.rates;

% eventCache is not shared across replications
if isfield(repoptions,'config') && isfield(repoptions.config,'eventcache')
    eventCache = EventCache.create(repoptions.config.eventcache, sn);
else
    eventCache = EventCache.create(false, sn);
end

% see _kb/06-solver-catalog.md for rationale (SSA utilization estimator)
userCap = sn.cap;
userClasscap = sn.classcap;
[probSysState,SSq,arvRates,depRates,~,~,sncl,~,startRates,preemptRates] = solver_ssa(sn, init_state, repoptions, eventCache);

XN = NaN*zeros(1,K);
UN = NaN*zeros(M,K);
QN = NaN*zeros(M,K);
RN = NaN*zeros(M,K);
TN = NaN*zeros(M,K);
CN = NaN*zeros(1,K);
StartN = zeros(M,K);
PreemptN = zeros(M,K);
for k=1:K
    refsf = sncl.stationToStateful(sncl.refstat(k));
    XN(k) = probSysState*depRates(:,refsf,k);
    for ist=1:M
        isf = sncl.stationToStateful(ist);
        TN(ist,k) = probSysState*depRates(:,isf,k);
        QN(ist,k) = probSysState*SSq(:,(ist-1)*K+k);
        StartN(ist,k) = probSysState*startRates(:,isf,k);
        PreemptN(ist,k) = probSysState*preemptRates(:,isf,k);
        switch sncl.sched(ist)
            case SchedStrategy.INF
                UN(ist,k) = QN(ist,k);
            otherwise
                % see _kb/06-solver-catalog.md for rationale (SSA utilization estimator)
                if ~isempty(PH{ist}{k})
                    % A LOAD-DEPENDENT STATION IS NORMALIZED BY ITS PEAK
                    % CAPACITY, max(c, max(alpha)), not by the server count: the
                    % scaling multiplies the nominal rate, so dividing by c
                    % alone reports the work delivered against a capacity the
                    % station has already exceeded, and gives a utilization
                    % ABOVE ONE (measured 1.6529 on a closed Delay+Queue, N=4,
                    % alpha = [1 1.5 2 2.5], where the answer is 0.6612). Same
                    % ceff as solver_ssa_analyzer_serial.m and SolverCTMC.
                    ceff = S(ist);
                    if ~isempty(sncl.lldscaling) && ist <= size(sncl.lldscaling,1)
                        ceff = max(ceff, max(sncl.lldscaling(ist,:)));
                    end
                    % see _kb/06-solver-catalog.md for rationale (SSA utilization estimator)
                    if isinf(sncl.njobs(k)) && (isfinite(userCap(ist)) || isfinite(userClasscap(ist,k)))
                        UN(ist,k) = TN(ist,k)/rates(ist,k)/ceff;
                    else
                        UN(ist,k) = probSysState*arvRates(:,ist,k)/rates(ist,k)/ceff;
                    end
                end
        end
    end
end

for k=1:K
    for ist=1:M
        if TN(ist,k)>0
            RN(ist,k) = QN(ist,k)./TN(ist,k);
        else
            RN(ist,k) = 0;
        end
    end
    CN(k) = NK(k)./XN(k);
end

QN(isnan(QN))=0;
CN(isnan(CN))=0;
RN(isnan(RN))=0;
UN(isnan(UN))=0;
XN(isnan(XN))=0;
TN(isnan(TN))=0;

% update cache actual hit and miss data for this replication
TNcache = zeros(sncl.nstateful,K);
for k=1:K
    for isf=1:sncl.nstateful
        if sncl.nodetype(isf) == NodeType.Cache
            TNcache(isf,k) = probSysState*depRates(:,isf,k);
        end
    end
end
for k=1:K
    for isf=1:sncl.nstateful
        if sncl.nodetype(isf) == NodeType.Cache
            ind = sncl.statefulToNode(isf);
            if length(sncl.nodeparam{ind}.hitclass)>=k
                h = sncl.nodeparam{ind}.hitclass(k);
                m = sncl.nodeparam{ind}.missclass(k);
                sncl.nodeparam{ind}.actualhitprob(k) = TNcache(isf,h)/sum(TNcache(isf,[h,m]));
                sncl.nodeparam{ind}.actualmissprob(k) = TNcache(isf,m)/sum(TNcache(isf,[h,m]));
            end
        end
    end
end
end
