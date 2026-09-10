function reason = nc_method_refusal(sn, method, options, forReport)
% REASON = NC_METHOD_REFUSAL(SN, METHOD, OPTIONS, FORREPORT)
%
% May METHOD run on this model? '' when it may, otherwise the reason it may
% not, in the words the analyzer refuses with.
%
% FORREPORT (default TRUE) says WHICH QUESTION IS BEING ASKED, and for two
% method names the two questions have different answers:
%
%   true  -- "should MODEL.HELP offer this pair?" A pair that comes back as a
%            table of zeros must not be offered, so the answer is no.
%   false -- "what does the reference DO when asked for it by name?" For
%            'mmint2' and 'gleint' outside their shape the reference deliberately
%            WARNS AND RETURNS A ZERO TABLE (pfqn_nc.m, case {'mmint2','gleint'}:
%            lG = [] and return, unconditionally), and a caller who names the
%            method keeps that answer.
%
% THE ASYMMETRY IS A RULING, NOT AN OVERSIGHT (2026-07-25, reaffirmed when this
% gate was added): the report answers "should this be offered" and the run
% answers "what does the reference do". 'comomld' is NOT in that bucket --
% PFQN_COMOMRM_LD raises 'The solver accepts at most a single queueing station.'
% natively -- so it is refused on both paths.
%
% ONE PREDICATE, TWO CALLERS. @SolverNC/runAnalyzer asks it once, ahead of the
% dispatch, and turns a non-empty answer into an error; SOLVERNC.SUPPORTSMODELMETHOD
% asks it so that MODEL.HELP and MODEL.FINDSOLVER never offer a (solver,method)
% pair that would raise, and so that SOLVERAUTO never delegates to one. Two
% copies of these rules is precisely how the report and the run drift apart,
% which is the failure this function exists to prevent, so a new rule goes here
% and not at a call site.
%
% ONLY WHAT THE FEATURE REGISTRY CANNOT NAME LIVES HERE. A feature set says "I
% ACCEPT this construct", so it can refuse a model for HAVING something and
% never for LACKING it: "no think time" and "closed population only" are
% therefore expressed in SOLVERNC.GETMETHODFEATURESET by dropping
% SchedStrategy_INF and OpenClass, while "requires a cache", "requires
% state-dependent routing", "requires a loss network", "requires exactly two
% stations" and "requires normal usage" have no such form and are decided here.
% So is "THIS ROUTE READS NO METHOD NAME": the loss network, the M/M/1/K closed
% form, the paired LCFS convolution, the SDR product form and the cache
% analyzers each switch on a few method names and answer every other name with their
% own default, which would come back under the caller's label.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

reason = '';
if nargin < 3
    options = [];
end
if nargin < 4 || isempty(forReport)
    forReport = true;
end
if isempty(method)
    method = 'default';
end

% The discrete-time route answers for itself: NC_IS_DT_MODEL decides
% admissibility on the slot lattice, and every gate below is written about a
% continuous-time queueing network. What it does not decide is the NAME:
% SOLVER_NC_DT_ANALYZER evaluates DPFQN_NC / DQSYS_BERNOULLI1 whatever the
% method says, so only the two names that mean "the exact route" have one.
if isstruct(options) && isfield(options,'config') && isstruct(options.config) ...
        && isfield(options.config,'slotted') && ~isempty(options.config.slotted) ...
        && logical(options.config.slotted)
    if ~any(strcmpi(method,{'default','exact'}))
        reason = sprintf(['Method ''%s'' has no route on the discrete-time (slotted) product ' ...
            'form: solver_nc_dt_analyzer evaluates dpfqn_nc / dqsys_bernoulli1 whatever the ' ...
            'method says. Use ''default'' or ''exact''.'], method);
    end
    return
end

% -- discriminatory processor sharing ------------------------------------
% Morrison's heavy-usage expansion is the ONLY NC route that can see the DPS
% weights; every other method builds a product-form normalizing constant that
% silently drops them and answers with the egalitarian-PS network, which is a
% wrong number rather than a coarse one.
if nc_is_dps_model(sn)
    if ~any(strcmpi(method,{'default','morrison'}))
        reason = sprintf(['Method ''%s'' cannot represent the DPS weights of a ' ...
            'discriminatory processor-sharing station; it would return the egalitarian-PS ' ...
            'network. Use method ''default'' or ''morrison'' (npfqn_dps_morrison), SolverMVA, ' ...
            'SolverFLD or SolverCTMC.'], method);
    end
    return
elseif sn_has_dps(sn)
    % A DPS station outside Morrison's shape. SchedStrategy_DPS is declared in
    % the feature set because a boolean feature cannot express "this shape
    % only"; this is that imperative half.
    reason = ['SolverNC analyzes a discriminatory processor-sharing station ' ...
        'only in the shape Morrison''s expansion is derived for: a CLOSED network of exactly ' ...
        'two stations, one infinite-server (think) station and one single-server DPS station, ' ...
        'exponential service, each class visiting the two equally often. Use SolverMVA, ' ...
        'SolverFLD or SolverCTMC for any other DPS model.'];
    return
elseif strcmpi(method,'morrison')
    % The method named on a model that is not the shape at all -- not even a DPS
    % station in it. Left ungated it reaches no route of its own and falls
    % through to the ordinary normalizing-constant path, which would answer the
    % product-form model UNDER THE CALLER'S LABEL.
    reason = ['Method ''morrison'' is the heavy-usage expansion of a CLOSED ' ...
        'network of exactly two stations, one infinite-server (think) station and one ' ...
        'single-server DPS station with exponential service, which this model is not. ' ...
        'Remove the ''method'' option to let SolverNC choose, or use SolverMVA, SolverFLD ' ...
        'or SolverCTMC.'];
    return
end

% -- Krzesinski state-dependent routing ----------------------------------
% An SDR model is intercepted by SOLVER_NC_SDR_ANALYZER whatever the method
% says, so a name other than its own two would be answered by the product form
% of eq. (16) under that name; and the analyzer's own premises (closed, no
% class switching, every stateful node a station, one FCFS rate per chain)
% live in NC_SDR_REFUSAL, which it asks too. Reaching the second test means
% the model declares none.
if isfield(sn,'sdr') && ~isempty(sn.sdr)
    if ~any(strcmpi(method,{'default','sdr','sdr.mva'}))
        reason = sprintf(['Method ''%s'' has no route on a model with state-dependent routing, ' ...
            'which SOLVER_NC_SDR_ANALYZER intercepts whatever the method says. Use ''default'', ' ...
            '''sdr'' or ''sdr.mva''.'], method);
        return
    end
    reason = nc_sdr_refusal(sn);
    return
end
if any(strcmpi(method,{'sdr','sdr.mva'}))
    reason = sprintf('Method %s requires state-dependent routing, which this model does not declare.', method);
    return
end

% -- stochastic Petri net ------------------------------------------------
% A net is served only by the MDD-rec route, and none of the gates below --
% written about stations, capacities and the queueing-network product form --
% says anything about a net. SPN_PF decides its product-form class, by name.
if any(sn.nodetype == NodeType.Place)
    if ~any(strcmpi(method,{'default','rec'}))
        reason = sprintf(['a stochastic Petri net is solved by the MDD-rec route; ' ...
            'method ''%s'' is a normalizing-constant algorithm for queueing ' ...
            'networks. Use ''rec'' or ''default'''], method);
    end
    return
end

% -- order-independent and pass-and-swap stations --------------------------
% Every method other than the four listed reads sn.rates, which holds only the
% single-job rate mu([r]) of an OI station: the rank rate mu(n) is silently
% dropped and the answer is that of an ordinary queue. The same is true of
% EVERY method, 'default' included, once the model is outside the two shapes
% that have an OI route -- the closed BCMP network with OI stations of
% SOLVER_NC_OI_ANALYZER and the two-station P&S tandem of
% SOLVER_NC_PAS_IS_ANALYZER -- because the dispatch then falls through to
% SOLVER_NC. What each of the two analyzers refuses on its own shape is in
% NC_OI_REFUSAL, which they ask too.
hasOIorPAS = any(sn.sched == SchedStrategy.OI | sn.sched == SchedStrategy.PAS);
if nc_is_oi_model(sn)
    if ~any(strcmpi(method,{'default','exact','is','sampling'}))
        reason = sprintf(['Method ''%s'' cannot represent the rank rate mu(n) of an ' ...
            'order-independent station; use method ''default'' or ''exact'' (pfqn_ncoi), ' ...
            '''is'', SolverMVA, or SolverCTMC.'], method);
        return
    end
    if any(strcmpi(method,{'is','sampling'}))
        reason = nc_oi_refusal(sn, 'pas');
    else
        reason = nc_oi_refusal(sn, 'oi');
    end
    return
elseif nc_is_pas_model(sn)
    % a genuine swap graph: the ordered-state chain is reducible and only the
    % importance sampler carries its recurrent class (Comte & Dorsman 2021)
    if ~any(strcmpi(method,{'default','is','sampling'}))
        reason = sprintf(['Method ''%s'' cannot represent the swap dynamics of a pass-and-swap ' ...
            'station; the two-station P&S tandem is served by importance sampling alone ' ...
            '(pfqn_pas_is). Use method ''default'', ''is'' or ''sampling'', or SolverCTMC.'], method);
        return
    end
    reason = nc_oi_refusal(sn, 'pas');
    return
elseif hasOIorPAS
    reason = ['SolverNC serves an order-independent or pass-and-swap station only in a CLOSED ' ...
        'network whose other stations are BCMP product-form (exact OI convolution, ' ...
        'solver_nc_oi_analyzer) or in the two-station pass-and-swap tandem (importance sampling, ' ...
        'solver_nc_pas_is_analyzer); elsewhere the normalizing-constant routes read only the ' ...
        'single-job rate mu([r]) and would answer for an ordinary queue. Use SolverMVA or SolverCTMC.'];
    return
end

% -- the paired LCFS network ---------------------------------------------
% A non-preemptive LCFS station is served by SOLVER_NC_LCFSQN alone, on the
% closed two-station LCFS + LCFS-PR shape of NC_LCFS_REFUSAL (which SOLVER_NC
% asks too), and that convolution reads no method name.
if any(sn.sched == SchedStrategy.LCFS)
    reason = nc_lcfs_refusal(sn);
    if isempty(reason) && ~any(strcmpi(method,{'default','exact'}))
        reason = sprintf(['Method ''%s'' has no route on the paired LCFS / LCFS-PR network, which ' ...
            'SOLVER_NC_LCFSQN answers by its own convolution (pfqn_lcfsqn_ca) whatever the method ' ...
            'says. Use ''default'' or ''exact''.'], method);
    end
    return
end

% -- caches --------------------------------------------------------------
% 'rayint' and 'spm' both name the SPM saddle point of a cache (and, on a
% retrieval model, the ray/WKB delayed-hit expansion), so they are admissible
% here and nowhere else. Three cache routes, three name sets:
%   Source-Cache-Sink     SOLVER_NC_CACHE_ANALYZER switches on exact / sampling
%                         and sends EVERY other name to the SPM saddle point;
%   open retrieval cache  SOLVER_NC_RETRIEVAL_ANALYZER reads 'rayint' and
%                         answers every other name with the exact recurrences;
%   cache-queueing        DA_CACHEQN takes exact or SPM for the cache and hands
%                         the name to the network solve, where PFQN_NC has no
%                         arm for 'rayint'/'spm' and raises on a closed chain
%                         (an open one is answered before the method switch).
% CACHE_PROB_EREC, which 'exact' uses on the first and the third, is exact for
% the exchangeable (RR/FIFO) family only.
if any(sn.nodetype == NodeType.Cache)
    ci = find(sn.nodetype == NodeType.Cache, 1);
    hasRetrieval = isstruct(sn.nodeparam{ci}) && isfield(sn.nodeparam{ci},'retrievalSystemCapacity') ...
        && sn.nodeparam{ci}.retrievalSystemCapacity > 0;
    if nc_is_noreentrant_cache(sn)
        if ~any(strcmpi(method,{'default','exact','sampling','rayint','spm'}))
            reason = sprintf(['Method ''%s'' has no route on a Source-Cache-Sink model: ' ...
                'solver_nc_cache_analyzer would answer it with the SPM saddle point under that ' ...
                'name. Use ''default'' (or ''spm''/''rayint''), ''exact'' or ''sampling''.'], method);
            return
        end
    elseif hasRetrieval && any(sn.nodetype == NodeType.Source)
        if ~any(strcmpi(method,{'default','exact','rayint'}))
            reason = sprintf(['Method ''%s'' has no route on a delayed-hit retrieval cache: ' ...
                'solver_nc_retrieval_analyzer would answer it with the exact recurrences under ' ...
                'that name. Use ''default'', ''exact'' or ''rayint''.'], method);
            return
        end
    elseif any(strcmpi(method,{'rayint','spm'})) && sn.nclosedjobs > 0
        reason = sprintf(['Method ''%s'' names the SPM saddle point of an isolated cache; on an ' ...
            'integrated cache-queueing model the name is handed to the network solve, where ' ...
            'pfqn_nc has no arm for it. Use ''default'' (SPM cache, default network) or ''exact''.'], method);
        return
    end
    if strcmpi(method,'exact') && ~hasRetrieval
        if isstruct(sn.nodeparam{ci}) && isfield(sn.nodeparam{ci},'replacestrat')
            rs = sn.nodeparam{ci}.replacestrat;
            % CACHE_PROB_EREC is exact for the exchangeable (RR/FIFO) family
            % only; anything else has to take the approximate route.
            if ~any(rs == [ReplacementStrategy.RR, ReplacementStrategy.FIFO])
                reason = ['NC does not support exact solution of the specified cache ' ...
                    'replacement policy; use the default (approximate) method or SolverCTMC.'];
            end
        end
    end
    return
end
if any(strcmpi(method,{'rayint','spm'}))
    reason = sprintf(['Method %s names the SPM saddle point of a cache and, on a retrieval ' ...
        'model, the ray/WKB delayed-hit expansion; this model declares no Cache node.'], method);
    return
end

% -- single-station M/M/1/K with tail drop -------------------------------
% Answered exactly by the probability-based QSYS_MM1K_LOSS branch under
% 'default' and 'exact', and by the censored GE/GE/1/N block of the maximum
% entropy route under 'mem' (@SolverNC/runAnalyzer lets that name past the
% closed form). No other name has a route: the closed form reads none.
if sn_is_mm1k_loss(sn)
    if ~any(strcmpi(method,{'default','exact','mem'}))
        reason = sprintf(['Method ''%s'' has no route on a single-station M/M/1/K with tail ' ...
            'drop, which is answered by the closed form of qsys_mm1k_loss under ''default'' and ' ...
            '''exact'' and by the censored GE/GE/1/N block under ''mem''.'], method);
    end
    return
end

% -- loss networks and finite capacity regions ---------------------------
[isLossn, hasLossnShape] = nc_is_lossn_model(sn);
if isLossn
    % SOLVER_NC_LOSSN_ANALYZER switches on its six method names and answers any other
    % name with the residue transform ('exact') under that name.
    if ~any(strcmpi(method,{'default','exact','ms','erlangfp','rec','mci'}))
        reason = sprintf(['Method ''%s'' has no route on a loss network, which ' ...
            'solver_nc_lossn_analyzer would answer with the residue transform under that name. ' ...
            'Use ''default'', ''exact'' (or ''ms''), ''rec'', ''erlangfp'' or ''mci''.'], method);
    end
    return
end
if hasLossnShape
    reason = ['SolverNC does not support finite capacity regions with WAITQ ' ...
        '(blocking) policy. Use DROP policy instead.'];
    return
end
if any(strcmpi(method,{'ms','erlangfp'}))
    reason = sprintf(['Method %s is admissible only on a loss network (open model, one DROP ' ...
        'region holding a single Delay).'], method);
    return
end
if strcmpi(method,'rec')
    reason = ['Method rec is the MDD-rec route, admissible on a stochastic Petri net or on ' ...
        'a loss network (open model, one DROP region holding a single Delay); this model is neither.'];
    return
end
if isfield(sn,'nregions') && ~isempty(sn.nregions) && sn.nregions > 0
    % NC does not enforce an aggregate region limit on queueing stations;
    % refuse rather than silently return the unconstrained answer.
    reason = ['This model uses a Finite Capacity Region (addRegion) on queueing stations, ' ...
        'which is not supported by SolverNC. Use SolverJMT, or setCapacity for a ' ...
        'single-station limit.'];
    return
end

% -- class- and joint-dependent rates beside a server count or a lattice --
% A model with sn.cdscaling or sn.jdscaling is diverted, on every route, to
% SOLVER_NC_CONV, Sauer's multichain convolution, which reads sn.nservers only
% to tell a delay from a queue and never reads sn.lldscaling: a c-server or
% load-dependent station beside a class-dependent one would be solved as a
% single fixed-rate server. Refused rather than answered wrongly, the same
% pattern as the residual FCR above.
if ~isempty(sn.cdscaling) || ~isempty(sn.jdscaling)
    if any(isfinite(sn.nservers) & sn.nservers > 1) || ~isempty(sn.lldscaling)
        reason = ['SolverNC solves class- or joint-dependent rates by the convolution of ' ...
            'solver_nc_conv, which reads no server count above one and no load-dependent rate ' ...
            'lattice; fold the multiserver capacity into the class-dependent handle, or use ' ...
            'SolverMVA or SolverCTMC.'];
        return
    end
end

% -- PANACEA's domain ----------------------------------------------------
% Normal usage is a property of the demands rather than of a declared
% construct, so it has no feature name; an open chain is refused earlier by the
% closed-population feature set of the load-dependent evaluators.
%
% BOTH TOKENS ARE GATED, because PFQN_NCLD evaluates 'pana' and 'panald'
% with the SAME PFQN_PANACEALD -- its case label is {'pana','panald'} --
% so on a model carrying a rate lattice the load-INDEPENDENT name reaches the
% load-dependent expansion and raises with it. Off that lattice 'pana' takes
% its own PFQN_NC arm, which warns and returns an empty constant rather than
% raising, so it is left alone there. Class- or joint-dependent scaling diverts
% the whole model to SOLVER_NC_CONV, which never reads the method at all.
if any(strcmpi(method,{'pana','panald'})) && ~sn_has_open_classes(sn)
    divertedToConv = ~isempty(sn.cdscaling) || ~isempty(sn.jdscaling);
    reachesLdKernel = strcmpi(method,'panald') || ~isempty(sn.lldscaling);
    if ~divertedToConv && reachesLdKernel && ~nc_is_normal_usage(sn)
        reason = ['The model is not in normal usage, so the ''panald'' asymptotic expansion ' ...
            'does not apply. Use ''exact'', ''clw'' or an approximate load-dependent method instead.'];
        if strcmpi(method,'pana')
            reason = ['Method ''pana'' reaches the load-dependent kernel on this model, ' ...
                'where pfqn_ncld evaluates it as ''panald''. ' reason];
        end
        return
    end
    % Off the lattice 'pana' is gated for the REPORT ONLY, like 'mmint2':
    % PFQN_PANACEA returns NaN outside normal usage and PFQN_NC answers with an
    % empty constant, so the run renders an empty table under that name. The
    % question is asked in the Seidmann form SOLVER_NC hands PFQN_NC.
    if forReport && ~divertedToConv && ~reachesLdKernel && ~nc_is_normal_usage(sn, 'seidmann')
        reason = ['The model is not in normal usage, so the ''pana'' asymptotic expansion ' ...
            'does not apply (pfqn_panacea returns NaN and the analyzer an empty table). Use ' ...
            '''exact'', ''clw'' or ''default'' instead.'];
        return
    end
end

% -- the single-queueing-station recursions -------------------------------
% Two families are stated for a model with a delay and ONE queueing station, and
% neither can say so with a feature name: it is a COUNT, and a feature set has
% no arithmetic. PFQN_NC states it for 'mmint2'/'gleint' in those words and
% PFQN_COMOMRM_LD raises 'The solver accepts at most a single queueing station.'
% 'comom' is the same recursion at fixed rates: PFQN_NC raises for it on a
% multiclass model with more than one queue, and on a single-class one the
% queue-length step augments the model with an auxiliary class and raises
% there, so the rule is the station count alone.
%
% The count is taken over the CLOSED chains only, and the rule is inactive
% without a closed population, because PFQN_NC answers an open network with the
% exact open formulas BEFORE its method switch -- the method name is never read there,
% so a purely open model with three queues runs these names correctly today and
% must go on doing so.
% 'mmint2' and 'gleint' are gated for the REPORT ONLY: PFQN_NC answers them with
% an empty constant and the caller renders a table of zeros, which is a pair the
% report must not offer and a run the reference nonetheless performs. See
% FORREPORT above.
if any(strcmpi(method,{'comom','comomld'})) || (forReport && any(strcmpi(method,{'mmint2','gleint'})))
    nq = nc_closed_queueing_stations(sn);
    if nq > 1
        if strcmpi(method,'comomld')
            reason = sprintf(['Method ''comomld'' is the load-dependent CoMoM recursion, and ' ...
                'pfqn_comomrm_ld accepts at most a single queueing station; this model has %d.'], nq);
        elseif strcmpi(method,'comom')
            reason = sprintf(['Method ''comom'' is the CoMoM recursion of pfqn_comomrm, which ' ...
                'accepts at most a single queueing station; this model has %d. Use ''default'' ' ...
                'or ''ca'' for an exact normalizing constant, or SolverJMT with method ''jmva.recal''.'], nq);
        else
            reason = sprintf(['The ''%s'' method requires a model with a delay and a single ' ...
                'queueing station; this model has %d.'], method, nq);
        end
        return
    end
end

% -- the closed forms that read fixed rates only ---------------------------
% 'divdiff' is stated for single-server load-independent queues. A c-server
% station reaches PFQN_NC as Seidmann's split, demand L/c plus a surrogate
% delay L(c-1)/c, and the delay is what the closed form has no integral for
% (Cor. 3.4), so PFQN_NC raises there. A server count has no feature name.
if strcmpi(method,'divdiff') && any(isfinite(sn.nservers) & sn.nservers > 1)
    reason = ['Method ''divdiff'' is the divided-difference closed form for single-server ' ...
        'queues; a multiserver station enters the normalizing constant as Seidmann''s surrogate ' ...
        'delay, which needs the integral form of Corollary 3.4. Use ''ca'', ''exact'' or ''default''.'];
    return
end
% The three fixed-rate names PFQN_NCLD also serves are reached by SOLVER_NCLD
% on its CLOSED branch only: an open chain beside a rate lattice sends the
% model to the mixed solver (pfqn_mvaldmx, reported as 'ncldmx'), which never
% reads the method, so the answer would come back under the caller's name. A
% conjunction (open AND load-dependent) is not a feature-set delta.
if any(strcmpi(method,{'divdiff','clw','pana'})) && sn_has_open_classes(sn) && ~isempty(sn.lldscaling)
    reason = sprintf(['Method ''%s'' has no route on an open or mixed load-dependent model: ' ...
        'solver_ncld sends an open chain to the mixed solver (pfqn_mvaldmx, reported as ' ...
        '''ncldmx''), which never reads the method name. Use ''default'' or ''exact''.'], method);
    return
end

% -- 'exact' outside its domain ------------------------------------------
if strcmpi(method,'exact')
    nservers = sn.nservers;
    hasMultiserver = any(nservers(isfinite(nservers)) > 1);
    if hasMultiserver && any(isinf(sn.njobs))
        reason = ['NC solver cannot provide exact solutions for open or mixed queueing ' ...
            'networks. Remove the ''exact'' option.'];
        return
    end
    if hasMultiserver && any(sn.nodetype == NodeType.Fork)
        % fjFixedPoint hands ncDispatch the MMT image, whose parallelism rides
        % on OPEN auxiliary classes: the inner model is mixed, and
        % SOLVER_NC_ANALYZER refuses 'exact' there in the words above.
        reason = ['NC solver cannot provide exact solutions for a fork-join model with a ' ...
            'multiserver station: the fork-join transformation carries the parallelism as open ' ...
            'auxiliary classes, so the inner model is mixed. Remove the ''exact'' option.'];
        return
    end
    hasScaling = ~isempty(sn.lldscaling) || ~isempty(sn.cdscaling) || ~isempty(sn.jdscaling);
    Nfin = sn.njobs(isfinite(sn.njobs));   % an open class has no lattice position
    if (hasScaling || hasMultiserver) && any(abs(Nfin - floor(Nfin)) > GlobalConstants.FineTol)
        % The load-dependent analyzer interpolates a fractional population
        % between the two integer neighbours, which is an approximation, so it
        % refuses the exactness the caller asked for by name.
        reason = ['NC load-dependent solver cannot provide exact solutions for fractional ' ...
            'populations.'];
        return
    end
end
end

function nq = nc_closed_queueing_stations(sn)
% NQ = NC_CLOSED_QUEUEING_STATIONS(SN)
% How many queueing (non-infinite-server) stations carry demand from a CLOSED
% chain, which is the row count L reaches PFQN_NC and PFQN_COMOMRM_LD with once
% the delay rows have been folded into Z and the zero-demand rows dropped.
% Zero when the model has no closed population at all.
nq = 0;
[Lchain,~,~,~,Nchain] = sn_get_demands_chain(sn);
closed = isfinite(Nchain(:)') & Nchain(:)' > 0;
if ~any(closed)
    return
end
for ist = 1:sn.nstations
    if isinf(sn.nservers(ist))
        continue
    end
    if any(abs(Lchain(ist,closed)) > GlobalConstants.FineTol)
        nq = nq + 1;
    end
end
end

function tf = nc_is_noreentrant_cache(sn)
% TF = NC_IS_NOREENTRANT_CACHE(SN)
% The Source-Cache-Sink model SOLVER_NC_CACHE_ANALYZER serves, as
% @SolverNC/runAnalyzer identifies it.
tf = sn.nclosedjobs == 0 && numel(sn.nodetype) == 3 ...
    && isequal(sort(sn.nodetype(:))', sort([NodeType.Source, NodeType.Cache, NodeType.Sink]));
end
