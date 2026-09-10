classdef SolverMVA < NetworkSolver
    % Mean Value Analysis solver for queueing networks
    %
    % Implements MVA algorithms for analyzing closed and open queueing networks.
    % Method 'amva.mapqn' is the horizontal-cut mean value analysis of a closed
    % model with one exponential delay and one FCFS queue whose service is a
    % MAP per class (mapqn_amva), the mean-value counterpart of the MAP-AMVA
    % bounds of SolverBA.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = SolverMVA(model,varargin)
            % SOLVERMVA Create an MVA solver instance
            %
            % @brief Creates a Mean Value Analysis solver for the given model
            % @param model Network model to be analyzed
            % @param varargin Optional solver options (method, tolerance, etc.)
            % @return self SolverMVA instance configured with specified options

            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, SolverMVA.defaultOptions));
            self.setLang();
        end

        function sn = getStruct(self)
            % GETSTRUCT Get model data structure for analysis
            %
            % @brief Returns the internal data structure representing the model
            % @return sn Structured data representing the queueing network
            sn = self.model.getStruct(false);
        end

        function tf = supportsExactSensitivity(self) %#ok<MANU>
            % TF = SUPPORTSEXACTSENSITIVITY()
            % MVA differentiates its own recursion: getSensitivityTable uses
            % the analytic branch (pfqn_sens) wherever the model is in scope.
            tf = true;
        end

        [runtime, analyzer] = runAnalyzer(self, options);
        [lNormConst] = getProbNormConstAggr(self);
        [Pnir,logPnir] = getProbAggr(self, ist);
        [Pnir,logPn] = getProbSysAggr(self);

        function [allMethods] = listValidMethods(self)
            % LISTVALIDMETHODS Get all valid MVA solution methods
            %
            % @brief Returns cell array of valid MVA methods for the current model
            % @return allMethods Cell array of method names available for this model

            sn = self.model.getStruct;
            % base set of methods
            allMethods = {'default',...
                'mva','exact','amva', ...
                'sum','esum', ...
                'qdlin','amva.qdlin', ...
                'bs','amva.bs', ...
                'aql','amva.aql', ...
                'qsa','amva.qsa', ...
                'sqni', ...
                'qd','amva.qd', ...
                'qli','amva.qli', ...
                'fli','amva.fli', ...
                'ab','amva.ab', ...
                'schmidt','amva.schmidt', ...
                'schmidt-ext','amva.schmidt-ext', ...
                'tay','amva.tay', ...
                'scat','amva.scat', ...
                'lcp','amva.lcp', ...
                'chow','amva.chow', ...
                'pamb','amva.pamb', ...
                'pami','amva.pami', ...
                'pamt','amva.pamt', ...
                'clust','amva.clust', ...
                'dmlin','amva.dmlin'};

            % SQNI (pfqn_sqni) is a closed form for one queueing station with a
            % delay; listing it elsewhere named a method that cannot run.
            if ~(sn.nstations==2 && sum(sn.sched==SchedStrategy.INF)==1)
                allMethods(strcmp(allMethods,'sqni')) = [];
            end

            % AQL (pfqn_aql), QSA (pfqn_qsa) and Tay (pfqn_tay) reject
            % multiserver stations at solve time (solver_amva.m), so they are
            % only advertised for single-server models.
            if sn_has_multi_server(sn)
                allMethods(strcmp(allMethods,'aql')) = [];
                allMethods(strcmp(allMethods,'amva.aql')) = [];
                allMethods(strcmp(allMethods,'qsa')) = [];
                allMethods(strcmp(allMethods,'amva.qsa')) = [];
                allMethods(strcmp(allMethods,'tay')) = [];
                allMethods(strcmp(allMethods,'amva.tay')) = [];
            end

            % QNA/RQNA/RQT advertised for fully open models only; see
            % _kb/06-solver-catalog.md (MVA section) for the default-dispatch rules
            if sn_is_open_model(sn)
                allMethods = [allMethods(1:4), {'qna','rqna','rqt'}, allMethods(5:end)];
            end

            % SQD (Smith Queue Decomposition) is only valid for closed
            % single-chain Blocking-After-Service networks; kept at this
            % position to preserve BAS regression-baseline ordering.
            if sn_is_bas_model(sn)
                allMethods{end+1} = 'sqd'; %#ok<AGROW>
            end

            % MVAC (exact mean value analysis by chain, pfqn_mvac): closed
            % single-server product-form networks only; rejects open/mixed and
            % multiserver at solve time.
            if ~any(isinf(sn.njobs))
                allMethods{end+1} = 'mvac'; %#ok<AGROW>
            end

            % SJN (shortest-job-next, pfqn_mvasjn / pfqn_amvasjn): the
            % conditional waiting time equation is a population recursion, so
            % the family runs on a CLOSED model with an SJF station and nowhere
            % else -- mvaDispatch rejects an open one by name.
            % Advertised only there, for the reason 'sqni' is gated above.
            if ~any(isinf(sn.njobs)) && any(sn.sched == SchedStrategy.SJF)
                allMethods = [allMethods, {'sjn.mva','sjn.amva'}]; %#ok<AGROW>
            end

            % Marie withheld for open models and for class-dependent routing;
            % see _kb/06-solver-catalog.md (MVA section) for the rationale
            if ~sn_is_open_model(sn) && ~sn_has_classdep_routing(sn)
                allMethods = [allMethods, {'marie','amva.marie'}]; %#ok<AGROW>
            end

            % amva.mapqn: the horizontal-cut MVA for one exponential delay and one
            % FCFS MAP queue (mapqn_amva); offered only on that shape, which
            % mva_mapqn_reason judges for the list, the report and the run
            if isempty(mva_mapqn_reason(sn))
                allMethods = [allMethods, {'amva.mapqn'}]; %#ok<AGROW>
            end

            % priomva: preemptive-resume priority arm (Chandy-Lakshmi [ChaL83]),
            % offered only when a station actually uses FCFSPRPRIO.
            if any(sn.sched == SchedStrategy.FCFSPRPRIO)
                allMethods = [allMethods, {'priomva','amva.priomva'}]; %#ok<AGROW>
            end

            allMethods = {allMethods{:}, ...
                'lin','egflin','gflin','amva.lin'}; %#ok<CCAT>

            % Bound methods (aba/bjb/gb/pb/sb/mwba) are NOT listed here: they
            % moved to SolverBA, and runAnalyzer raises line_error for the whole
            % family. Listing them made this a false claim, since every listed
            % name errored when run. See SolverBA.listValidMethods.

            if sn_is_open_model(sn) && sn.nstations == 2 && sn.nclasses == 1
                % methods to add for queueing systems
                qsys = {'mm1','mmk','mg1','mgi1','gm1','gig1','gim1','gig1.kingman', ...
                    'gigk','gigk.kingman_approx', ...
                    'gig1.gelenbe','gig1.heyman','gig1.kimura','gig1.allen', ...
                    'gig1.kobayashi','gig1.klb','gig1.marchal', ...
                    'gigk.whitt','qed','gig1.extremal','gigk.diffusion'};
                % The two abandonment methods are listed only when the station
                % actually reneges: they have nothing to say about a queue
                % nobody walks away from, and listing them there would name a
                % method that cannot run.
                queue_ist = sn.nodeToStation(sn.nodetype == NodeType.Queue);
                if ~isempty(queue_ist) && ~isempty(sn_patience_handles(sn, queue_ist, 1))
                    qsys = {qsys{:}, 'erlanga', 'mgisrgi'}; %#ok<CCAT>
                end
                % append, keeping original order and avoiding duplicates
                allMethods = {allMethods{:}, qsys{:}}; %#ok<CCAT>
            end
        end

        function method = resolveMethod(self, options)
            % Feature-driven resolution of options.method='default'. A bursty
            % single-class open network has a non-renewal (MAP/MMPP) arrival
            % process whose autocorrelation a two-moment method cannot capture,
            % so the default dispatch selects RQNA (robust queueing network
            % analyzer, indices of dispersion). All other cases keep 'default',
            % which the analyzer expands with its own heuristics. Expressed as
            % RQNA's MAP-family coverage rather than a bespoke gate.
            method = options.method;
            if strcmp(options.method, 'default')
                sn = self.model.getStruct();
                if (sn.nclasses == 1) && all(isinf(sn.njobs)) && sn_has_bursty_arrival(sn)
                    method = 'rqna';
                elseif (sn.nclasses == 1) && all(isinf(sn.njobs)) && sn.nstations == 2
                    % A single-class open station customers ABANDON is a
                    % different model, not a correction to a G/G/k one: the
                    % resolution has to happen here as well as in the analyzer,
                    % because the feature gate runs on the resolved name and
                    % Reneging is admitted for these two methods only.
                    queue_ist = sn.nodeToStation(sn.nodetype == NodeType.Queue);
                    if ~isempty(queue_ist)
                        hpat = sn_patience_handles(sn, queue_ist, 1);
                        if ~isempty(hpat)
                            if hpat.isExponential
                                method = 'erlanga';
                            else
                                method = 'mgisrgi';
                            end
                        end
                    end
                end
            end
        end

        function featSupported = getMethodFeatureSet(self, method)
            % Per-method feature deltas applied to the base MVA envelope.
            % RQNA natively consumes non-renewal MAP/MMPP/RAP arrival and
            % service processes (open only); QNA is a two-moment open-network
            % method. The queueing-system and bounds methods are already
            % structurally restricted by listValidMethods, so they inherit the
            % base envelope unchanged.
            featSupported = SolverMVA.getFeatureSet();
            method = regexprep(method, '^amva\.', '');
            if SolverMVA.isClosedPopulationMethod(method)
                % The closed-population AMVA family estimates the
                % arrival-instant queue length as a function of the population
                % vector N and is handed (L,N,Z) alone, so an open chain gives
                % it nothing to recur on: solver_amva has no arm for any of
                % these outside its closed product-form branch, and falling
                % through returned the qd-family answer, or a table of zeros,
                % under their name. The remaining preconditions of that branch
                % (load dependence, product form) have no registry name and are
                % applied by SUPPORTSCLOSEDPOPULATION instead.
                featSupported.setFalse({'OpenClass'});
            end
            % solver_mvald_analyzer serves a load-, class- or joint-dependent
            % model through 'exact'/'mva' (load dependence only, it has no
            % class- or joint-dependent recursion) and through the
            % default/amva/qd/lin/qdlin arms, and refuses every other name by
            % name. The queueing-system closed forms are intercepted by
            % mvaDispatch upstream of that analyzer and keep the base envelope.
            if any(strcmp(method, [SolverMVA.closedPopulationMethods(), ...
                    {'sum','esum','mvac','qli','fli','gflin','egflin','qna','rqna','rqt', ...
                    'priomva','sqd','marie','sjn.mva','sjn.amva'}]))
                featSupported.setFalse({'LoadDependence','ClassDependence','JointDependence'});
            elseif any(strcmp(method, {'mva','exact'}))
                featSupported.setFalse({'ClassDependence','JointDependence'});
            end
            if ~any(strcmp(method, {'default','exact'}))
                % An order-independent or pass-and-swap station is served by
                % solver_mva_oi_analyzer alone, which mvaDispatch reaches only
                % under 'default' or 'exact'; every other name is refused there
                % by name, so it must not be advertised for such a model.
                featSupported.setFalse({'SchedStrategy_OI','SchedStrategy_PAS'});
            end
            if strcmp(method, 'mva')
                % solver_mva has arms for INF, PS, FCFS, SIRO and LCFS-PR (plus
                % the paired LCFS recursion) and refuses every other discipline
                % by name; the single-station shapes mvaDispatch claims ahead of
                % it (Cobham, size-based M/G/1, M/M/1-DPS) are not the MVA
                % recursion and are reached through 'default' or 'exact'.
                % LCFS and SJF stay: SUPPORTSLCFS and SUPPORTSSJN gate them.
                featSupported.setFalse(setdiff(SolverMVA.nonBcmpSchedFeatures(), ...
                    {'SchedStrategy_LCFS','SchedStrategy_SJF'}));
            end
            switch method
                case {'sum','esum'}
                    % solver_mva_sum passes each station to sum_closed /
                    % sum_closing as an INF, PS, LCFS-PR, FCFS or SIRO centre
                    % and refuses every other discipline by name.
                    featSupported.setFalse(SolverMVA.nonBcmpSchedFeatures());
                case 'mvac'
                    % pfqn_mvac recurs on the closed chains over single-server
                    % fixed-rate (SSFR) and infinite-server centres; solver_mvac
                    % refuses every other discipline by name.
                    featSupported.setFalse([{'OpenClass'}, SolverMVA.nonBcmpSchedFeatures()]);
                case 'exact'
                    % the exact cache recursion (cache_mva) is the RR/FIFO one:
                    % solver_mva_cache_analyzer refuses LRU and h-LRU by name
                    % and the cache-in-network miss callback runs it whatever
                    % the policy, so neither may be admitted under 'exact'
                    featSupported.setFalse({'ReplacementStrategy_LRU','ReplacementStrategy_HLRU'});
                case {'schmidt','schmidt-ext'}
                    % pfqn_schmidt has INF, PS and FCFS arms only: an LCFS-PR
                    % station fell through its switch and kept a zero wait
                    featSupported.setFalse({'SchedStrategy_LCFSPR'});
                case 'rqna'
                    featSupported.setTrue({'MAP','MMPP2','MMAP','RAP'});
                    featSupported.setFalse({'ClosedClass','SelfLoopingClass'});
                case 'marie'
                    % solver_mva_marie_analyzer folds INF stations into the
                    % think time and isolates FCFS, PS and LCFS-PR queues from
                    % chain demands and SCVs alone: closed networks, no cache,
                    % no fork-join, no scaling. The rest of its premise (class
                    % routing, multiserver with several chains) is structural,
                    % see mva_marie_reason.
                    featSupported.setFalse([{'OpenClass','Source','Sink', ...
                        'Fork','Forker','Join','Joiner','JoinPartial', ...
                        'Cache','CacheClassSwitcher','CacheRetrieval', ...
                        'SchedStrategy_SIRO'}, SolverMVA.nonBcmpSchedFeatures()]);
                case {'sjn.mva','sjn.amva'}
                    % solver_mva_sjn_analyzer: a closed model whose SJF station
                    % and other single-server BCMP stations enter Kant's
                    % recursion through chain demands and SCVs; see SUPPORTSSJN
                    % for the server-count half of the rule.
                    featSupported.setFalse([{'OpenClass','Source','Sink', ...
                        'Fork','Forker','Join','Joiner','JoinPartial', ...
                        'Cache','CacheClassSwitcher','CacheRetrieval'}, ...
                        setdiff(SolverMVA.nonBcmpSchedFeatures(), {'SchedStrategy_SJF'})]);
                case {'mm1','mmk','qed'}
                    % M/M/1, M/M/k and the Halfin-Whitt limit read the two
                    % means alone; under any other law they answered the
                    % exponential system under the caller's name. BMAP is
                    % dropped too: the MX/M/1 branch would preempt them.
                    featSupported.setFalse(SolverMVA.nonExponentialFeatures());
                case 'mapqn'
                    % the horizontal-cut MVA consumes a MAP service natively (a
                    % closed delay + FCFS queue model, see mva_mapqn_reason);
                    % declaring MAP here is what keeps NetworkSolver.needsMapEnv
                    % from routing the model through its random-environment image
                    featSupported.setTrue({'MAP','MMPP2'});
                    featSupported.setFalse([{'OpenClass','Source','Sink','Fork','Forker','Join','Joiner','JoinPartial', ...
                        'ClassSwitch','StatelessClassSwitcher','Cache','CacheClassSwitcher','CacheRetrieval', ...
                        'LoadDependence','ClassDependence','JointDependence', ...
                        'SchedStrategy_PS','SchedStrategy_SIRO','SchedStrategy_LCFSPR', ...
                        'SchedStrategy_SRPT','SchedStrategy_PSJF','SchedStrategy_FB','SchedStrategy_LRPT', ...
                        'SchedStrategy_SETF','SchedStrategy_OI','SchedStrategy_PAS'}, SolverMVA.nonBcmpSchedFeatures()]);
                case 'rqt'
                    % robust queueing theory: single-class open networks, the
                    % primitives entering the uncertainty sets are two moments
                    featSupported.setFalse({'ClosedClass','SelfLoopingClass'});
                case {'erlanga','mgisrgi'}
                    % The only analytical methods in LINE that carry an
                    % abandonment rate. Reneging stays OUT of the base MVA
                    % envelope: every other method here would silently ignore
                    % the patience law and report the no-abandonment answer.
                    featSupported.setTrue({'Reneging'});
                    featSupported.setFalse({'ClosedClass','SelfLoopingClass'});
                    % The MX/M/1 branch of the analyzer preempts every name on
                    % a BMAP source and carries no abandonment. Erlang A is the
                    % M/M/s+M chain and reads the two means alone; M/GI/s/r+GI
                    % keeps the service law (the Poisson source is structural,
                    % see SUPPORTSQUEUEINGSYSTEM).
                    if strcmp(method, 'erlanga')
                        featSupported.setFalse(SolverMVA.nonExponentialFeatures());
                    else
                        featSupported.setFalse({'BMAP'});
                    end
                case 'qna'
                    % round-robin dispatching enters as a deterministic
                    % traffic split (npfqn_traffic_split_rr), which the
                    % exact-MVA paths have no counterpart for
                    featSupported.setTrue({'RoutingStrategy_RROBIN'});
                    featSupported.setFalse({'ClosedClass','SelfLoopingClass'});
                    % solver_qna's station loop has an arm for INF, PS and FCFS
                    % and none for anything else, so on a SIRO, LCFS-PR, HOL or
                    % priority station it left that row of Q, U, R and T at
                    % zero and reported the table as a solution.
                    featSupported.setFalse([{'SchedStrategy_SIRO','SchedStrategy_LCFSPR'}, ...
                        SolverMVA.nonBcmpSchedFeatures()]);
            end
            if any(strcmp(method, {'rqna','rqt'}))
                % A Join is a synchronisation node, not a queue: it carries no
                % service process, so the index-of-dispersion curve these two
                % read off every station does not exist for it, and neither
                % analyzer has a synchronisation term to put in its place. QNA
                % keeps Fork/Join -- its station loop has an explicit Join arm.
                featSupported.setFalse({'Fork','Forker','Join','Joiner','JoinPartial'});
                % Both decompose the network into FCFS queues (GI/G/1 workload
                % for RQNA, the worst-case G/G/k system time for RQT) and read
                % no discipline: a PS, SIRO, LCFS-PR or priority station was
                % answered as FCFS. Neither has a cache arm either; the
                % cache-in-network branch of mvaDispatch would hand them the
                % relabelled network and they would refuse it by class count.
                featSupported.setFalse([{'SchedStrategy_PS','SchedStrategy_SIRO','SchedStrategy_LCFSPR', ...
                    'Cache','CacheClassSwitcher','CacheRetrieval'}, SolverMVA.nonBcmpSchedFeatures()]);
            end
            % MULTISERVER (registry name since 2026-09-05). The single-server
            % recursions: AQL, QSA and Tay (SUPPORTSCLOSEDPOPULATION), MVAC's
            % SSFR chain recursion (SUPPORTSMVAC), RQNA's GI/G/1 workload
            % (SUPPORTSSINGLECLASSOPEN), Kant's SJN recursion (SUPPORTSSJN) and
            % the single-server closed forms of the queueing-system analyzer,
            % every M/G/1, G/M/1 and G/G/1 name (SUPPORTSQUEUEINGSYSTEM). Each
            % predicate stays, wording the refusal for the run; the delta is
            % what makes it nameable. RQT, QNA, M/M/k, G/G/k and the rest of
            % the envelope carry a server count.
            if any(strcmp(method, {'aql','qsa','tay','mvac','rqna','sjn.mva','sjn.amva', ...
                    'mm1','mg1','mgi1','gm1','gim1'})) || strncmp(method, 'gig1', 4)
                featSupported.setFalse({'MultiServer'});
            end
            % FINITECAPACITY (registry name since 2026-09-05) is NOT in the
            % base envelope: the product-form recursions solve a buffer away,
            % which is what SUPPORTSFINITECAPACITY refuses. The names that
            % honour one are granted it here, and the structural predicate
            % keeps the shape half of each rule. 'default' and 'sqd' reach
            % solver_sqd, the one Blocking-After-Service arm. The single-station
            % M/M/1/K with tail drop is served by the qsys_mg1k_loss_mgs branch
            % of solver_mva_qsys_analyzer, which mvaDispatch reaches under EVERY
            % name on that shape (it ORs sn_is_mm1k_loss into the single-station
            % test), 'exact' excepted since the branch is exact at scv=1 only;
            % that grant is judged on the model because no name can carry it.
            if any(strcmp(method, {'default','sqd'}))
                featSupported.setTrue({'FiniteCapacity'});
            elseif ~strcmp(method, 'exact') && isa(self.model, 'Network') ...
                    && sn_is_mm1k_loss(self.model.getStruct())
                featSupported.setTrue({'FiniteCapacity'});
            end
        end

        function [bool, reason] = supportsModelMethod(self, method)
            % Finite station/class capacity has no registry feature name, so
            % the coarse per-method feature gate cannot see it. Apply the
            % structural capacity check on top of it, otherwise MVA silently
            % returns the unconstrained product-form answer for models built
            % with setCapacity / a finite classCap (BUG-39).
            [bool, reason] = supportsModelMethod@NetworkSolver(self, method);
            if bool && isa(self.model, 'Network')
                [bool, reason] = SolverMVA.supportsFiniteCapacity(self.model, method);
            end
            if bool && isa(self.model, 'Network')
                [bool, reason] = SolverMVA.supportsExactness(self.model, method);
            end
            if bool && isa(self.model, 'Network')
                sn = self.model.getStruct();
                [bool, reason] = SolverMVA.supportsClosedPopulation(sn, method);
                if bool
                    [bool, reason] = SolverMVA.supportsSingleClassOpen(sn, method);
                end
                if bool
                    [bool, reason] = SolverMVA.supportsMvac(sn, method);
                end
                if bool
                    [bool, reason] = SolverMVA.supportsMapqn(sn, method);
                end
                if bool
                    [N, isFcfs] = SolverMVA.schmidtArmInputs(sn, method);
                    [bool, reason] = SolverMVA.supportsSchmidtExt(N, isFcfs, method);
                end
                % The single-analyzer disciplines and shapes: each predicate is
                % the one its analyzer (or mvaDispatch) calls at run time.
                if bool
                    [bool, reason] = SolverMVA.supportsLcfs(sn, method);
                end
                if bool
                    [bool, reason] = SolverMVA.supportsSjn(sn, method);
                end
                if bool
                    [bool, reason] = SolverMVA.supportsOrderIndependent(sn, method);
                end
                if bool
                    [bool, reason] = SolverMVA.supportsPolling(sn, method);
                end
                if bool
                    [bool, reason] = SolverMVA.supportsSizeBased(sn, method);
                end
                if bool
                    [bool, reason] = SolverMVA.supportsQueueingSystem(sn, method);
                end
                if bool
                    [bool, reason] = SolverMVA.supportsPreemptivePriority(sn, method);
                end
                if bool
                    [bool, reason] = SolverMVA.supportsMarie(sn, method);
                end
            end
        end

    end

    methods(Static)
        function names = closedPopulationMethods()
            % NAMES = CLOSEDPOPULATIONMETHODS()
            % The AMVA algorithms whose recursion is over a CLOSED population
            % vector: each approximates the arrival-instant queue length
            % E[Q(N-1_r)] from E[Q(N)] and is handed (L,N,Z) alone, with no
            % arrival rate and no load-dependent rate function. Canonical
            % spellings only; the 'amva.' aliases are stripped before this list
            % is consulted, exactly as solver_amva strips them before it
            % dispatches. ONE list, read by getMethodFeatureSet (which drops
            % OpenClass for them) and by solver_amva (which refuses them by
            % name), so the gate and the run cannot drift apart.
            names = {'bs','aql','qsa','sqni','tay','scat','lcp','chow', ...
                'pamb','pami','pamt','clust','dmlin','ab','schmidt','schmidt-ext'};
        end

        function tf = isClosedPopulationMethod(method)
            % TF = ISCLOSEDPOPULATIONMETHOD(METHOD)
            % True when METHOD, with any leading 'amva.' removed, names one of
            % the closed-population AMVA algorithms.
            tf = any(strcmpi(regexprep(method, '^amva\.', ''), ...
                SolverMVA.closedPopulationMethods()));
        end

        function names = nonExponentialFeatures()
            % NAMES = NONEXPONENTIALFEATURES()
            % The distribution names of the base MVA envelope other than Exp.
            % The M/M closed forms (mm1, mmk, qed, erlanga) drop them all: the
            % registry cannot tell an arrival law from a service law, and an
            % M/M formula reads only the means of either.
            names = {'APH','Coxian','Erlang','HyperExp','BMAP','MMAP', ...
                'Pareto','Weibull','Lognormal','Uniform','Det','Replayer'};
        end

        function names = nonBcmpSchedFeatures()
            % NAMES = NONBCMPSCHEDFEATURES()
            % The scheduling feature names OUTSIDE the BCMP set {INF, PS, FCFS,
            % SIRO, LCFS-PR} that the base MVA envelope declares. The chain
            % algorithms that walk the stations one by one (solver_mva_sum,
            % solver_mvac, solver_qna) accept the BCMP set and refuse the rest,
            % so each drops these from its own envelope.
            names = {'SchedStrategy_HOL','SchedStrategy_DPS','SchedStrategy_FCFSPRPRIO', ...
                'SchedStrategy_LCFS','SchedStrategy_POLLING','SchedStrategy_SJF', ...
                'SchedStrategy_SRPT','SchedStrategy_PSJF','SchedStrategy_FB', ...
                'SchedStrategy_LRPT','SchedStrategy_SETF', ...
                'SchedStrategy_OI','SchedStrategy_PAS'};
        end

        function [bool, reason] = supportsClosedPopulation(sn, method)
            % [BOOL, REASON] = SUPPORTSCLOSEDPOPULATION(SN, METHOD)
            % Open chains and strict product form are the preconditions of
            % the closed-population AMVA family that the analyzer must refuse by
            % name: solver_amva implements every one
            % of these algorithms in its product-form branch alone, and outside
            % it there is no arm to fall through to. Product form has no
            % registry feature name, so the check cannot live in
            % getMethodFeatureSet and belongs here. SOLVER_AMVA calls this same
            % predicate, so a model the report offers is one the analyzer runs.
            bool = true;
            reason = '';
            if ~SolverMVA.isClosedPopulationMethod(method)
                return
            end
            base = regexprep(method, '^amva\.', '');
            if sn_has_open_classes(sn)
                bool = false;
                reason = sprintf(['the ''%s'' method approximates the arrival-instant queue ' ...
                    'length as a function of the closed population vector N, so it is defined ' ...
                    'for closed models only; use ''default'', ''lin'', ''qd'' or ''qna'' for a ' ...
                    'model with open classes'], base);
                return
            end
            % Server counts and station counts have no registry name. AQL, QSA
            % and Tay are single-server recursions (solver_amva used to refuse
            % them inline, after the report had already offered them) and SQNI
            % is a closed form for one queueing station with a delay.
            if any(strcmpi(base, {'aql','qsa','tay'})) && sn_has_multi_server(sn)
                bool = false;
                reason = sprintf(['the ''%s'' method is defined for single-server stations; ' ...
                    'use ''default'' or ''lin'' for a model with a multiserver station'], base);
                return
            end
            if strcmpi(base, 'sqni') && ~(sn.nstations == 2 && sum(sn.sched == SchedStrategy.INF) == 1)
                bool = false;
                reason = ['SQNI is defined for a single queueing station with a delay. ' ...
                    'Try with the ''default'' or ''lin'' methods.'];
                return
            end
            % ab, schmidt and schmidt-ext ARE the class-dependent FCFS
            % algorithms, so heterogeneous FCFS service means are their subject
            % matter rather than a disqualification.
            checkMeans = ~any(strcmpi(base, {'ab','schmidt','schmidt-ext'}));
            if sn_has_product_form_not_het_fcfs(sn, checkMeans)
                return
            end
            bool = false;
            reason = sprintf(['the ''%s'' method is defined for strict product-form, ' ...
                'load-independent models; use ''default'', ''lin'' or ''qd'' for this model'], base);
        end

        function [bool, reason] = supportsSingleClassOpen(sn, method)
            % [BOOL, REASON] = SUPPORTSSINGLECLASSOPEN(SN, METHOD)
            % RQNA and RQT decompose an open network into GI/G/1 queues and
            % build one uncertainty set per flow from the first two moments of a
            % SINGLE stream, so a multiclass model has no counterpart in their
            % equations. There is no registry feature name for a class count, so
            % that half of the rule is structural. A FORK-JOIN model is refused
            % too, and that half IS nameable, so getMethodFeatureSet drops
            % Fork/Join for these two as well and this is the analyzer's half of
            % it. SOLVER_RQNA and SOLVER_RQT call this same predicate, so the
            % message the report gives is the message the run gives.
            bool = true;
            reason = '';
            if ~any(strcmpi(method, {'rqna','rqt'}))
                return
            end
            if any(sn.nodetype == NodeType.Fork | sn.nodetype == NodeType.Join)
                bool = false;
                reason = sprintf(['%s decomposes an open network into GI/G/1 queues and has ' ...
                    'no synchronisation term; a Join carries no service process for its index ' ...
                    'of dispersion to be read from. Use SolverMVA''s ''default'' method for a ' ...
                    'fork-join model.'], upper(method));
                return
            end
            if sn.nclasses ~= 1
                bool = false;
                reason = sprintf(['%s supports single-class open networks only. ' ...
                    'Use the ''qna'' method for multiclass models.'], upper(method));
                return
            end
            % RQNA's robust workload is the GI/G/1 one (qsys_gig1_rq, in the
            % network analyzer and in the single-station arm alike) and reads
            % no server count: a multiserver station was answered as a single
            % server of the same rate. RQT carries the count (qsys_gigk_rqt).
            if strcmpi(method, 'rqna')
                qset = find(sn.sched ~= SchedStrategy.INF & sn.sched ~= SchedStrategy.EXT);
                if any(isfinite(sn.nservers(qset)) & sn.nservers(qset) > 1)
                    bool = false;
                    reason = ['RQNA decomposes the network into GI/G/1 queues and has no ' ...
                        'multiserver workload; use ''rqt'', ''qna'' or ''default'' for a ' ...
                        'model with a multiserver station.'];
                end
            end
        end

        function [N, isFcfs] = schmidtArmInputs(sn, method)
            % [N, ISFCFS] = SCHMIDTARMINPUTS(SN, METHOD)
            % The population vector and per-row discipline that SOLVER_AMVA's
            % schmidt-ext arm hands the kernel: CHAIN populations, and which of
            % the queueing-station rows is served FCFS. Empty when METHOD is not
            % 'schmidt-ext', so the product-form parameters are built only for
            % the one method that reads them. The delay rows the arm stacks on
            % top are left out: an INF row is never an FCFS station.
            N = []; isFcfs = [];
            if ~strcmpi(regexprep(method, '^amva\.', ''), 'schmidt-ext')
                return
            end
            [~,~,N] = sn_get_product_form_chain_params(sn);
            queueIdx = sn.nodetype == NodeType.Queue;
            isFcfs = sn.sched(sn.nodeToStation(queueIdx)) == SchedStrategy.FCFS;
        end

        function [bool, reason] = supportsSchmidtExt(N, isFcfs, method)
            % [BOOL, REASON] = SUPPORTSSCHMIDTEXT(N, ISFCFS, METHOD)
            % Schmidt's EXTENSION over plain Schmidt is an alpha correction
            % applied at an FCFS station, and the correction is computed from
            % the network with ONE class-r customer TAGGED, that is at
            % population N - 1_r. A class holding no customer has none to tag:
            % the sub-problem is formed at a negative population, whose state
            % lattice prod(N+1) collapses to zero and the recursion indexes an
            % empty array. Plain 'schmidt' forms no such sub-problem, which is
            % why the requirement is the -ext arm's alone.
            % THE TEST IS STATED AT THE FCFS STATION AND NOT AT A CLASS-DEPENDENT
            % ONE, because the four kernels differ on when they form the
            % correction: this one, the C++ port and native python form it only
            % where the demands differ by class, the JAR forms it at every FCFS
            % station. Stating the union is what keeps one rule safe for all
            % four; the case it costs -- an FCFS station whose demands are
            % identical across classes, one of them empty -- is one where the
            % extension reduces to plain 'schmidt', which stays offered.
            % N and ISFCFS are the numbers the CALLER'S OWN arm passes, so
            % SOLVER_AMVA and the report ask the same question. They are
            % CHAIN-indexed here and in the C++ port, CLASS-indexed in the JAR
            % and in native python, whose arms aggregate no chains; that
            % difference belongs to those arms, not to this rule.
            bool = true;
            reason = '';
            if ~strcmpi(regexprep(method, '^amva\.', ''), 'schmidt-ext')
                return
            end
            if isempty(N) || isempty(isFcfs) || ~any(isFcfs)
                return
            end
            r = find(isfinite(N) & N < 1, 1);
            if isempty(r)
                return
            end
            bool = false;
            reason = sprintf(['the ''schmidt-ext'' method corrects an FCFS station from the ' ...
                'network with one customer of that class tagged, so it needs every class to ' ...
                'hold at least one customer; class %d holds none. Use ''schmidt'' for the ' ...
                'uncorrected recursion.'], r);
        end

        function [bool, reason] = supportsMvac(sn, method)
            % [BOOL, REASON] = SUPPORTSMVAC(SN, METHOD)
            % MVAC (Conway-de Souza e Silva-Lavenberg) is the exact chain
            % recursion over single-server fixed-rate (SSFR) queues and
            % infinite-server centres of a product-form network. Neither the
            % server count nor product form has a registry feature name, so both
            % are structural; the scheduling restriction IS nameable and lives
            % in getMethodFeatureSet. SOLVER_MVAC calls this same predicate.
            bool = true;
            reason = '';
            if ~strcmpi(method, 'mvac')
                return
            end
            if ~sn_has_product_form(sn)
                bool = false;
                reason = 'MVAC requires a product-form model.';
                return
            end
            qset = find(sn.sched ~= SchedStrategy.INF & sn.sched ~= SchedStrategy.EXT);
            if any(sn.nservers(qset) ~= 1)
                bool = false;
                reason = ['MVAC supports single-server (SSFR) queues only; ' ...
                    'use method ''exact'' for multiserver stations.'];
            end
        end

        function [bool, reason] = supportsMapqn(sn, method)
            % [BOOL, REASON] = SUPPORTSMAPQN(SN, METHOD)
            % The structural premise of 'amva.mapqn' (one exponential delay, one
            % FCFS single-server queue with Markovian service, serial
            % routing), judged by mva_mapqn_reason, the predicate the list
            % and the analyzer read as well.
            bool = true;
            reason = '';
            if ~strcmpi(regexprep(method, '^amva\.', ''), 'mapqn')
                return
            end
            reason = mva_mapqn_reason(sn);
            bool = isempty(reason);
        end

        function [bool, reason] = supportsMarie(sn, method)
            % [BOOL, REASON] = SUPPORTSMARIE(SN, METHOD)
            % The structural premise of 'marie' (closed model, classes routed
            % alike, FCFS/PS/LCFS-PR/delay stations, a multiserver station only
            % with one chain), judged by mva_marie_reason, the predicate the
            % analyzer reads as well. listValidMethods withholds the name on
            % the two halves that make the whole family pointless (open model,
            % class-dependent routing); the rest is reported here.
            bool = true;
            reason = '';
            if ~strcmpi(regexprep(method, '^amva\.', ''), 'marie')
                return
            end
            reason = mva_marie_reason(sn);
            bool = isempty(reason);
        end

        function [bool, reason] = supportsLcfs(sn, method)
            % [BOOL, REASON] = SUPPORTSLCFS(SN, METHOD)
            % A non-preemptive LCFS station is served by one arm alone, the
            % two-station LCFS/LCFS-PR recursion (pfqn_lcfsqn_mva) that
            % solver_mva reaches under 'exact', 'mva' and the default ladder.
            % No AMVA kernel has an LCFS arm: solver_amvald_forward left the
            % station's wait at zero and reported the table. The pairing, the
            % station and server counts and the self-loop rule have no registry
            % name. SOLVER_MVA and MVADISPATCH call this same predicate.
            bool = true;
            reason = '';
            lcfs = find(sn.sched == SchedStrategy.LCFS);
            if isempty(lcfs)
                return
            end
            base = regexprep(method, '^amva\.', '');
            bool = false;
            if ~any(strcmp(base, {'default','exact','mva'}))
                reason = sprintf(['an LCFS station is served by the exact LCFS/LCFS-PR recursion alone, ' ...
                    'which the ''%s'' method does not reach; use ''default'', ''exact'' or ''mva'''], base);
                return
            end
            lcfspr = find(sn.sched == SchedStrategy.LCFSPR);
            if numel(lcfs) ~= 1 || numel(lcfspr) ~= 1 || sn.nstations ~= 2
                reason = 'LCFS MVA requires exactly one LCFS and one LCFS-PR station and no other station.';
                return
            end
            if any(isinf(sn.njobs))
                reason = 'LCFS MVA requires a closed queueing network.';
                return
            end
            if any(sn.nservers([lcfs, lcfspr]) ~= 1)
                reason = 'LCFS MVA is a single-server recursion at both stations.';
                return
            end
            if ~isempty(sn.lldscaling) || ~isempty(sn.cdscaling) || ~isempty(sn.jdscaling)
                reason = 'LCFS MVA carries no load-, class- or joint-dependent scaling.';
                return
            end
            if strcmp(base, 'default') && ~sn_has_product_form(sn)
                % the ladder hands the model to solver_mva, which refuses an
                % implicit request outside product form; 'mva' is the override
                reason = ['the LCFS recursion is the exact one and this model has no product-form ' ...
                    'solution; use ''mva'' to run it as an approximation'];
                return
            end
            K = sn.nclasses;
            for ist = [lcfs, lcfspr]
                isf = sn.stationToStateful(ist);
                for r = 1:K
                    if sn.rt((isf-1)*K+r, (isf-1)*K+r) > 0
                        reason = 'LCFS MVA does not support self-loops at stations.';
                        return
                    end
                end
            end
            bool = true;
        end

        function [bool, reason] = supportsSjn(sn, method)
            % [BOOL, REASON] = SUPPORTSSJN(SN, METHOD)
            % A shortest-job-next (SJF) station is served by
            % solver_mva_sjn_analyzer alone: Kant's conditional waiting-time
            % equation is a population recursion over chain demands and SCVs,
            % so the model is closed, every station single-server with a BCMP
            % discipline, and only the names the analyzer dispatches may ask
            % for it (any other name reached it too and got the default choice
            % under its own label). Counts have no registry name. MVADISPATCH
            % and the analyzer call this same predicate.
            bool = true;
            reason = '';
            if ~any(sn.sched == SchedStrategy.SJF)
                return
            end
            base = regexprep(method, '^amva\.', '');
            bool = false;
            if any(isinf(sn.njobs))
                reason = ['SolverMVA supports shortest-job-next (SJF) scheduling only in closed models, ' ...
                    'the conditional waiting time equation being a population recursion. Use SolverLDES, ' ...
                    'or SolverMVA with SRPT or PSJF for the preemptive size-based open queue.'];
                return
            end
            if ~any(strcmp(base, {'default','exact','mva','amva','sjn.mva','sjn.amva'}))
                reason = sprintf(['a closed model with an SJF station is served by the SJN analyzer alone, ' ...
                    'which answers ''default'', ''exact'', ''mva'', ''amva'', ''sjn.mva'' and ''sjn.amva''; ' ...
                    'the ''%s'' method has no arm for it'], base);
                return
            end
            if ~isempty(sn.lldscaling) || ~isempty(sn.cdscaling) || ~isempty(sn.jdscaling)
                reason = ['the SJN analyzer reads mean demands and SCVs only and carries no ' ...
                    'load-, class- or joint-dependent scaling.'];
                return
            end
            if any(sn.nodetype == NodeType.Cache) || sn_has_fork_join(sn)
                reason = 'the SJN analyzer has no cache or fork-join arm; use SolverLDES for this model.';
                return
            end
            for ist = 1:sn.nstations
                switch sn.sched(ist)
                    case {SchedStrategy.EXT, SchedStrategy.INF}
                        % no-op
                    case SchedStrategy.SJF
                        if sn.nservers(ist) ~= 1
                            reason = sprintf(['SJN scheduling at station %d requires a single server, ' ...
                                'the response time equation is a single-server one.'], ist);
                            return
                        end
                    case {SchedStrategy.PS, SchedStrategy.LCFSPR, SchedStrategy.FCFS, SchedStrategy.SIRO}
                        if sn.nservers(ist) ~= 1
                            reason = sprintf(['station %d has %d servers, the SJN analyzer solves the ' ...
                                'remaining stations with the single-server MVA equation.'], ist, sn.nservers(ist));
                            return
                        end
                    otherwise
                        reason = sprintf('The SJN analyzer does not support %s scheduling at the other stations.', ...
                            SchedStrategy.toText(sn.sched(ist)));
                        return
                end
            end
            bool = true;
        end

        function [bool, reason] = supportsOrderIndependent(sn, method)
            % [BOOL, REASON] = SUPPORTSORDERINDEPENDENT(SN, METHOD)
            % An order-independent (OI) or pass-and-swap (PAS) station is
            % served by solver_mva_oi_analyzer alone, the exact CMVA that
            % mvaDispatch reaches under 'default' and 'exact'. Its premise (a
            % closed model, a non-empty all-zero swap graph and a service-rate
            % function at some OI/PAS station, product-form companions, one
            % class per chain) has no registry name. MVADISPATCH calls this
            % same predicate.
            bool = true;
            reason = '';
            if ~any(sn.sched == SchedStrategy.OI | sn.sched == SchedStrategy.PAS)
                return
            end
            bool = false;
            base = regexprep(method, '^amva\.', '');
            if ~any(strcmp(base, {'default','exact'}))
                reason = sprintf(['SolverMVA supports order-independent (OI) and pass-and-swap (PAS) stations only\n' ...
                    'through its exact order-independent analyzer, which requires method ''default'' or ''exact''\n' ...
                    '(got ''%s''), an empty/zero swap graph at every OI/PAS station, a closed model, and every\n' ...
                    'other station to be product-form (INF, PS, LCFS-PR, SIRO, or class-independent-rate FCFS).\n' ...
                    'Use SolverCTMC or SolverLDES for this model.'], base);
                return
            end
            hasRateFun = false;
            if ~any(isinf(sn.njobs))
                for ist = 1:sn.nstations
                    if sn.sched(ist) ~= SchedStrategy.PAS && sn.sched(ist) ~= SchedStrategy.OI
                        continue
                    end
                    ind = sn.stationToNode(ist);
                    if ind < 1 || ind > numel(sn.nodeparam) || ~isstruct(sn.nodeparam{ind}) ...
                            || ~isfield(sn.nodeparam{ind}, 'swapGraph') ...
                            || ~isfield(sn.nodeparam{ind}, 'svcRateFun') || isempty(sn.nodeparam{ind}.svcRateFun)
                        continue
                    end
                    sg = sn.nodeparam{ind}.swapGraph;
                    if isempty(sg) || any(sg(:) ~= 0)
                        continue
                    end
                    hasRateFun = true;
                    break
                end
            end
            if ~hasRateFun || ~nc_is_oi_model(sn)
                reason = ['the order-independent analyzer requires a closed model, a zero swap graph and a ' ...
                    'service-rate function at every OI/PAS station, and every other station to be product-form ' ...
                    '(INF, PS, LCFS-PR, SIRO, or class-independent-rate FCFS). Use SolverCTMC or SolverLDES for this model.'];
                return
            end
            for c = 1:sn.nchains
                if numel(sn.inchain{c}) > 1
                    reason = 'the order-independent analyzer requires one class per chain (no class switching).';
                    return
                end
            end
            bool = true;
        end

        function [bool, reason] = supportsPolling(sn, method)
            % [BOOL, REASON] = SUPPORTSPOLLING(SN, METHOD)
            % A POLLING station is served by solver_mva_polling_analyzer alone,
            % the station-time analysis of an open multiclass single-station
            % system (Source -> Queue -> Sink) that mvaDispatch claims by shape;
            % anywhere else the station reached solver_amvald_forward, which has
            % no polling arm and left its wait at zero. The analyzer answers
            % 'default' and 'exact' (exhaustive or gated, Poisson arrivals, one
            % server) and refuses every other name. MVADISPATCH and the analyzer
            % call this same predicate.
            bool = true;
            reason = '';
            if ~any(sn.sched == SchedStrategy.POLLING)
                return
            end
            bool = false;
            base = regexprep(method, '^amva\.', '');
            if ~(sn.nclosedjobs == 0 && sn.nclasses > 1 && numel(sn.nodetype) == 3 && ...
                    all(sort(sn.nodetype(:))' == sort([NodeType.Source, NodeType.Queue, NodeType.Sink])))
                reason = ['SolverMVA serves a POLLING station only as an open multiclass single-station ' ...
                    'system (Source -> Queue -> Sink); use SolverLDES or SolverJMT for this model.'];
                return
            end
            if ~any(strcmp(base, {'default','exact'}))
                reason = sprintf(['the polling analyzer answers ''default'' and ''exact'' only; ' ...
                    'the ''%s'' method has no polling arm'], base);
                return
            end
            queue_ist = sn.nodeToStation(sn.nodetype == NodeType.Queue);
            source_ist = sn.nodeToStation(sn.nodetype == NodeType.Source);
            ind = sn.stationToNode(queue_ist);
            ptype = NaN(1, sn.nclasses);
            for r = 1:sn.nclasses
                if ~isnan(sn.rates(queue_ist, r)) && iscell(sn.nodeparam{ind}) && numel(sn.nodeparam{ind}) >= r ...
                        && isfield(sn.nodeparam{ind}{r}, 'pollingType')
                    ptype(r) = PollingType.toId(sn.nodeparam{ind}{r}.pollingType);
                end
            end
            ptype = ptype(~isnan(ptype));
            if isempty(ptype) || any(ptype ~= ptype(1))
                reason = 'the polling analyzer requires one polling type, declared for every class the station serves.';
                return
            end
            if ptype(1) == PollingType.KLIMITED && isfield(sn.nodeparam{ind}{1}, 'pollingPar') ...
                    && ~isempty(sn.nodeparam{ind}{1}.pollingPar) && sn.nodeparam{ind}{1}.pollingPar ~= 1
                reason = 'MVA method unavailable for K-limited polling with K>1.';
                return
            end
            if strcmp(base, 'exact')
                ca = sqrt(sn.scv(source_ist, :));
                if ~(all(ca == 1) && sn.nservers(queue_ist) == 1 && ...
                        (ptype(1) == PollingType.EXHAUSTIVE || ptype(1) == PollingType.GATED))
                    reason = 'MVA exact method unavailable for this model.';
                    return
                end
            end
            bool = true;
        end

        function [bool, reason] = supportsSizeBased(sn, method) %#ok<INUSD>
            % [BOOL, REASON] = SUPPORTSSIZEBASED(SN, METHOD)
            % SRPT, PSJF, FB, LRPT and SETF are served by
            % solver_mva_qsys_sizebased_analyzer alone, the M/GI/1 formulas of
            % Wierman and Harchol-Balter for an open single-station system
            % (Source -> Queue -> Sink) with Poisson arrivals and one server,
            % which mvaDispatch claims by shape under every method name.
            % Anywhere else such a station reached solver_amvald_forward, which
            % has no size-based arm and left its wait at zero. MVADISPATCH calls
            % this same predicate.
            bool = true;
            reason = '';
            sizeBased = [SchedStrategy.SRPT, SchedStrategy.PSJF, SchedStrategy.FB, SchedStrategy.LRPT, SchedStrategy.SETF];
            if ~any(ismember(sn.sched, sizeBased))
                return
            end
            bool = false;
            if ~(sn.nclosedjobs == 0 && numel(sn.nodetype) == 3 && ...
                    all(sort(sn.nodetype(:))' == sort([NodeType.Source, NodeType.Queue, NodeType.Sink])))
                reason = ['SolverMVA serves the size-based disciplines (SRPT, PSJF, FB, LRPT, SETF) only as ' ...
                    'an open single-station M/G/1 system (Source -> Queue -> Sink); use SolverLDES or SolverJMT ' ...
                    'for this model.'];
                return
            end
            queue_ist = sn.nodeToStation(sn.nodetype == NodeType.Queue);
            source_ist = sn.nodeToStation(sn.nodetype == NodeType.Source);
            if sn.nservers(queue_ist) ~= 1
                reason = ['the size-based M/G/1 formulas are single-server ones; use SolverLDES or SolverJMT ' ...
                    'for a multiserver station.'];
                return
            end
            ca2 = sn.scv(source_ist, :);
            if ~all(abs(ca2 - 1) < 1e-6 | ~isfinite(ca2))
                reason = ['the size-based M/G/1 formulas need Poisson arrivals; use SolverLDES or SolverJMT ' ...
                    'for a non-Poisson source.'];
                return
            end
            bool = true;
        end

        function [bool, reason] = supportsQueueingSystem(sn, method)
            % [BOOL, REASON] = SUPPORTSQUEUEINGSYSTEM(SN, METHOD)
            % The by-name premises of the single-station closed forms that
            % solver_mva_qsys_analyzer answers on the single-class open
            % Source -> Queue -> Sink shape. Which laws a model uses is a
            % registry matter (getMethodFeatureSet keeps 'mm1', 'mmk', 'qed'
            % and 'erlanga' exponential), but the registry cannot tell the
            % ARRIVAL law from the SERVICE law, nor see a server count: the
            % M/G/1 names read a Poisson source, the G/M/1 names an exponential
            % server, every G/G/1 and M/G/1 name one server, and 'exact' has a
            % closed form only where one of the two laws is exponential. Each
            % of these answered the other system under the caller's name. The
            % analyzer calls this same predicate.
            bool = true;
            reason = '';
            if ~(sn.nclasses == 1 && sn.nclosedjobs == 0 && numel(sn.nodetype) == 3 && ...
                    all(sort(sn.nodetype(:))' == sort([NodeType.Source, NodeType.Queue, NodeType.Sink])))
                return
            end
            if ~qsys_serves_method(method)
                return
            end
            queue_ist = sn.nodeToStation(sn.nodetype == NodeType.Queue);
            source_ist = sn.nodeToStation(sn.nodetype == NodeType.Source);
            k = sn.nservers(queue_ist);
            ca = sqrt(sn.scv(source_ist, 1));
            cs = sqrt(sn.scv(queue_ist, 1));
            arrivalExp = sn.procid(source_ist, 1) == ProcessType.EXP;
            serviceExp = sn.procid(queue_ist, 1) == ProcessType.EXP;
            singleServer = {'mm1','mg1','mgi1','gm1','gim1','gig1','gig1.kingman','gig1.gelenbe', ...
                'gig1.heyman','gig1.kimura','gig1.allen','gig1.kobayashi','gig1.klb','gig1.marchal','gig1.extremal'};
            bool = false;
            if any(strcmp(method, singleServer)) && k ~= 1
                reason = sprintf(['the ''%s'' closed form is a single-server one and the station has %d servers; ' ...
                    'use ''mmk'', ''gigk'' or one of the ''gigk.*'' methods'], method, k);
                return
            end
            if any(strcmp(method, {'mg1','mgi1','mgisrgi'})) && ~arrivalExp
                reason = sprintf(['the ''%s'' closed form needs a Poisson source and the arrival process is %s; ' ...
                    'use ''gig1'' (or ''gigk.diffusion'')'], method, ProcessType.toText(sn.procid(source_ist, 1)));
                return
            end
            if any(strcmp(method, {'gm1','gim1'})) && ~serviceExp
                reason = sprintf(['the ''%s'' closed form needs exponential service and the service process is %s; ' ...
                    'use ''gig1'''], method, ProcessType.toText(sn.procid(queue_ist, 1)));
                return
            end
            if any(strcmp(method, {'erlanga','mgisrgi'})) && isempty(sn_patience_handles(sn, queue_ist, 1))
                reason = sprintf('method ''%s'' needs a reneging patience law on the queue.', method);
                return
            end
            % 'exact' resolves to M/M/1, M/M/k, M/G/1 or G/M/1 (or the MX/M/1
            % branch on a batch Poisson source) and has nothing for a G/G/k, or
            % a G/G/1 with both laws non-exponential.
            batchPoisson = sn.procid(source_ist, 1) == ProcessType.BMAP && serviceExp;
            if strcmp(method, 'exact') && ~batchPoisson && ...
                    ~((ca == 1 && cs == 1) || (ca == 1 && k == 1) || (cs == 1 && k == 1))
                reason = ['MVA exact method unavailable for this model: no exact closed form for a G/G/k, ' ...
                    'or a G/G/1 with both laws non-exponential; use ''default'''];
                return
            end
            bool = true;
        end

        function [bool, reason] = supportsPreemptivePriority(sn, method) %#ok<INUSD>
            % [BOOL, REASON] = SUPPORTSPREEMPTIVEPRIORITY(SN, METHOD)
            % The preemptive-resume priority arm of solver_amvald_forward
            % (Chandy-Lakshmi; 'priomva' and every AMVA name that reaches the
            % forward step) is a single-server one, and the server count has no
            % registry name. SOLVER_AMVA calls this same predicate before the
            % forward step, which used to be where the refusal surfaced.
            bool = true;
            reason = '';
            k = find(sn.sched(:) == SchedStrategy.FCFSPRPRIO & isfinite(sn.nservers(:)) & sn.nservers(:) > 1, 1);
            if isempty(k)
                return
            end
            bool = false;
            reason = sprintf(['Station %d uses FCFSPRPRIO with %d servers. The preemptive-resume priority arm ' ...
                '(priomva) is implemented for single-server stations only; use SolverCTMC or SolverSSA for ' ...
                'multiserver PRS.'], k, sn.nservers(k));
        end

        function [bool, reason] = supportsExactness(model, method)
            % [BOOL, REASON] = SUPPORTSEXACTNESS(MODEL, METHOD)
            % Method 'exact' requires a product-form solution, the same rule
            % mvaDispatch enforces at solve time. Order-independent and
            % pass-and-swap stations are exempt: solver_mva_oi_analyzer is
            % exact for them regardless of the product-form test. Single-station
            % open systems are exempt too: they go to a queueing-system formula
            % (M/G/1 PK, M/M/k, Cobham, matrix-geometric, ...) that holds
            % outside product form, never to the MVA recursion. Product form
            % has no registry feature name, so the check cannot live in
            % getMethodFeatureSet and belongs here.
            bool = true;
            reason = '';
            if ~strcmp(method, 'exact')
                return
            end
            if model.hasProductFormSolution()
                return
            end
            sn = model.getStruct();
            if any(sn.sched == SchedStrategy.OI | sn.sched == SchedStrategy.PAS)
                return
            end
            isSingleStation = sn.nclosedjobs == 0 && length(sn.nodetype)==3;
            if isSingleStation && all(sort(sn.nodetype)' == sort([NodeType.Source,NodeType.Cache,NodeType.Sink]))
                return
            end
            if isSingleStation && all(sort(sn.nodetype)' == sort([NodeType.Source,NodeType.Queue,NodeType.Sink]))
                % Exempt only where mvaDispatch has a closed form for the
                % shape: the single-class analyzer (SUPPORTSQUEUEINGSYSTEM
                % decides which), the size-based M/G/1 formulas, the polling
                % station-time analysis (SUPPORTSPOLLING decides), Cobham's HOL
                % M/G/1 and the truncated M/M/1-DPS chain, each on the model
                % its branch tests for. Any other multiclass single station
                % falls to solver_mva, which refuses 'exact' outside product form.
                queue_ist = sn.nodeToStation(sn.nodetype == NodeType.Queue);
                source_ist = sn.nodeToStation(sn.nodetype == NodeType.Source);
                qsched = sn.sched(queue_ist);
                scva = sn.scv(source_ist, :);
                scvs = sn.scv(queue_ist, :);
                poissonArrivals = all(abs(scva - 1) < 1e-6 | ~isfinite(scva));
                expService = all(abs(scvs - 1) < 1e-6 | ~isfinite(scvs));
                sizeBased = [SchedStrategy.SRPT, SchedStrategy.PSJF, SchedStrategy.FB, SchedStrategy.LRPT, SchedStrategy.SETF];
                if sn.nclasses == 1 || any(qsched == sizeBased) || qsched == SchedStrategy.POLLING
                    return
                end
                if qsched == SchedStrategy.HOL && sn.nservers(queue_ist) == 1 && poissonArrivals
                    return
                end
                if qsched == SchedStrategy.DPS && sn.nclasses <= 3 && sn.nservers(queue_ist) == 1 ...
                        && poissonArrivals && expService
                    return
                end
                bool = false;
                reason = ['method ''exact'' has no closed form for this multiclass single-station system ' ...
                    'outside product form; use ''default'''];
                return
            end
            bool = false;
            reason = ['method ''exact'' requires a product-form solution; ' ...
                'use ''mva'' for the approximation based on the exact MVA algorithm'];
        end

        function [bool, reason] = supportsFiniteCapacity(model, method)
            % [BOOL, REASON] = SUPPORTSFINITECAPACITY(MODEL, METHOD)
            % MVA-specific finite-capacity gate: a Blocking-After-Service model
            % is exempt where the finite buffers ARE honoured, which is
            % solver_sqd alone: method 'sqd' names it and the default ladder of
            % solver_mva_analyzer selects it. Every other name (the exact
            % recursion under 'mva', every AMVA kernel, sum, mvac) solves the
            % buffers away, and so does 'default' on a load-dependent model,
            % whose ladder (solver_mvald_analyzer) has no sqd arm. Everything
            % else defers to the shared product-form gate. METHOD '' is the
            % static solver-level question and is answered as 'default'.
            if nargin < 2 || isempty(method)
                method = 'default';
            end
            bool = true;
            reason = '';
            if ~isa(model, 'Network')
                return
            end
            snbas = model.getStruct();
            if sn_is_bas_model(snbas) && any(strcmp(method, {'default','sqd'})) ...
                    && isempty(snbas.lldscaling) && isempty(snbas.cdscaling) && isempty(snbas.jdscaling)
                return
            end
            % Single-station M/M/1/K with tail drop is handled by the
            % moment-based (MacGregor Smith) qsys_mg1k_loss_mgs branch in
            % solver_mva_qsys_analyzer, exact only at scv=1. It is an
            % approximation in general, so method='exact' is NOT exempted (it
            % must reject); every other method is.
            if ~strcmp(method, 'exact') && sn_is_mm1k_loss(model.getStruct())
                return
            end
            [bool, reason] = NetworkSolver.checkBindingCapacity(model, 'SolverMVA');
        end

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()

            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink','Source',...
                'ClassSwitch','Delay','DelayStation','Queue',...
                'APH','Coxian','Erlang','Exp','HyperExp','BMAP',...
                'Pareto','Weibull','Lognormal','Uniform','Det', ...
                'StatelessClassSwitcher','InfiniteServer','SharedServer','Buffer','Dispatcher',...
                'CacheClassSwitcher','Cache', ...
                'CacheRetrieval', ...
                'Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS','SchedStrategy_FCFSPRPRIO',...
                'SchedStrategy_DPS','SchedStrategy_FCFS','SchedStrategy_SIRO','SchedStrategy_HOL',...
                'SchedStrategy_LCFS','SchedStrategy_LCFSPR','SchedStrategy_POLLING',...
                ...% size-based M/G/1 disciplines, served by solver_mva_qsys_sizebased_analyzer
                ...% (Wierman and Harchol-Balter, SIGMETRICS 2003). Listed only now that the
                ...% analyzer runs: it read sn.visits by STATION index where sn.visits is
                ...% indexed by CHAIN, so it rejected every multiclass model it exists for.
                'SchedStrategy_SRPT','SchedStrategy_PSJF','SchedStrategy_FB',...
                'SchedStrategy_LRPT','SchedStrategy_SETF',...
                'SchedStrategy_OI','SchedStrategy_PAS',...  % exact order-independent path only (solver_mva_oi_analyzer)
                'SchedStrategy_SJF',...  % closed models only (solver_mva_sjn_analyzer)
                'Fork','Forker','Join','Joiner',...
                'JoinPartial',... % quorum join: the MMT fixed point charges the k-th branch completion (fj_ordstat_exp)
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO', 'ReplacementStrategy_LRU',...
                'ReplacementStrategy_HLRU',...
                'MMAP',...  % marked MAP sources (cache LRU via cache_ttl_lrum_map)
                'ClosedClass','SelfLoopingClass','OpenClass','Replayer',...
                'LoadDependence','ClassDependence','JointDependence',...
                ... % c-server stations: the exact recursion, every AMVA kernel,
                ... % qna/rqt and the M/M/k and G/G/k closed forms carry the
                ... % count; getMethodFeatureSet withdraws it from the
                ... % single-server names. FiniteCapacity is deliberately NOT
                ... % here (see getMethodFeatureSet and supportsFiniteCapacity).
                'MultiServer'});
        end

        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)

            featUsed = model.getUsedLangFeatures();
            featSupported = SolverMVA.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
            if bool
                % Registry inclusion cannot see finite capacity (no feature
                % name); apply the structural gate too, so that the static
                % supports() agrees with the runAnalyzer feature gate.
                [bool, reason] = SolverMVA.supportsFiniteCapacity(model);
                if ~bool
                    line_warning(mfilename, '%s\n', reason);
                end
            end
        end

        function options = defaultOptions
            % OPTIONS = DEFAULTOPTIONS()

            options = SolverOptions('MVA');
        end

        function libs = getLibrariesUsed(sn, options)
            % GETLIBRARIESUSED Get list of external libraries used by MVA solver
            % MVA uses internal algorithms, no external library attribution needed
            libs = {};
        end

    end
end
