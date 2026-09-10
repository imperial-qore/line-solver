classdef SolverNC < NetworkSolver
    % SolverNC Normalizing Constant solver for product-form networks
    %
    % SolverNC implements normalizing constant algorithms for analyzing closed
    % product-form queueing networks. It computes the normalizing constant and
    % associated performance measures efficiently without explicitly enumerating
    % all network states, making it suitable for medium to large closed networks.
    %
    % @brief Normalizing constant solver for efficient closed network analysis
    %
    % Key characteristics:
    % - Normalizing constant computation for product-form networks
    % - Avoids explicit state enumeration
    % - Efficient algorithms for closed networks
    % - Multiple computational methods (exact, approximation)
    % - State probability computation via normalization
    %
    % NC solver methods:
    % - Exact normalizing constant computation
    % - IMCI (Improved Modular Computer Implementation)
    % - Linearizer methods (LS, LE)
    % - Interpolation methods (MMINT2, GLEINT)
    % - Approximation methods (CA, Panacea)
    %
    % SolverNC is ideal for:
    % - Closed product-form networks
    % - Medium to large population networks
    % - Systems requiring efficient exact solutions
    % - Networks with complex routing patterns
    % - Performance analysis requiring state probabilities
    %
    % Example:
    % @code
    % solver = SolverNC(model, 'method', 'exact');
    % solver.getProbAggr();          % State probabilities
    % solver.getNormalizingConstant(); % Normalizing constant
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties (Access = private)
        % TRUE only while @SolverNC/runAnalyzer is inside the shared
        % NETWORKSOLVER.RUNANALYZERCHECKS gate, which is a RUN-path gate that
        % reaches SUPPORTSMODELMETHOD through the same two-argument call
        % MODEL.FINDSOLVER uses. The two callers need different answers for
        % 'mmint2'/'gleint' -- the report must not offer a pair that returns a
        % zero table, the run must still perform it, see NC_METHOD_REFUSAL --
        % and the base gate gives no other way to tell them apart.
        gateAsksTheRunQuestion = false;
    end

    methods
        function self = SolverNC(model,varargin)
            % SOLVERNC Create a Normalizing Constant solver instance
            %
            % @brief Creates an NC solver for product-form network analysis
            % @param model Network model to be analyzed via normalizing constant methods
            % @param varargin Optional parameters (method, tolerance, etc.)
            % @return self SolverNC instance configured for NC analysis
            
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            self.setLang();
        end
        
        runtime = runAnalyzer(self, options)
        Pnir = getProb(self, node, state)
        Pnir = getProbAggr(self, node, state_a)
        Pn   = getProbSys(self)
        Pn   = getProbSysAggr(self)
        RD = getCdfRespT(self, R);
        
        function [normConst,lNormConst] = getNormalizingConstant(self)
            normConst = exp(getProbNormConstAggr(self));
            lNormConst = getProbNormConstAggr(self);
        end
        
        [lNormConst] = getProbNormConstAggr(self)

        [b, lG, lH] = getAvgBusyPeriod(self, stations, n)
        
        function sn = getStruct(self)
            % QN = GETSTRUCT()
            
            % Get data structure summarizing the model
            sn = self.model.getStruct(false); %no need for initial state
        end
        
        function tf = supportsExactSensitivity(self) %#ok<MANU>
            % TF = SUPPORTSEXACTSENSITIVITY()
            % The normalizing-constant solver is exact on the same
            % product-form class that pfqn_sens differentiates, so
            % getSensitivityTable uses the analytic branch.
            tf = true;
        end

        function [allMethods] = listValidMethods(self)
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            sn = self.model.getStruct();
            % 'ms' names the Manjunath-Sikdar transform of the loss-network
            % analyzer, which is the only place it is admissible; it is listed
            % because solver_nc_lossn_analyzer branches on the token and would
            % otherwise be unreachable, unlike 'exact' which the product-form
            % gate below refuses on a region model.
            % 'rayint' and 'spm' both name the SPM saddle point on a cache, which
            % serves cache_spm_size once the items carry storage costs. On a
            % retrieval model 'rayint' is instead the ray/WKB delayed-hit
            % expansion, admissible only with an infinite-server fetch system;
            % solver_nc_retrieval_analyzer branches on the method name and warns and
            % falls back to 'exact' anywhere else.
            % 'rec' is the MDD-rec route: on a loss network it is the exact
            % normalizing constant without the residue transform's integrality
            % demand, and on a Petri net it is the only admissible method.
            % 'divdiff' is the divided-difference closed form of Casale
            % (SIGMETRICS 2017); it needs no think time, since a delay would ask
            % for the integral form of Cor. 3.4, and pfqn_nc refuses one by name.
            % Load-dependent rates ARE served: pfqn_ncld substitutes the limited
            % load-dependent kernel of Casale-Harrison-Ong (Perform. Eval. 2021),
            % Thm. 1, and reports itself as 'divdiff.ld/...'.
            allMethods = {'default','exact','divdiff','rayint','spm','ms','erlangfp','mci','imci','ls',...
                'le','ble','aghq','mmint2','gleint','pana','panald','ca',...
                'clw','kt','bkt','lekt','bk','bkue','lc','lc.ue','sampling','is',...
                'mcmc',...
                'propfair','comom','comomld','cub','rgf','ger',...
                'rd', 'nrp','nrl','nre','gm','mem','sdr','sdr.mva',...
                ... % 'morrison' is the heavy-usage asymptotic expansion of the
                ... % generating function for a closed think+DPS network
                ... % (npfqn_dps_morrison, solver_nc_dps_analyzer). It is the
                ... % DEFAULT on that shape and inadmissible anywhere else, where
                ... % runAnalyzer refuses it: nothing else in NC can see the DPS
                ... % weights. Non-product-form, so it returns no lG.
                'morrison',...
                ... % 'rec' is the MDD-rec route for product-form Petri nets; it is
                ... % the only method solver_nc_spn_analyzer serves, and the only one
                ... % admissible on a net (see runAnalyzer).
                'rec'};
        end

        function method = resolveMethod(self, options) %#ok<INUSL>
            % 'default' is NOT resolved ahead of the run. Every route the
            % default dispatch can take (Morrison, the OI convolution, the
            % cache and loss-network families, comom, the cubature) is gated
            % as 'default' itself: the base envelope plus
            % NC_METHOD_REFUSAL('default'), which is what runAnalyzer asks.
            %
            % Until 2026-09-05 this resolved an open non-Markovian model to
            % 'mem', mirroring a dispatch the run had stopped taking on
            % 2026-07-13 (1df9d1df6: mem on explicit request only). The gate
            % then asked mem's questions about a run that took the
            % product-form path: an open FCFS model with a finite buffer and a
            % non-exponential service passed as mem.blocking and was solved
            % buffer-free by solver_nc.
            method = options.method;
        end

        function featSupported = getMethodFeatureSet(self, method)
            % Per-method feature deltas applied to the base NC envelope.
            %
            % Defining this is what lets NetworkSolver.supportsModelMethod name
            % the offending features: with no method feature set it falls back
            % to the coarse supports(model) and returns an empty reason, so the
            % gate could only report "features not supported" without saying
            % which ones.
            %
            % ONLY THE RESTRICTIONS A FEATURE NAME CAN CARRY LIVE HERE. A
            % feature set declares what the method ACCEPTS, so it can refuse a
            % model for HAVING a construct and never for lacking one: "closed
            % population only" and "no think time" are expressible by dropping
            % OpenClass and SchedStrategy_INF, while "requires a cache" or
            % "requires a loss network" are not and belong to
            % NC_METHOD_REFUSAL, which SUPPORTSMODELMETHOD consults next.
            %
            % A non-Network model (e.g. a LayeredNetwork) has no
            % getUsedLangFeatures, so it keeps the coarse path and the
            % structural checks/redirects that operate on such models.
            if ~isa(self.model, 'Network')
                featSupported = [];
                return;
            end
            featSupported = SolverNC.getFeatureSet();
            method = lower(method);
            % A rate lattice (setLoadDependence) sends the model to SOLVER_NCLD
            % on every route, and PFQN_NCLD has an arm for these names only:
            % any other name raises 'Unrecognized method for solving
            % load-dependent models' there. The SDR analyzer folds the lattice
            % into its own alpha_i(n) and is intercepted upstream.
            if ~any(strcmp(method, {'default','exact','is','clw','pana','panald', ...
                    'divdiff','rd','nrp','nrl','nre','comomld','sdr','sdr.mva'}))
                featSupported.setFalse({'LoadDependence'});
            end
            % Class- or joint-dependent rates divert the model, on every route,
            % to SOLVER_NC_CONV, Sauer's multichain convolution, which never
            % reads the method: any other name would be answered by that
            % recursion. It is exact for the product-form beta_{i,r}(n) but a
            % joint eta_i(n) breaks the BCMP recurrence it runs on (see
            % _kb/04-networkstruct.md), so 'exact' keeps the former only.
            if ~any(strcmp(method, {'default','exact'}))
                featSupported.setFalse({'ClassDependence','JointDependence'});
            elseif strcmp(method, 'exact')
                featSupported.setFalse({'JointDependence'});
            end
            switch method
                case 'divdiff'
                    % The divided-difference closed form of Casale (SIGMETRICS
                    % 2017), Eqs. (15)-(16), covers load-independent queues; a
                    % think time would ask for the integral form of Cor. 3.4,
                    % which is not implemented, so PFQN_NC and PFQN_NCLD both
                    % refuse one by name. An infinite server is where a think
                    % time comes from, so the envelope drops it. A c-server
                    % station is a second source of one (Seidmann's surrogate
                    % delay) and has no feature name: NC_METHOD_REFUSAL.
                    featSupported.setFalse({'SchedStrategy_INF'});
                case {'rd','nrp','nrl','nre','comomld','panald'}
                    % The load-dependent normalizing-constant evaluators are
                    % reached by SOLVER_NCLD only on its CLOSED branch, where
                    % PFQN_NCLD reads the method name. An open chain sends the
                    % model to the mixed route (pfqn_mvaldmx), which never reads
                    % it, so every one of these names silently became 'ncldmx'.
                    % A fork-join model is the same case one step later:
                    % fjFixedPoint hands ncDispatch the MMT image, whose
                    % parallelism rides on OPEN auxiliary classes, so the inner
                    % model is mixed (ncldmx) or, off a lattice, reaches
                    % PFQN_NC, which has no arm for these names.
                    featSupported.setFalse({'OpenClass','Fork','Forker','Join','Joiner','JoinPartial'});
                case 'is'
                    % The sample-an-ordering estimator of PFQN_IS integrates
                    % over a closed population simplex; there is no open-class
                    % form of it, and SOLVER_NC_ANALYZER refuses one by name.
                    % Use 'sampling' (pfqn_mci/pfqn_ls) for an open or mixed
                    % model. The fork-join image is mixed (see above), so a
                    % Fork is refused for the same reason.
                    featSupported.setFalse({'OpenClass','Fork','Forker','Join','Joiner','JoinPartial'});
            end
            % MULTISERVER (registry name since 2026-09-05): the divided-
            % difference closed form of 'divdiff' covers load-independent
            % single-server queues, and NC_METHOD_REFUSAL keeps wording why (a
            % c-server station enters the constant as Seidmann's surrogate
            % delay); every other route folds the count into its own kernel.
            if strcmp(method, 'divdiff')
                featSupported.setFalse({'MultiServer'});
            end
            % FINITECAPACITY (registry name since 2026-09-05) is NOT in the base
            % envelope: the product-form routes solve a buffer away, which is
            % what the binding-capacity gate below refuses. Two arms honour
            % one: 'mem' represents it as a GE/GE/c/0;N queue (mem.blocking,
            % SOLVER_NC_MEM_SUPPORTS), and the single-station M/M/1/K with tail
            % drop is solved in closed form (qsys_mm1k_loss) under 'default' and
            % 'exact'. The shape half of each rule stays structural.
            if any(strcmp(method, {'mem','default','exact'}))
                featSupported.setTrue({'FiniteCapacity'});
            end
        end

        function [bool, reason] = supportsModelMethod(self, method)
            % MEM (Kouvatsos maximum entropy) has structural applicability
            % rules beyond a flat feature set (open-only, no class switching,
            % non-priority scheduling); delegate to solver_nc_mem_supports,
            % which returns a precise reason. All other NC methods inherit the
            % coarse product-form feature gate.
            % Discrete-time (slotted) route: a finite buffer on a Bernoulli
            % server is the loss system of Daduna's corollary 2.8, which
            % solver_nc_dt_analyzer solves exactly, so the structural
            % capacity gate below must not fire. nc_is_dt_model performs the
            % real admissibility check and reports a precise reason.
            if self.isSlotted()
                dt = nc_is_dt_model(self.model.getStruct(), self.getOptions);
                bool = ~strcmp(dt.kind, 'none');
                if bool
                    % the shape admits the route; the NAME is a second question
                    % (the route reads none), asked of NC_METHOD_REFUSAL like
                    % every other name rule so the run refuses the same pairs
                    reason = nc_method_refusal(self.model.getStruct(), method, self.getOptions, ...
                        ~self.gateAsksTheRunQuestion);
                    bool = isempty(reason);
                else
                    reason = sprintf('options.config.slotted is set but %s', dt.reason);
                end
                return
            end
            % A Petri net is served only by the MDD-rec route, and none of the
            % gates below -- which are written about stations, capacities and
            % the queueing-network product form -- say anything about a net.
            % SPN_PF decides the product-form class of a net, and does it by
            % name; see solver_nc_spn_analyzer. NC_METHOD_REFUSAL holds the rule
            % itself, so runAnalyzer refuses the same pairs this gate does.
            if isa(self.model, 'Network') && any(self.model.getStruct().nodetype == NodeType.Place)
                reason = nc_method_refusal(self.model.getStruct(), method, self.getOptions, ...
                    ~self.gateAsksTheRunQuestion);
                bool = isempty(reason);
                return
            end
            memblocking = false;
            if strcmp(method, 'mem')
                sn = self.model.getStruct();
                [bool, reason, memblocking] = solver_nc_mem_supports(sn);
                if bool
                    % a mem-admissible model still has to clear the feature gate;
                    % see _kb/06-solver-catalog.md (mem.blocking note)
                    [bool, reason] = supportsModelMethod@NetworkSolver(self, method);
                end
            else
                [bool, reason] = supportsModelMethod@NetworkSolver(self, method);
            end
            % structural finite-capacity gate (no registry feature name); mem is
            % the one exception (it represents the buffer as a GE/GE/c/0;N queue).
            % see _kb/06-solver-catalog.md (finite capacity gate + mem.blocking)
            if bool && isa(self.model, 'Network') && ~memblocking
                % Single-station M/M/1/K with tail drop is handled exactly by the
                % probability-based qsys_mm1k_loss branch in runAnalyzer; exempt
                % it from the product-form capacity gate.
                if ~sn_is_mm1k_loss(self.model.getStruct())
                    [bool, reason] = NetworkSolver.checkBindingCapacity(self.model, 'SolverNC');
                end
            end
            % 'exact', 'is' and 'panald' need a product-form solution, the
            % same rule runAnalyzer enforces at solve time; 'is' on a
            % pass-and-swap model is the exception (pfqn_pas_is). The other
            % methods fall back to Seidmann's comom on a non-product-form
            % model, so they stay admissible. Product form has no registry
            % feature name, so the check cannot live in getMethodFeatureSet.
            % A loss network (open, one DROP region holding a single Delay) IS
            % product form -- the truncated Poisson law the residue transform
            % of SOLVER_NC_LOSSN_ANALYZER evaluates exactly under 'exact' -- but
            % SN_HAS_BLOCKING reads any region as blocking, so the shape is
            % exempted here and at the 'exact' arm of runAnalyzer alike,
            % through the one predicate both ask, NC_IS_LOSSN_MODEL.
            if bool && isa(self.model, 'Network') ...
                    && any(strcmpi(method, {'exact','is','panald'})) ...
                    && ~self.model.hasProductFormSolution() ...
                    && ~nc_is_lossn_model(self.model.getStruct())
                sn = self.model.getStruct();
                % An OI/PAS station carries a rank rate mu(n) of the whole
                % occupancy vector, so HASPRODUCTFORMSOLUTION reads false, but
                % RUNANALYZER routes such a model to SOLVER_NC_OI_ANALYZER, which
                % is EXACT (pfqn_ncoi). Without this exemption the guard refuses
                % 'exact' before the OI routing is ever reached, which is the same
                % exemption MVADISPATCH makes with its hasOIorPAS test.
                hasOIorPAS = any(sn.sched == SchedStrategy.OI | sn.sched == SchedStrategy.PAS);
                % 'is' is admissible on ANY OI/PAS model, not only a P&S one:
                % SOLVER_NC_PAS_IS_ANALYZER samples the rank rate either way,
                % with PFQN_OI_IS for an empty swap graph and PFQN_PAS_IS
                % otherwise. Gating it on nc_is_pas_model refused the sampler on
                % exactly the models PFQN_OI_IS exists to serve.
                if ~(strcmpi(method, 'is') && hasOIorPAS) ...
                        && ~(any(strcmpi(method, {'exact'})) && hasOIorPAS)
                    bool = false;
                    if hasOIorPAS
                        % NEVER recommend comom here. On an OI/PAS station it is
                        % not a coarser answer, it is a WRONG one: comom reads
                        % sn.rates, which holds only the single-job rate mu([r]),
                        % so the rank rate mu(n) is dropped and the model solved
                        % is an ordinary queue. RUNANALYZER refuses it by name for
                        % exactly that reason, so advising it here would send the
                        % caller into a second refusal, or worse past one.
                        reason = sprintf(['method ''%s'' cannot represent the rank rate ' ...
                            'mu(n) of an order-independent station. Use ''default'' or ' ...
                            '''exact'' (pfqn_ncoi, exact), SolverMVA (exact MVA-OI), or ' ...
                            'SolverCTMC'], method);
                    else
                        reason = sprintf(['method ''%s'' requires a product-form solution; ' ...
                            'use SolverCTMC for an exact answer, or ''comom'', which ' ...
                            'applies Seidmann''s approximation, for an approximate one'], method);
                    end
                end
            end
            % The structural per-method rules the feature registry cannot name:
            % which route a method has on THIS model, and whether it exists.
            % NC_METHOD_REFUSAL is the single copy of them, asked here and by
            % runAnalyzer, so a pair this gate offers is a pair the run accepts.
            if bool && isa(self.model, 'Network')
                reason = nc_method_refusal(self.model.getStruct(), method, self.getOptions, ...
                    ~self.gateAsksTheRunQuestion);
                bool = isempty(reason);
            end
        end

        function bool = isSlotted(self)
            % BOOL = ISSLOTTED()
            % True when the caller asked for the discrete-time route with
            % options.config.slotted. Never inferred from the model: a
            % Geometric service time is an ordinary continuous-time model
            % unless the caller declares the slot lattice.
            bool = false;
            opts = self.getOptions;
            if isstruct(opts) && isfield(opts,'config') && isstruct(opts.config) ...
                    && isfield(opts.config,'slotted') && ~isempty(opts.config.slotted)
                bool = logical(opts.config.slotted);
            end
        end

        function bool = isStochasticMethod(self, method) %#ok<INUSL>
            % BOOL = ISSTOCHASTICMETHOD(METHOD)
            % NC is deterministic except for the Monte Carlo integration
            % methods (mci/imci), logistic sampling (ls), the importance
            % sampling method (is), the Chen-O'Cinneide Markov chain Monte
            % Carlo method (mcmc), and the sampling method, whose estimates
            % depend on the random seed. Method names are tokenized so that
            % runtime-resolved names such as 'default/imci' and prefixed names
            % such as 'nc.ls' classify correctly.
            tokens = regexp(lower(method), '[./]', 'split');
            bool = any(ismember(tokens, {'mci','imci','ls','sampling','is','mcmc'}));
        end
    end
    
    methods (Static)
        
        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()
            
            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink','Source',...
                'ClassSwitch','Delay','DelayStation','Queue',...
                'APH','Coxian','Erlang','Det','Exp','HyperExp',...
                ... % Geometric is admitted for the discrete-time route only
                ... % (options.config.slotted, solver_nc_dt_analyzer); on the
                ... % continuous-time routes it is treated by its mean and SCV
                'Geometric',...
                'StatelessClassSwitcher','InfiniteServer',...
                'SharedServer','Buffer','Dispatcher',...
                ... % Finite capacity regions: NC solves the OPEN single-Delay
                ... % loss-network case exactly (Erlang fixed point,
                ... % solver_nc_lossn_analyzer). It cannot do an FCR on queueing
                ... % stations -- a boolean feature cannot express that split, so
                ... % runAnalyzer keeps the imperative check for the queueing case,
                ... % the same pattern as SolverMVA.supportsFiniteCapacity.
                'Region', ...
                'Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS','SchedStrategy_SIRO',...
                'SchedStrategy_LCFS','SchedStrategy_LCFSPR',...
                ... % DPS is served ONLY in Morrison's closed think+DPS shape
                ... % (nc_is_dps_model -> solver_nc_dps_analyzer). A boolean
                ... % feature cannot express that restriction, so runAnalyzer
                ... % keeps an imperative check for every other DPS model, the
                ... % same pattern as 'Region' above.
                'SchedStrategy_DPS',...
                'RoutingStrategy_SDR',...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'SchedStrategy_FCFS','SchedStrategy_OI','SchedStrategy_PAS',...
                'Fork','Join','Forker','Joiner',... % fork-join via the MMT transformation (fjFixedPoint)
                'JoinPartial',... % quorum join: the MMT fixed point charges the k-th branch completion (fj_ordstat_exp)
                'ClosedClass','SelfLoopingClass',...
                'Cache','CacheClassSwitcher','OpenClass',            ...
                'CacheRetrieval', ...
                'CacheItemSize', ... % per-list storage cost caps, exact and sampling methods
                'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO',...
                'ReplacementStrategy_HLRU',...
                ... % Petri nets: the 'rec' route (solver_nc_spn_analyzer) walks the
                ... % reachable set in a decision diagram, so a Place is a token
                ... % container rather than a station with a service process. A
                ... % queueing Place is refused by SPN_PF, which is where the
                ... % product-form class is decided.
                'Place','Transition','Linkage','Enabling','Inhibiting','Timing','Firing','Storage',...
                'LoadDependence','ClassDependence','JointDependence',...
                ... % c-server stations: every route but 'divdiff' carries the
                ... % count (see getMethodFeatureSet). FiniteCapacity is
                ... % deliberately NOT here: mem, default and exact are granted
                ... % it per method, the rest solve a buffer away.
                'MultiServer'});
            %'OpenClass',...
        end
        
        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverNC.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end        
        
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = SolverOptions('NC');
        end

        function libs = getLibrariesUsed(sn, options)
            % GETLIBRARIESUSED Get list of external libraries used by NC solver
            % NC uses internal normalizing constant algorithms, no external libraries needed
            libs = {};
        end
    end
end
