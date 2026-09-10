classdef SolverBA < NetworkSolver
    % Bound Analysis solver for closed queueing networks.
    %
    % SolverBA is a dedicated home for asymptotic and hierarchical throughput/
    % queue-length BOUNDS. Unlike SolverMVA (point estimates), each method
    % returns an optimistic or pessimistic bound; the .upper/.lower pair for a
    % family brackets the exact solution. Bounds need only demands (visits x
    % service time) and populations -- no service-time distributions -- so the
    % featset is deliberately narrow (closed product-form-parameterized models).
    %
    % FINITE-BUFFER BLOCKING IS REFUSED, NOT BOUNDED. Needing only demands and a
    % population is the BCMP parameterization, which presumes UNBOUNDED buffers;
    % a buffer that binds couples the station occupancies and the resulting
    % numbers do not bracket the blocked model. runAnalyzer therefore gates on
    % sn_has_blocking and listValidMethods drops every blocking-blind method.
    % The exceptions are qrf.bas*/qrf.rsrd, which carry the blocking tables
    % explicitly. Use SolverMVA method 'sqd' for a point estimate.
    %
    % Method families:
    %   Composite:
    %     auto.* runs every feasible noniterative family and keeps the
    %            tightest upper and the tightest lower; single-class closed
    %            only. Bare 'auto' means 'auto.upper'.
    %   Noniterative (SolverBA-native, solver_ba_analyzer):
    %     aba.*  Asymptotic Bound Analysis (Denning-Buzen)
    %     bjb.*  Balanced Job Bounds (Zahorjan et al.)
    %     pb.*   Proportional Bounds (Eager-Sevcik)
    %     gb.*   Geometric Bounds (Casale-Muntz-Serazzi)
    %     sb.*   Simple bounds from power sums (Harel-Namn-Sturm); NOT a
    %            geometric variant despite sitting beside gb. Rejects delay
    %            stations, as does lr; gb accepts them.
    %     mwba.* Multiclass worst-case balanced bound
    %   Hierarchical / iterative (SolverBA-native, level-parameterized):
    %     pbh    Performance Bound Hierarchy (Eager-Sevcik 1983), level option
    %     sib    Successively Improving Bounds (Srinivasan 1985), level option
    %     cbh    Convolutional Bound Hierarchies (Dowdy et al. 1984)
    %     bjbk / pbk  iterative BJB(k)/PB(k) with delay (Casale et al. 2008)
    %     cub    Composite Upper Bound, multiclass upper-only (Kerola 1986)
    %     looping Eager Looping (1984), the multiclass bracket that seeds the
    %            multiple-class PBH; iterates queue-length lower bounds and
    %            the unaccounted-congestion HEAPS to a fixed point
    %     mbjb   multiclass Balanced Job Bounds lower (Kerola eq. 10; seeds cub)
    %     ssd    multiServer Disaggregation bounds (Suri-Dallery 1986)
    %     ldbcmp LD-BCMP closed-open equivalence bound (Anselmi-Cremonesi 2008)
    %     bpt    Achievable-region LP relaxation (Bertsimas-Paschalidis-
    %            Tsitsiklis 1994). THE ONLY OPEN-NETWORK FAMILY: it lower
    %            bounds the mean response times attainable by ANY non-idling
    %            policy in a multiclass open Markovian network, so it is
    %            refused on the closed models every other family requires.
    %            Lower-only, exponential service, single-server, no delay.
    %     bgt    Piecewise-linear Lyapunov bound (Bertsimas-Gamarnik-
    %            Tsitsiklis 2001), the OPEN-network upper side. Solves the
    %            Down-Meyn global-stability LP; a feasible gamma > 0 both
    %            certifies that every work-conserving policy is stable and
    %            yields a finite (loose) bound on the mean queue lengths.
    %            Upper-only, and needs deterministic non-merging routes.
    %     snc    Stochastic network calculus (Fidler-Rizk 2015), the third
    %            OPEN-network family and the only one whose native object is a
    %            TAIL: MGF arrival/service envelopes are propagated hop by hop
    %            and the delay bound is integrated into a mean. Upper-only,
    %            valid for every work-conserving policy, and it needs a
    %            feed-forward station graph with deterministic routing
    %            downstream of the Source. The quantile itself is reached with
    %            getDelayPerc/getBacklogPerc/getPercTable, not through the
    %            mean columns, which are deliberately loose.
    %     scb    Single-Class Bounds (Dowdy et al. 1992). THE ONLY FAMILY THAT
    %            DOES NOT BRACKET THIS MODEL'S OWN SOLUTION: it brackets the
    %            multiclass system the given single-class model aggregates, so
    %            scb.lower IS the exact single-class throughput. Kept out of
    %            auto.* for that reason. Z=0, single-server.
    %   LP-based reduction bounds (QRF library, solver_ba_qrf_analyzer):
    %     qr / qrf.mmi        Quadratic Reduction Framework (single-class PH)
    %     lr / qrf.mmi.linear Linear Reduction variant
    %     mapamva.*           MAP-AMVA (Casale-Smirni, DSN 2009), the ONLY family
    %                         that consumes the CORRELATION between successive
    %                         services rather than the service mean alone. Its LP
    %                         is written over the exact mean-value balances of a
    %                         closed MAP queueing network in the per-phase
    %                         variables QN(i,k)/UN(i,k), so a workload whose
    %                         burstiness moves the bottleneck between stations is
    %                         bounded rather than averaged into a renewal process.
    %                         'MAP'/'MMPP2' reach this family's feature set, the
    %                         QRF chain's (which carries the same (D0,D1) pair at
    %                         every station) and snc's (source law), and no other.
    %                         Single-class closed, single-server, no delay, FCFS
    %                         at the MAP station, and phases at ONE station
    %                         (permuted to last).
    %     qrf.bethe           Tree-reweighted (Bethe) free entropy, lambda=1/M
    %     qrf.bas.bethe       Bethe free entropy over the BAS polytope
    %     qrf.mem, qrf.mmi.ld, qrf.bas.*, qrf.rsrd  further QRF sub-methods
    %   Petri-net moment-relaxation bounds (solver_ba_spnlp_analyzer):
    %     spnlp.upper / spnlp.lower        Markovian LP over the marking
    %                                      moments (Liu 1998, Table 1).
    %                                      Exponential firing times.
    %     spnlp.op.upper / spnlp.op.lower  the same LP without the
    %                                      second-moment, covariance and
    %                                      Little's-law families, so it needs
    %                                      only a mean firing time and admits
    %                                      any phase-type law. Much looser.
    %   The spnlp family is offered ONLY on a model holding Transition nodes,
    %   and is the only family offered there: every other one is parameterized
    %   by demands and a population, which a marking is not.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = SolverBA(model, varargin)
            % SOLVERBA Create a bound-analysis solver instance
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, SolverBA.defaultOptions));
            self.setLang();
        end

        function sn = getStruct(self)
            % GETSTRUCT Get model data structure for analysis
            sn = self.model.getStruct(false);
        end

        [runtime, analyzer] = runAnalyzer(self, options);

        function bounds = getBounds(self)
            % GETBOUNDS Return the {lower,upper} bracket for the requested
            % bound family as a struct with fields Xlower/Xupper (throughput)
            % and Qlower/Qupper (queue lengths). The current method's family
            % prefix (before the first '.') selects the family; a hierarchical
            % method uses its own level option.
            % One-sided families (cub upper-only, mbjb/ldbcmp lower-only)
            % return NaN on the missing side.
            fam = self.options.method;
            dot = strfind(fam,'.');
            if ~isempty(dot)
                fam = fam(1:dot(1)-1);
            end
            % A blocking-blind family on a blocked model is dropped from
            % listValidMethods, so both sides below would silently come back
            % NaN. Refuse by name instead, with the reason runAnalyzer gives.
            % Judged on the RESOLVED full method, not on FAM: the family prefix
            % of 'qrf.bas' is the bare 'qrf', which is not a blocking bound.
            %
            % runAnalyzer routes 'default'/'auto' to 'qrf.bas' on a blocked
            % model, and getBounds has to say so rather than contradict it --
            % but it still cannot BRACKET, because the analyzer solves qrf.bas
            % in the 'max' direction alone. So the routed case is refused on
            % its own terms: one-sided, take the upper bound from getAvgTable.
            if ba_ignores_blocking(ba_resolve_method(self.options.method)) && ...
                    sn_has_blocking(self.getStruct())
                routed = ba_blocking_default(self.getStruct());
                if ~isempty(routed) && any(strcmp(self.options.method, ...
                        {'default','auto','auto.upper'}))
                    line_error(mfilename, ['''%s'' resolves to ''%s'' on this model, which ' ...
                        'has a binding finite buffer, and that bound is UPPER-only: there is ' ...
                        'no bracket to return. Call getAvgTable/getAvg for the upper bound, ' ...
                        'or SolverMVA with method ''sqd'' for a point estimate.'], ...
                        self.options.method, routed);
                end
                line_error(mfilename, ['Family ''%s'' does not support finite-buffer ' ...
                    'blocking; see the SolverBA method gate. Use ''qrf.bas''/''qrf.rsrd'', ' ...
                    'or SolverMVA with method ''sqd'' for a point estimate.'], fam);
            end
            valid = self.listValidMethods();
            Tl = NaN; Ql = NaN; Tu = NaN; Qu = NaN;
            if any(strcmp([fam,'.lower'], valid))
                [Ql,~,~,Tl] = self.solveSide([fam,'.lower']);
            end
            if any(strcmp([fam,'.upper'], valid))
                [Qu,~,~,Tu] = self.solveSide([fam,'.upper']);
            end
            bounds = struct('Tlower',Tl,'Tupper',Tu,'Qlower',Ql,'Qupper',Qu);
        end

        function [Q,U,R,T] = solveSide(self, method)
            % SOLVESIDE Re-run the solver for one side of the bracket.
            % The re-run instance inherits the CALLER'S FULL OPTION SET (level,
            % verbose, tol, ...) and only overrides the method. Constructing it
            % with just 'method' would silently reset options.level to its
            % default of 2, so a hierarchical family (pbh/cbh/pbk/bjbk/sib)
            % reached through getBounds would never tighten as level is raised.
            s = SolverBA(self.model);
            opts = self.options;
            opts.method = method;
            s.setOptions(opts);
            [Q,U,R,T] = s.getAvg();
        end

        function BoundsTable = getBoundsTable(self, keepDisabled)
            % GETBOUNDSTABLE Table of the {lower,upper} bracket per station and
            % class, in the layout of GETAVGTABLE. Columns:
            %   Station, JobClass, Qlower, Qupper, Tlower, Tupper
            % One-sided families (cub upper-only, mbjb/ldbcmp lower-only) carry
            % NaN on the missing side; NaN is preserved, never replaced by zero.
            if nargin < 2
                keepDisabled = false;
            end
            b = self.getBounds();
            sn = self.getStruct();
            M = sn.nstations;
            K = sn.nclasses;
            Ql = SolverBA.expandBound(b.Qlower, M, K);
            Qu = SolverBA.expandBound(b.Qupper, M, K);
            Tl = SolverBA.expandBound(b.Tlower, M, K);
            Tu = SolverBA.expandBound(b.Tupper, M, K);
            [Qlval, Quval, Tlval, Tuval] = deal([]);
            JobClass = {};
            Station = {};
            for ist = 1:M
                for k = 1:K
                    vals = [Ql(ist,k), Qu(ist,k), Tl(ist,k), Tu(ist,k)];
                    finite = vals(~isnan(vals));
                    % Mirror getAvgTable's drop of disabled station-class pairs,
                    % but NaN-safe: a row is kept when any value that is present
                    % is nonzero, so an all-NaN side never removes the row.
                    if keepDisabled || isempty(finite) || any(finite ~= 0)
                        JobClass{end+1,1} = sn.classnames{k}; %#ok<AGROW>
                        Station{end+1,1} = sn.nodenames{sn.stationToNode(ist)}; %#ok<AGROW>
                        Qlval(end+1) = Ql(ist,k); %#ok<AGROW>
                        Quval(end+1) = Qu(ist,k); %#ok<AGROW>
                        Tlval(end+1) = Tl(ist,k); %#ok<AGROW>
                        Tuval(end+1) = Tu(ist,k); %#ok<AGROW>
                    end
                end
            end
            Station = label(Station);
            JobClass = label(JobClass);
            Qlower = Qlval(:); Qupper = Quval(:);
            Tlower = Tlval(:); Tupper = Tuval(:);
            BoundsTable = Table(Station, JobClass, Qlower, Qupper, Tlower, Tupper);
            BoundsTable = IndexedTable(BoundsTable);
        end

        function D = getDelayPerc(self, eps)
            % GETDELAYPERC Per-station response-time QUANTILE at violation
            % probability EPS, from the stochastic network calculus bound:
            % D(i,r) is the smallest d for which P{D_ir > d} <= eps is
            % certified. Default eps = 1e-3.
            %
            % This is the native output of the 'snc' family, so the accessor
            % runs the SNC analyzer directly whatever options.method says; a
            % station-class pair carrying no traffic stays NaN. Every other
            % family bounds means only and has no counterpart.
            if nargin < 2
                eps = 1e-3;
            end
            [D, ~] = self.sncPerc(eps);
        end

        function B = getBacklogPerc(self, eps)
            % GETBACKLOGPERC Per-station queue-length QUANTILE, in jobs, at
            % violation probability EPS: B(i,r) is the smallest b for which
            % P{Q_ir > b} <= eps is certified. Default eps = 1e-3.
            % The counterpart of getDelayPerc; see it for the conventions.
            if nargin < 2
                eps = 1e-3;
            end
            [~, B] = self.sncPerc(eps);
        end

        function PercTable = getPercTable(self, eps)
            % GETPERCTABLE Table of the response-time and queue-length
            % quantiles at violation probability EPS, in the layout of
            % GETAVGTABLE. Columns: Station, JobClass, RespTPerc, QLenPerc.
            % Rows are the station-class pairs that carry traffic.
            if nargin < 2
                eps = 1e-3;
            end
            [D, B] = self.sncPerc(eps);
            sn = self.getStruct();
            Station = {}; JobClass = {};
            Dval = []; Bval = [];
            for ist = 1:sn.nstations
                for k = 1:sn.nclasses
                    if isnan(D(ist,k)) && isnan(B(ist,k))
                        continue
                    end
                    Station{end+1,1} = sn.nodenames{sn.stationToNode(ist)}; %#ok<AGROW>
                    JobClass{end+1,1} = sn.classnames{k}; %#ok<AGROW>
                    Dval(end+1) = D(ist,k); %#ok<AGROW>
                    Bval(end+1) = B(ist,k); %#ok<AGROW>
                end
            end
            Station = label(Station);
            JobClass = label(JobClass);
            RespTPerc = Dval(:);
            QLenPerc = Bval(:);
            PercTable = Table(Station, JobClass, RespTPerc, QLenPerc);
            PercTable = IndexedTable(PercTable);
        end

        function [D, B] = sncPerc(self, eps)
            % SNCPERC Both quantile matrices from one envelope propagation.
            % The analyzer returns the (arrival, service) envelope handles it
            % built, so the two quantiles cost one pass over the network and
            % one Chernoff search per pair and metric.
            sn = self.getStruct();
            opts = self.options;
            opts.method = 'snc.upper';
            [~,~,~,~,~,~,~,~,env] = solver_ba_snc_analyzer(sn, opts);
            D = NaN(sn.nstations, sn.nclasses);
            B = NaN(sn.nstations, sn.nclasses);
            for ist = 1:sn.nstations
                for k = 1:sn.nclasses
                    if env.carries(ist,k)
                        D(ist,k) = snc_perc_delay(env.arv{ist,k}, env.srv{ist,k}, eps);
                        B(ist,k) = snc_perc_backlog(env.arv{ist,k}, env.srv{ist,k}, eps);
                    end
                end
            end
        end

        function [allMethods] = listValidMethods(self)
            % LISTVALIDMETHODS Bound methods this MODEL can run.
            % A narrowing of listAllMethods: a name that would always be
            % refused on this model is not offered, so a caller enumerating
            % the list never asks for one. Ask for it by name anyway and
            % runAnalyzer still dispatches, so the analyzer's own reason is
            % what comes back.
            allMethods = SolverBA.listAllMethods();

            % The STRUCTURAL premises -- single-class closed, fully closed,
            % single-server, product form, the phase-type chains, blocking, the
            % Petri split and the open-network routing -- come from
            % BA_METHOD_REFUSAL, the same predicate runAnalyzer and the
            % analyzers raise on and supportsModelMethod reports. It is asked
            % first so that every later narrowing works on names this model
            % could actually run. Before it, a two-class closed network was
            % offered all 36 demand-parameterized bounds and 30 of them raised
            % the moment they were run. The blocks below it are PROJECTIONS
            % that only narrow further (they also drop names the feature gate
            % would refuse, e.g. 'lr' on a delay model), kept so the list stays
            % the one the python, JAR and C++ twins produce.
            % ... and BA_METHOD_DEGENERATE withholds the second kind of name:
            % one whose premises this model MEETS but whose formula says nothing
            % here. 'ldbcmp.lower' at N == Qhat is the only such case: it
            % reports the trivial X >= 0, which propagates into an all-zero
            % table a caller cannot tell from an answer. Offering is what stops;
            % asking for it by name still runs, since a vacuous bound is a valid
            % one and the analyzer is entitled to publish it.
            sn = self.model.getStruct();
            keepStruct = false(size(allMethods));
            for mi = 1:numel(allMethods)
                keepStruct(mi) = isempty(ba_method_refusal(sn, allMethods{mi}, self.options)) && ...
                    isempty(ba_method_degenerate(sn, allMethods{mi}));
            end
            allMethods = allMethods(keepStruct);

            % The QR/LR/QRF reduction bounds share one premise: a single-class
            % closed network of single-server stations (solver_ba_qrf_analyzer
            % gates on it, and solver_ba_analyzer refuses 'lr' on a delay).
            % Naming them on a model they cannot run turns a rejection into a
            % method a caller is invited to ask for, the same reason 'sqni' is
            % gated in SolverMVA.
            reducible = ~any(isinf(sn.njobs)) && sn.nclasses == 1 && ...
                ~any(sn.sched == SchedStrategy.INF) && ...
                ~any(sn.nservers(sn.sched ~= SchedStrategy.INF) > 1);
            % The LOAD-DEPENDENT arms survive where the rest of the family
            % cannot run: alpha(i,n) is the rate law of a delay (alpha = n), of
            % a c-server station (alpha = min(n,c)) and of limited load
            % dependence alike, so 'qrf.mmi.ld' and 'qrf.mmi.linear' answer
            % those models on the model's own chain. sn_to_qrf_alpha owns the
            % one restriction that survives, exponential service wherever a
            % station serves several jobs at once. Dropping them with the rest
            % would hide from a caller enumerating the list the only two bound
            % methods such a model has.
            ldReducible = false;
            if ~reducible && ~any(isinf(sn.njobs)) && sn.nclasses == 1
                [~, alphaMsg] = sn_to_qrf_alpha(sn);
                ldReducible = isempty(alphaMsg);
            end
            if ~reducible
                isQrf = strcmp(allMethods,'qr') | strcmp(allMethods,'lr') | ...
                    startsWith(allMethods,'lr.') | startsWith(allMethods,'qrf.');
                isLdArm = strcmp(allMethods,'qrf.mmi.ld') | ...
                    strcmp(allMethods,'qrf.mmi.linear');
                allMethods = allMethods(~isQrf | (ldReducible & isLdArm));
            end

            % 'mapamva' shares the LP families' premise exactly -- single-class
            % closed, single-server, no delay -- so it narrows with them rather
            % than being offered on a model it would refuse on contact. It is
            % NOT one of the load-dependent arms: its q carries no population
            % index, so a delay or a c-server station has nowhere to go.
            if ~reducible
                allMethods = allMethods(~startsWith(allMethods,'mapamva'));
            end

            % 'bpt', 'bgt' and 'snc' are the mirror image: all three are
            % derived for an open network of single-server exponential
            % stations, so every closed model, every delay station and every
            % multiserver station rules them out. Every other family rules OUT
            % the open model, so on an open network the list narrows to those
            % three. 'bgt.upper' additionally needs deterministic non-merging
            % routes and 'snc.upper' a feed-forward station graph, both of
            % which their analyzers check by walking the routing matrix -- too
            % expensive to repeat here, so they stay listed and refuse by name.
            bptOk = all(isinf(sn.njobs)) && ...
                ~any(sn.sched == SchedStrategy.INF) && ...
                ~any(sn.nservers(sn.sched ~= SchedStrategy.INF) > 1);
            if ~bptOk
                allMethods = allMethods(~strcmp(allMethods,'bpt.lower') & ...
                    ~strcmp(allMethods,'bgt.upper') & ...
                    ~strcmp(allMethods,'snc.upper'));
            end
            if all(isinf(sn.njobs))
                keep = strcmp(allMethods,'bpt.lower') | strcmp(allMethods,'bgt.upper') | ...
                    strcmp(allMethods,'snc.upper');
                allMethods = allMethods(keep);
            end

            % 'spnlp.*' is the only family indexed by a MARKING rather than by
            % demands and a population, and the split is total in both
            % directions: on a Petri net nothing else has a representation of
            % the model, and off one spnlp has nothing to read. Two of the
            % gates above already half-cover this by accident -- a Place is an
            % INF station, so 'reducible' and 'bptOk' are both false on any
            % Petri net and the QR/LR/QRF and open families are gone before
            % this point -- but the demand-parameterized families ABA, BJB, PB,
            % GB and the rest survive them and must be dropped by name.
            isPetri = any(sn.nodetype == NodeType.Transition);
            spnlpNames = startsWith(allMethods,'spnlp');
            if isPetri
                allMethods = allMethods(spnlpNames);
            else
                allMethods = allMethods(~spnlpNames);
            end

            % A binding finite buffer rules out everything but the QRF
            % blocking bounds: the other families presume unbounded buffers,
            % and runAnalyzer refuses them by name on such a model. The list
            % can legitimately come back EMPTY here -- a blocked model that is
            % not single-class closed single-server has no bound method at
            % all, and offering one would be the mis-selection this gate
            % exists to prevent.
            if sn_has_blocking(sn)
                keep = false(size(allMethods));
                for mi = 1:numel(allMethods)
                    keep(mi) = ~ba_ignores_blocking(ba_resolve_method(allMethods{mi}));
                end
                allMethods = allMethods(keep);
                % 'default' is offered back when it now MEANS one of the
                % survivors: runAnalyzer routes it to 'qrf.bas' on a blocked
                % model of the right shape, so a caller enumerating the list
                % would otherwise be told the model's own default is invalid.
                if ~isempty(ba_blocking_default(sn))
                    allMethods = [{'default'}, allMethods(:)'];
                end
            end
        end

        function featSupported = getMethodFeatureSet(self, method)
            % FEATSUPPORTED = GETMETHODFEATURESET(METHOD)
            %
            % The base envelope with the per-method deltas the registry CAN
            % name. A feature set says "I accept this construct", so it can
            % refuse a model for HAVING one and never for lacking one; that is
            % exactly the shape of the delay-station and open-class premises
            % below, and exactly not the shape of "one class" or "one server",
            % which have no feature name and live in BA_METHOD_REFUSAL instead.
            %
            % Judged on the RESOLVED name so that 'default' carries the envelope
            % of the gb.upper it runs as.
            %
            % DELAY STATIONS. 'sb' and 'lr' reject an infinite-server station
            % outright, and 'harel', 'sib' and 'scb' reject a nonzero think
            % time, which on these models is the same station: harel
            % extrapolates the exact normalizing constant of a delay-free
            % network, SIB Section 3.2 is the extension that would carry Z and
            % is not implemented, and the SCB Theorem 3 rests on the delay-free
            % balanced-network throughput. The three OPEN families reject one
            % too: the achievable region, the Lyapunov LP and the SNC envelopes
            % are each derived for one server per station.
            %
            % CLASS TYPES. The three OPEN families drop ClosedClass, which is
            % the whole of their class premise. The MIRROR delta -- dropping
            % OpenClass from every demand-parameterized family -- is
            % deliberately NOT applied: 'supports single-class closed networks
            % only' is ONE rule, its single-class half has no feature name, and
            % splitting it across the two mechanisms would report the closed
            % half here and the single-class half in BA_METHOD_REFUSAL for the
            % same model. It is stated once, structurally. 'spnlp' takes no
            % delta at all: it is indexed by the marking, and whether that
            % marking is bounded is a question about the P-invariants of the
            % net, which spn_lpbnd answers.
            featSupported = SolverBA.getFeatureSet();
            if isa(self.model, 'Network')
                % Model-aware: on a blocked model 'default' carries the
                % envelope of the qrf.bas it is routed to.
                resolved = ba_resolve_model_method(self.model.getStruct(), method);
            else
                resolved = ba_resolve_method(method);
            end
            dot = strfind(resolved,'.');
            if isempty(dot)
                fam = resolved;
            else
                fam = resolved(1:dot(1)-1);
            end
            if any(strcmp(fam, {'bpt','bgt','snc'}))
                featSupported.setFalse({'ClosedClass', ...
                    'Delay','DelayStation','SchedStrategy_INF'});
                if ~strcmp(fam, 'snc')
                    % SERVICE AND ARRIVAL LAWS. 'bpt' and 'bgt' are derived for
                    % a MARKOVIAN open network and read the mean alone, so a
                    % non-exponential law anywhere is not refused by them at run
                    % time -- it is silently bounded as if it were Poisson.
                    % Measured on the M/M/1 shape: replacing the Exp(1) source
                    % by an Erlang of the same mean leaves bgt.upper at QLen
                    % 32.6667 and bpt.lower at 1.0, digit for digit. That is a
                    % bound on a DIFFERENT system, so the laws are dropped here
                    % rather than left to a run-time check the analyzers do not
                    % make. Their own procid test covers the queueing stations
                    % only and would miss exactly the source case.
                    %
                    % 'snc' is excluded: it CONSUMES the arrival law (the same
                    % substitution moves it from 3.8244 to 3.0092) and its
                    % analyzer branches on a non-exponential source deliberately.
                    % Its rule is about the SERVICE only, which no feature name
                    % can say, so it lives in BA_METHOD_REFUSAL instead.
                    featSupported.setFalse({'APH','Coxian','Cox2','Erlang', ...
                        'HyperExp','PH','Det','Lognormal','Pareto','Uniform','Weibull'});
                end
            elseif any(strcmp(fam, {'sb','harel','sib','scb','lr'}))
                featSupported.setFalse({'Delay','DelayStation','SchedStrategy_INF'});
            end
            % MODULATED SERVICE, THE ONE DELTA THAT RUNS THE OTHER WAY. 'MAP'
            % and 'MMPP2' sit in the base envelope so that the PHASE-TYPE
            % CHAINS can accept them -- 'mapamva', whose LP is written over the
            % (D0,D1) pair of its MAP station, and 'qr'/'qrf.*', whose q is
            % built from D1 (completions) and D0 (background phase changes) at
            % every station, i.e. the closed MAP queueing network the QRF
            % library is derived for. 'snc' keeps them too, because it CONSUMES
            % the law of a Markovian Source (makeMap in its analyzer) and its
            % service rule is structural (BA_OPEN_REFUSAL). Every OTHER family
            % has to give them back: each reads the service MEAN alone and would
            % bound a correlated model as though its services were independent,
            % quietly returning a bracket for a different system. Written as a
            % delta on the complement rather than as a grant because a feature
            % set can refuse a model for HAVING a construct and never for
            % lacking one.
            %
            % 'mapamva' takes the delay delta instead: its LP is a network of
            % queues, and Casale-Smirni name the delay extension as open work.
            if strcmp(fam, 'mapamva')
                featSupported.setFalse({'Delay','DelayStation','SchedStrategy_INF'});
            elseif ~any(strcmp(fam, {'qrf','snc'}))
                featSupported.setFalse({'MAP','MMPP2'});
            end
            % LOAD DEPENDENCE. The two load-dependent QRF arms read
            % sn.lldscaling through sn_to_qrf_alpha (alpha(i,n) IS the rate law
            % of limited load dependence, composed with the delay and
            % multiserver ones), so they are the two bound methods a model
            % declaring setLoadDependence has. Nothing else here reads it:
            % ldbcmp evaluates its fixed-rate form (Heffes c = 0) from the base
            % rates, so it keeps refusing the feature rather than answering the
            % unscaled network under the caller's label.
            if any(strcmp(resolved, {'qrf.mmi.ld','qrf.mmi.linear'}))
                featSupported.setTrue('LoadDependence');
            end
            % PRIORITY DISCIPLINES. 'mwba' is the one family whose lower bound
            % is discipline-aware: mwrbb_disc_code maps HOL (= FCFSPRIO) to the
            % non-preemptive lemma and the preemptive-resume pair to Lemma 3 of
            % Majumdar-Woodside, so those three are its own. DPS/GPS and the
            % PS-priority variants are NOT granted: the code maps them to the
            % egalitarian-PS lemma, which is not derived for weighted sharing.
            if strcmp(fam, 'mwba')
                featSupported.setTrue({'SchedStrategy_HOL', ...
                    'SchedStrategy_FCFSPRPRIO','SchedStrategy_LCFSPRPRIO'});
            end
            % MULTISERVER (registry name since 2026-09-05) is out of the base
            % envelope: every demand-parameterized family reads one server per
            % station (BA_METHOD_REFUSAL's singleServerFams and fullyClosedFams
            % name 'ssd' as the alternative), the alpha-free QRF arms refuse a
            % c-server station through sn_to_qrf_alpha and the open families
            % through BA_OPEN_REFUSAL. What carries the count is granted here:
            % 'ssd' (the multiserver bound), 'ldbcmp' (its fixed-rate form runs
            % on the c-server rate law), 'auto' (which picks among them) and the
            % two load-dependent QRF arms, whose alpha(i,n) IS min(n,c).
            if any(strcmp(fam, {'ssd','ldbcmp','auto'})) || ...
                    any(strcmp(resolved, {'qrf.mmi.ld','qrf.mmi.linear'}))
                featSupported.setTrue({'MultiServer'});
            end
            % FINITECAPACITY (registry name since 2026-09-05): only the QRF
            % blocking bounds carry the buffer (the MM, MM1, ZZ, ZM, BB, F
            % tables), which is the same split BA_IGNORES_BLOCKING makes; the
            % structural refusal keeps naming them. 'default'/'auto'/'auto.upper'
            % resolve to 'qrf.bas' on a blocked model of the right shape, so
            % they are granted it through RESOLVED. 'spnlp' is NOT granted:
            % its polytope reads no Place capacity (solver_ba_spnlp_analyzer
            % reads sn.nodeparam only), so a capped place would be relaxed away.
            if startsWith(resolved, 'qrf.bas') || startsWith(resolved, 'qrf.rsrd')
                featSupported.setTrue({'FiniteCapacity'});
            end
        end

        function method = resolveMethod(self, options)
            % METHOD = RESOLVEMETHOD(OPTIONS)
            %
            % The method OPTIONS.METHOD runs as on this model, for
            % NetworkSolver.supportsResolvedMethod and hence for findSolver:
            % the aliases of BA_RESOLVE_METHOD, then the finite-buffer routing
            % that sends 'default'/'auto'/'auto.upper' to 'qrf.bas' on a
            % blocked model of the right shape -- the same resolution
            % runAnalyzer dispatches on, so the report gates the name the run
            % serves.
            method = options.method;
            if isa(self.model, 'Network')
                method = ba_resolve_model_method(self.model.getStruct(), method);
            else
                method = ba_resolve_method(method);
            end
        end

        function [bool, reason] = supportsModelMethod(self, method, forRun)
            % [BOOL, REASON] = SUPPORTSMODELMETHOD(METHOD, FORRUN)
            %
            % The feature gate above, plus the structural premises no feature
            % name can express. BA_METHOD_REFUSAL is the same predicate
            % solver_ba_analyzer raises on, so a caller gets one answer
            % whichever of the two it meets first -- which is the point: this
            % gate is what findSolver, listValidMethods, SolverAUTO's ranked
            % choice AND runAnalyzer read (the base envelope can be asked only
            % from inside this method, MATLAB allowing a superclass call from
            % the same-named method alone), and while it was silent about them
            % a two-class closed network was reported as able to run 36 bound
            % methods of which 30 raised on contact.
            %
            % FORRUN (default false) is runAnalyzer's flag: it skips
            % BA_METHOD_DEGENERATE, since a vacuous bound is a valid one that
            % is withheld from the report but published when named.
            %
            % Structural predicate FIRST, then the feature envelope: the same
            % order for the report and the run, so the report's Reason is the
            % sentence the run raises when a model fails both.
            if nargin < 3
                forRun = false;
            end
            if isa(self.model, 'Network')
                structReason = ba_method_refusal(self.model.getStruct(), method, self.options);
                if ~isempty(structReason)
                    bool = false;
                    reason = structReason;
                    return
                end
            end
            [bool, reason] = supportsModelMethod@NetworkSolver(self, method);
            if ~bool || forRun || ~isa(self.model, 'Network')
                return
            end
            % A bound that APPLIES but says nothing is not offered either; see
            % BA_METHOD_DEGENERATE for why that is a separate question and why
            % the analyzer is still allowed to answer it when asked by name.
            degenReason = ba_method_degenerate(self.model.getStruct(), method);
            if ~isempty(degenReason)
                bool = false;
                reason = degenReason;
            end
        end
    end

    methods (Static)
        function [allMethods] = listAllMethods()
            % LISTALLMETHODS Every bound method the solver implements,
            % independently of the model. Hierarchical methods (pbh/sib/cbh/
            % bjbk/pbk/cub/ssd/ldbcmp) are appended here as each backing
            % pfqn_* algorithm is integrated.
            allMethods = { ...
                'default', ...
                'auto.upper','auto.lower', ...
                'aba.upper','aba.lower', ...
                'bjb.upper','bjb.lower', ...
                'pb.upper','pb.lower', ...
                'gb.upper','gb.lower', ...
                'sb.upper','sb.lower', ...
                'harel.upper','harel.lower', ...
                'mwba.upper','mwba.lower', ...
                'pbh.upper','pbh.lower', ...
                'pbk.upper','pbk.lower', ...
                'bjbk.upper','bjbk.lower', ...
                'cbh.upper','cbh.lower', ...
                'ssd.upper','ssd.lower', ...
                'cub.upper','mbjb.lower', ...
                'looping.upper','looping.lower', ...
                'sib.upper','sib.lower', ...
                'scb.upper','scb.lower', ...
                'ldbcmp.lower', ...
                'bpt.lower','bgt.upper','snc.upper', ...
                'qr','lr','lr.upper','lr.lower', ...
                'mapamva.upper','mapamva.lower', ...
                'qrf.mmi','qrf.mem','qrf.bethe','qrf.mmi.ld','qrf.mmi.linear', ...
                'qrf.bas.mmi','qrf.bas.mem','qrf.bas.bethe','qrf.bas','qrf.rsrd', ...
                'spnlp.upper','spnlp.lower','spnlp.op.upper','spnlp.op.lower'};
        end

        function Y = expandBound(X, M, K)
            % EXPANDBOUND Normalize a bracket side to an (M,K) matrix. A
            % one-sided family leaves its missing side as the scalar NaN
            % returned by getBounds; expand it so the table keeps full shape
            % with NaN entries rather than erroring or dropping the column.
            if isempty(X)
                Y = NaN(M,K);
            elseif isscalar(X)
                Y = repmat(X, M, K);
            else
                Y = X;
            end
        end

        function libs = getLibrariesUsed(sn, options) %#ok<INUSL>
            % GETLIBRARIESUSED External libraries used by SolverBA.
            %
            % The Optimization Toolbox is the dependency in every case; which
            % of its solvers is named depends on the method, and getting that
            % wrong is what this used to do. 'qrf.bas', 'qrf.rsrd' and the
            % 'lr' family are LINEAR programs and call linprog ONLY -- not one
            % fmincon call between them -- yet every qrf* name was attributed to
            % fmincon, so a plain SolverBA(model) on a blocked model announced a
            % solver it never invokes.
            %
            % The split is by what the method actually calls: fmincon for the
            % nonlinear arms (qrf.mmi, qrf.mem, qrf.bethe, qrf.mmi.ld,
            % qrf.mmi.linear and
            % the qrf.bas.* entropy objectives), linprog for the LP ones. The
            % nonlinear arms additionally reach linprog and quadprog through
            % qrf_noblo_start's phase-1 feasible start, and 'qrf.mmi.linear'
            % uses quadprog as well -- the toolbox attribution covers those, and
            % naming the primary optimizer is the point here.
            %
            % Judged on the RESOLVED name so the aliases are attributed as what
            % they run: 'qr' is qrf.mmi (fmincon) and bare 'lr' is lr.upper
            % (linprog). The old test read the raw name, so an explicit
            % 'lr.upper' matched neither branch and got no attribution at all.
            libs = {};
            if isempty(options) || ~isfield(options,'method')
                return
            end
            resolved = ba_resolve_method(options.method);
            if startsWith(resolved,'spnlp')
                libs{end+1} = 'MATLAB Optimization Toolbox (linprog)';
            end
            if startsWith(resolved,'qrf') || startsWith(resolved,'lr')
                if any(strcmp(resolved,{'qrf.bas','qrf.rsrd'})) || startsWith(resolved,'lr')
                    libs{end+1} = 'MATLAB Optimization Toolbox (linprog)';
                else
                    libs{end+1} = 'MATLAB Optimization Toolbox (fmincon)';
                end
            end
            if startsWith(resolved,'bpt') || startsWith(resolved,'bgt') || ...
                    startsWith(resolved,'mapamva')
                % mapamva runs mapqn_bnd_lr_mva, a linear program too.
                libs{end+1} = 'MATLAB Optimization Toolbox (linprog)';
            end
        end

        function options = defaultOptions
            % OPTIONS = DEFAULTOPTIONS()
            options = SolverOptions('MVA');
            options.method = 'default';
            % Hierarchy level for pbh/cbh (MVA steps / exact-convolved servers)
            % and iteration count k for pbk/bjbk. Default 2.
            options.level = 2;
        end

        function [bool, featSupported] = supports(model)
            % SUPPORTS Whether SolverBA can bound the given model
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverBA.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()
            %
            % NO SERVICE DISTRIBUTION WAS DECLARED HERE, so SUPPORTS refused
            % EVERY model -- an ordinary closed exponential network came back
            % false with "feature: Exp". Nothing could then reach the solver
            % through a gate: SolverAUTO's ranked choice never proposed
            % SolverBA, and listValidMethods offered no 'ba.*' method name even on the
            % models the bounds are derived for. ('ClosedClass_multiclass', the
            % one name that looked like a distribution gate, is not in
            % SolverFeatureSet.fields at all, so setTrue ignored it silently.)
            %
            % WHAT MAKES A DISTRIBUTION ADMISSIBLE HERE IS ITS MEAN.
            % solver_ba_analyzer reads sn.rates and sn.visits and nothing else:
            % every bound in the ABA/BJB/PB/GB/SB/Harel/MWBA families is a
            % function of the demands D = V./rates and the think time, so any
            % renewal law with a finite mean is admissible whatever its higher
            % moments. solver_ba_qrf_analyzer is the one that needs more, and
            % what it needs is a PH representation -- it reads the {D0,D1} pair
            % out of sn.proc, which is what the phase-type families below carry.
            %
            % THE MODULATED LAWS 'MAP' AND 'MMPP2' ARE IN THE BASE ENVELOPE
            % FOR THE PHASE-TYPE CHAINS ('mapamva' AND 'qr'/'qrf.*') AND FOR
            % 'snc' (WHICH CONSUMES THE SOURCE LAW), AND GETMETHODFEATURESET
            % STRIPS THEM FROM EVERY OTHER. They were out entirely until mapamva landed, on
            % the correct ground that a renewal bound derived for a product-form
            % network says nothing about a correlated one: its mean rate exists,
            % so the utilization law still holds and the formula still returns a
            % number, but that number brackets a DIFFERENT system. MAP-AMVA
            % (Casale-Smirni, DSN 2009) is derived FOR the correlated model --
            % its variables are the per-phase QN(i,k) and UN(i,k) -- so the same
            % reasoning that refuses the others admits it.
            %
            % The direction matters and is forced: a feature set refuses a model
            % for HAVING a construct and never for lacking one, so the only way
            % to grant a law to one family is to put it in the base envelope and
            % take it away from the rest. MMAP and BMAP stay OUT everywhere --
            % the LP is single-class and has no marked or batch arrival variable
            % -- and Cache and Fork/Join stay out for the older reason, that no
            % analyzer here has a representation of either.
            %
            % THE PETRI-NET CONSTRUCTS ARE IN, since solver_ba_spnlp_analyzer.
            % They were out on exactly the same "no representation" ground and
            % that ground is gone: the spnlp relaxation is indexed by the
            % marking, reads the enabling, inhibiting and firing arcs out of
            % sn.nodeparam, and refuses by name the modes it cannot carry
            % (immediate, multi-server, marking-dependent, and phase-type on
            % its Markovian side). QueueingPlace stays OUT and is refused by
            % name: a place with an embedded queue has local state the
            % relaxation has no variable for. Same division SolverNC draws.
            featSupported = SolverFeatureSet;
            featSupported.setTrue({ ...
                'ClassSwitch','Delay','DelayStation','Queue', ...
                'Sink','Source','Router', ...
                'StatelessClassSwitcher', ...
                'ClosedClass','OpenClass', ...
                ...% A self-looping class is a closed chain of one station, which
                ...% the demand-parameterized bounds read as any other chain; the
                ...% qn sanity goldens pin mwba/cub/mbjb/looping rows on the
                ...% *slc* models.
                'SelfLoopingClass', ...
                ...% renewal service laws: the bounds need the mean, the QRF
                ...% reduction needs the PH form, and sn.proc carries both
                'APH','Coxian','Cox2','Erlang','Exp','HyperExp','PH', ...
                'Det','Lognormal','Pareto','Uniform','Weibull', ...
                ...% modulated service, for 'mapamva' alone; see the note above
                'MAP','MMPP2', ...
                'SchedStrategy_INF','SchedStrategy_PS', ...
                'SchedStrategy_FCFS','SchedStrategy_LCFSPR', ...
                'RoutingStrategy_PROB','RoutingStrategy_RAND', ...
                ...% Petri-net constructs, for the spnlp family. QueueingPlace
                ...% is deliberately absent; see the note above.
                'Place','Transition','Linkage','Enabling','Inhibiting', ...
                'Timing','Firing','Storage'});
        end
    end
end
