function reason = jmtMethodRefusal(sn, method, options, engine)
% REASON = JMTMETHODREFUSAL(SN, METHOD, OPTIONS, ENGINE)
%
% The structural half of SolverJMT's method gate: the rules that decide
% whether a JMT METHOD can run this model and that no feature name can state.
% Returns '' when the pair is admissible.
%
% ONE PREDICATE, TWO CALLERS. supportsModelMethod asks it, so findSolver and
% SolverAUTO never offer a pair that would die at run time, and the analyzer
% asks it again -- writeJMVA for the JMVA arm, runAnalyzer for 'replication' --
% so a caller naming the method by hand gets the same sentence rather than a
% JMT stack trace. A second copy of either rule is how the gate and the run
% drift into two different answers.
%
% THE RULES. A finite timespan for 'replication'; single-server stations for
% the eight closed-form JMVA algorithms; a multi-chain model for 'jmva.comom',
% whose JMVA engine answers with an unseeded random perturbation of the model;
% immediate feedback (sn.immfeed), which neither JMT document can state; for
% the JMVA document, a fork or join and a
% non-exponential law at a station whose discipline is not insensitive, since
% that document carries a mean demand and a visit count per chain and nothing
% else; and a binding finite buffer, which neither engine carries -- JSIM
% because no JMT drop strategy reproduces LINE's blocking, JMVA because its
% document has no capacity element at all.
%
% ENGINE ('jsim' or 'jmva') OVERRIDES the engine the METHOD name implies, and
% writeJMVA passes 'jmva' because it IS the JMVA writer whoever called it:
% SolverQNS reaches it with its own method names ('default', 'conway', ...), and
% keying the buffer rule on the name handed a QNS run JSIM's verdict, which threw
% before the .jmva file was written and left qnsolver reporting "Cannot open
% input file". Omit it and the name decides; a name that is neither engine's then
% gets no verdict at all, since the solver that owns it carries its own gate.
%
% WHAT IS NOT HERE. Everything a feature name CAN state lives in
% SolverJMT.getMethodFeatureSet instead: the JMVA envelope is narrower than the
% JSIM one (no cache, no fork-join, no Petri net, no finite capacity region, no
% impatience, and only the BCMP disciplines survive a writer that emits a
% station type, a demand and a visit count), and the closed-only JMVA
% algorithms additionally drop OpenClass and LoadDependence.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

reason = '';

% The transient arm integrates over [0,T]: a mean at an unstated horizon is not
% a quantity, and runAnalyzer auto-selects 'replication' only when the horizon
% IS finite. options.timespan is not a model feature, so the featset cannot see
% it and this is the only place the rule can live.
if strcmpi(method, 'replication')
    if nargin < 3 || ~isstruct(options) || ~isfield(options,'timespan') ...
            || numel(options.timespan) < 2 || ~isfinite(options.timespan(2))
        reason = ['The replication method needs a finite timespan, e.g. ' ...
            'SolverJMT(model,''timespan'',[0,10]).'];
        return
    end
end

% JMVA implements RECAL, CoMoM, Chow, Bard-Schweitzer, AQL, Linearizer and De
% Souza-Muntz Linearizer for SINGLE-SERVER stations only, which is why
% writeJMVA refuses the model rather than writing an <ldstation> the algorithm
% cannot read. Multiserverhood is a station count and not a declared feature,
% so it cannot ride in getMethodFeatureSet the way the load-dependent scaling
% of the same restriction does.
if jmtJmvaIsClosedOnly(method)
    ns = sn.nservers(isfinite(sn.nservers));
    if ~isempty(ns) && max(ns) > 1
        reason = sprintf('%s does not support multi-server stations.', method);
        return
    end
end

% JMVA'S MULTICLASS CoMoM ANSWERS A DIFFERENT MODEL EVERY RUN, under a name
% that promises an exact normalizing constant. Measured against the JMT.jar in
% common/ (1.2.x, built 2026-04-30): SolverMultiClosedCoMoM.solve wraps its
% whole primary solve in catch(Exception), and when CoMoMBTFSolver reports
% "LUP Decomposition failed: Singular Matrix" it PERTURBS every demand (held as
% an integer scaled by 100000) by (i+1) + round(100*U), U from a MersenneTwister
% that RandomEngine.makeDefault() seeds off Math.random(), re-solves, and emits
% the answer as <algorithm name="CoMoM">, ok="true" and successful="true" on
% every measure. Nothing distinguishes that document from an exact one.
% The perturbation is about 1e-3 of a demand, so the answer looks plausible and
% is neither exact nor repeatable: the -seed the JMT CLI takes is the SIMULATION
% seed and never reaches the analytical engine. On the .jmva of
% sanity_CQN_2q_psfcfs_1class_1slcateachqueue, five runs at seed 23000 spread
% 0.22% against a 0.1% CoarseTol, and each reported a self-looping class holding
% jobs at the station it never visits, which no tolerance can excuse. RECAL on
% the same document returns the exact CTMC values, bit-identical on every run.
% WHY THE RULE IS THE CHAIN COUNT and not the singularity itself: the singular
% case is a rank condition inside JMT's BTF basis, and a demand-matrix
% characterisation of it checked against 400 random models was wrong on 21% of
% them, so LINE cannot tell the exact path from the perturbed one either before
% or after the run. A SINGLE-CHAIN model is exempt and stays runnable:
% SolverDispatcher.solveSingle serves SolverAlgorithm.COMOM with
% SolverSingleClosedRECAL, which is exact and reproducible; only the multiclass
% entry point is the broken one.
if any(strcmpi(method, {'jmva.comom','jmt.jmva.comom'})) && sn.nchains > 1
    reason = sprintf(['%s is refused on a model with %d chains: JMT''s multiclass CoMoM ' ...
        'engine silently re-solves a randomly perturbed model whenever its exact linear ' ...
        'system is singular, and reports that as an exact CoMoM solution, so the answer is ' ...
        'wrong and different on every run. Use ''jmva.recal'', which is exact and ' ...
        'reproducible on the same document, ''jmva''/''jmva.mva'' for exact MVA, or ' ...
        'SolverNC for a LINE-native normalizing constant.'], method, sn.nchains);
    return
end

% A BINDING FINITE BUFFER, which neither engine can carry, though for opposite
% reasons -- so the test is shared and the verdict is not.
%
% What makes a buffer BIND is not that sn.cap is finite: refreshCapacity DERIVES
% a finite cap for every station nobody capped. It is that the cap is strictly
% below the population that can REACH the station, which is the writer's own
% test (jmtReachablePopulation), and an infinite-server station has no buffer at
% all. Both are taken from @JMTIO/saveBufferCapacity, so the gate binds exactly
% where the writer binds.
%
% The ENGINE argument decides the verdict; see the header. Without one the
% method name does, and a name that is neither engine's gets no verdict at all.
if nargin >= 4 && ~isempty(engine)
    isJmva = strcmpi(engine, 'jmva');
    isJsim = strcmpi(engine, 'jsim');
else
    isJmva = strncmpi(method, 'jmva', 4);
    isJsim = any(strcmpi(method, {'default','jsim','replication'}));
end
if ~isJmva && ~isJsim
    return
end

% IMMEDIATE FEEDBACK: a job that self-loops keeps its server instead of
% re-queueing, which neither JMT document can state. runAnalyzer used to WARN
% and return no solution, so the gate called the pair runnable and the table
% came back empty. Keyed on the CALLER: without an ENGINE argument the asker is
% SolverJMT's own gate or runAnalyzer. writeJMVA passes 'jmva' on behalf of
% SolverQNS too, which refuses the same feature in its own words
% (qns_immfeed_refusal), so this sentence must not be handed to it.
askedByJmt = nargin < 4 || isempty(engine);
if askedByJmt && isfield(sn,'immfeed') && ~isempty(sn.immfeed) && any(sn.immfeed(:))
    reason = ['SolverJMT does not support immediate feedback (sn.immfeed): neither the JSIM ' ...
        'nor the JMVA document can keep a self-looping job on its server. Use SolverCTMC ' ...
        'or SolverSSA, whose state space carries the self-loop.'];
    return
end

% AN EDD OR EDF STATION WHOSE CLASSES CARRY NO DUE DATE. JSIM orders such a
% buffer by a per-station deadline it reads from <classSoftDeadlines>, and a
% job that arrives without one aborts the run instead of being served FCFS.
% JSIM only: getJMVAFeatureSet withdraws SchedStrategy_EDD and _EDF outright,
% since that document carries no discipline beyond the BCMP ones.
if isJsim
    reason = jmtDeadlineRefusal(sn);
    if ~isempty(reason)
        return
    end
end

if isJmva
    % THE JMVA DOCUMENT CARRIES A MEAN DEMAND AND A VISIT COUNT PER CHAIN and
    % nothing else, so two rules no feature name can state are asked here, of
    % the engine and not of the name: SolverQNS writes this same document on
    % its qnsolver path while declaring Fork/Join for its lqns path.
    if sn_has_fork_join(sn)
        reason = ['The JMVA document has no fork or join element, so the analytical engine ' ...
            'would solve the model without the synchronisation. Use a simulator (SolverJMT ' ...
            '''jsim'', SolverLDES) or SolverMVA, which approximates fork-join.'];
        return
    end
    % A non-exponential law at a station whose discipline is NOT insensitive
    % (FCFS, SIRO, non-preemptive LCFS): the mean alone decides nothing there,
    % and the engine would answer for the exponential law under the user's
    % name. SolverMVA's AMVA corrections read the SCV, so it is not refused.
    for ist = 1:sn.nstations
        sc = sn.sched(ist);
        if sc ~= SchedStrategy.FCFS && sc ~= SchedStrategy.SIRO && sc ~= SchedStrategy.LCFS
            continue
        end
        for r = 1:sn.nclasses
            if isfinite(sn.rates(ist,r)) && sn.rates(ist,r) > 0 && isfinite(sn.scv(ist,r)) ...
                    && abs(sn.scv(ist,r) - 1) > GlobalConstants.FineTol
                reason = sprintf(['Station %s serves class %s under %s with a non-exponential ' ...
                    'law (SCV %.4g). The JMVA document carries the mean demand only, so the ' ...
                    'analytical engine would answer for the exponential law. Use a simulator ' ...
                    '(SolverJMT ''jsim'', SolverLDES), SolverMVA, whose AMVA corrections read ' ...
                    'the SCV, or SolverCTMC or SolverSSA.'], ...
                    sn.nodenames{sn.stationToNode(ist)}, sn.classnames{r}, upper(SchedStrategy.toText(sc)));
                return
            end
        end
    end
end

for ist = 1:sn.nstations
    % A SOURCE AND A SINK HAVE NO BUFFER THAT CAN BIND. The Source IS the
    % external world and the Sink absorbs, so neither ever holds a job that a
    % capacity could refuse, yet refreshCapacity writes them a row like any
    % other station. Excluded on NODE TYPE, as NetworkSolver.checkBindingCapacity
    % and the C++ qn::binding_capacity_reason exclude them, and not by name.
    ntype = sn.nodetype(sn.stationToNode(ist));
    if ntype == NodeType.Source || ntype == NodeType.Sink
        continue
    end
    % UNBOUNDED IS Inf HERE, and that is worth saying: the JAR cannot carry Inf
    % on sn.cap because Station.cap is an int whose "no bound" value is
    % Integer.MAX_VALUE, and refreshCapacity SUMS that sentinel across the
    % classes served at a station -- so a mixed station comes out as
    % 2147483647 + N there and needs SaveHandlers.jmtCapIsUnbounded. MATLAB,
    % native python and C++ all default station.cap to Inf, so isfinite is the
    % whole test.
    if ~isfinite(sn.cap(ist)) || sn.cap(ist) >= jmtReachablePopulation(sn, ist)
        continue
    end
    if isinf(sn.nservers(ist))
        continue
    end
    if isJmva
        % JMVA is refused OUTRIGHT: writeJMVA emits a station type, a per-chain
        % service demand and a per-chain visit count and nothing else, so the
        % document has no capacity element for the buffer to ride in and the
        % model would be solved as if it were unbounded. Measured on a closed
        % Delay+FCFS model, N=4, cap 2: every jmva method reported 2.19 jobs at
        % a station that can hold 2, against the exact 1.33. This is the rule
        % NetworkSolver.checkBindingCapacity already applies to SolverMVA,
        % SolverNC and SolverQNS -- and SolverQNS writes THIS SAME DOCUMENT, so
        % the jmva arm was the one hole in it.
        reason = sprintf(['Station %s carries a finite capacity %d that binds. The JMVA ' ...
            'document has no capacity element at all, so the analytical engine would solve ' ...
            'the model as if the buffer were unbounded and report that as the answer. Use ' ...
            'the ''jsim'' method, which exports the buffer with its drop rule when JMT can ' ...
            'express it, or SolverCTMC, SolverSSA or SolverLDES.'], ...
            sn.nodenames{sn.stationToNode(ist)}, sn.cap(ist));
        return
    end
    % JSIM exports the buffer, but only for the rules JMT can read: the shared
    % predicate saveBufferCapacity raises decides which, so an open loss buffer
    % and a declared BAS one stay runnable and only the cases JMT would answer
    % unconstrained are refused.
    reason = jmtStationCapRefusal(sn, ist);
    if ~isempty(reason)
        return
    end
end
end
