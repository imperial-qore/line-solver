/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import java.io.Serializable;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.line_warning;
import static jline.io.InputOutput.mfilename;

/**
 * Class representing the features of a particular solver
 */
public class FeatureSet implements Serializable {
    private final HashMap<String, Boolean> set;

    /**
     * Creates a new FeatureSet with all features initialized to false.
     * This represents the default state where no features are enabled.
     */
    public FeatureSet() {
        this.set = new HashMap<>();
        set.put("ClassSwitch", false);
        set.put("Cache", false);
        set.put("Delay", false);
        set.put("DelayStation", false);
        set.put("Fork", false);
        set.put("Join", false);
        set.put("Logger", false);
        set.put("Place", false);
        set.put("QueueingPlace", false);
        set.put("Queue", false);
        set.put("JobSink", false);
        set.put("Source", false);
        set.put("Sink", false);
        set.put("Router", false);
        set.put("Transition", false);
        set.put("Coxian", false);
        set.put("Cox2", false);
        set.put("APH", false);
        set.put("Det", false);
        set.put("Disabled", false);
        set.put("Erlang", false);
        set.put("Exp", false);
        set.put("Bernoulli", false);
        set.put("Binomial", false);
        set.put("Gamma", false);
        set.put("Geometric", false);
        set.put("Poisson", false);
        set.put("HyperExp", false);
        set.put("Immediate", false);
        set.put("Lognormal", false);
        set.put("MAP", false);
        set.put("DMAP", false);
        set.put("MMAP", false);
        set.put("BMAP", false);
        set.put("MMPP2", false);
        set.put("NHPP", false);
        set.put("MAPt", false);
        set.put("PHt", false);
        // Name declared by EmpiricalCDF.java, which is what getUsedLangFeatures
        // marks; the JSON wire type "EmpiricalCDF" is unrelated to it.
        set.put("EmpiricalCdf", false);
        set.put("Expolynomial", false);
        set.put("Normal", false);
        set.put("Pareto", false);
        set.put("PH", false);
        set.put("ME", false);
        set.put("RAP", false);
        set.put("Replayer", false);
        set.put("Trace", false);
        set.put("Uniform", false);
        set.put("Weibull", false);
        // Distributions that no solver supports as a service or arrival
        // process, but that Distribution.getName() can still emit into
        // getUsedLangFeatures. They must be registered: an unregistered name
        // makes setTrue line_error, i.e. a crash instead of the intended
        // "feature not supported by the chosen solver" rejection. No solver
        // declares them, so a model using one is rejected cleanly.
        set.put("DiscreteSampler", false);
        set.put("DiscreteUniform", false);
        set.put("Empirical", false);
        set.put("GMM", false);
        set.put("MMDP", false);
        set.put("MMDP2", false);
        set.put("MMPP", false);
        set.put("MultivariateNormal", false);
        set.put("NegBinomial", false);
        set.put("Prior", false);
        set.put("Zipf", false);
        set.put("StatelessClassSwitcher", false);
        set.put("CacheClassSwitcher", false);
        set.put("CacheRetrieval", false);
        set.put("CacheItemSize", false);
        set.put("InfiniteServer", false);
        set.put("Forker", false);
        set.put("Joiner", false);
        // quorum join: k of n siblings, k < n
        set.put("JoinPartial", false);
        set.put("LogTunnel", false);
        set.put("SharedServer", false);
        set.put("Buffer", false);
        // "Region" (finite capacity region) is DECLARATIVE ONLY in the JAR and
        // in MATLAB: SolverLDES/SolverJMT advertise it via setTrue, but
        // Network.getUsedLangFeatures never records it, so it never gates a
        // model. It is kept here because the registry is the synced canonical
        // name list (MATLAB SolverFeatureSet.m, python lang/__init__.py) and
        // because setTrue("Region") line_errors on an unrecognized name, so
        // removing it would break the LDES/JMT featsets.
        // Do NOT start recording it without first auditing the FCR-capable
        // featsets: MATLAB hard-gates through runAnalyzerChecks, and JMT/CTMC
        // /SSA all accept FCR models today while declaring no "Region", as does
        // the NC single-Delay loss-network path. Recording it would reject all
        // of those. FCR is instead rejected imperatively where it is wrong, in
        // SolverMVA.runAnalyzer, SolverNC.runAnalyzer and -- since 2026-07-17 --
        // SolverMAM.runAnalyzer and SolverFluid.runAnalyzer (MAM/FLD used to
        // accept FCR and silently return the UNCONSTRAINED answer: MAM 4.000,
        // i.e. exactly rho/(1-rho), and FLD 0.800, vs CTMC's exact 0.8525).
        // NOTE python IS the exception: it records "Region" in
        // get_used_lang_features, so its featset gate does reject FCR, and its
        // CTMC declares "Region" deliberately to stay accepting. python's MVA/NC
        // featset checks are no longer advisory at the solve path either (both
        // enforce via supportsModelMethod -> runAnalyzerChecks since 2026-07-17).
        // See BUGS.md BUG-39 / NEW-1 / NEW-2 and _kb/06-solver-catalog.md.
        set.put("Region", false);
        set.put("Linkage", false);
        set.put("Enabling", false);
        set.put("Inhibiting", false);
        set.put("Timing", false);
        set.put("Firing", false);
        set.put("Storage", false);
        set.put("RandomSource", false);
        set.put("Dispatcher", false);
        set.put("Server", false);
        set.put("ServiceTunnel", false);
        set.put("RoutingStrategy_PROB", false);
        set.put("RoutingStrategy_RAND", false);
        set.put("RoutingStrategy_RROBIN", false);
        set.put("RoutingStrategy_WRROBIN", false);
        set.put("RoutingStrategy_JSQ", false);
        set.put("RoutingStrategy_SQ", false);
        set.put("RoutingStrategy_SDR", false);
        set.put("SchedStrategy_INF", false);
        set.put("SchedStrategy_FCFS", false);
        set.put("SchedStrategy_FCFSPR", false);
        set.put("SchedStrategy_FCFSPI", false);
        set.put("SchedStrategy_FCFSPRIO", false);
        set.put("SchedStrategy_FCFSPRPRIO", false);
        set.put("SchedStrategy_FCFSPIPRIO", false);
        set.put("SchedStrategy_LCFS", false);
        set.put("SchedStrategy_LCFSPR", false);
        set.put("SchedStrategy_LCFSPI", false);
        set.put("SchedStrategy_LCFSPRIO", false);
        set.put("SchedStrategy_LCFSPRPRIO", false);
        set.put("SchedStrategy_LCFSPIPRIO", false);
        set.put("SchedStrategy_SEPT", false);
        set.put("SchedStrategy_LEPT", false);
        set.put("SchedStrategy_SJF", false);
        set.put("SchedStrategy_LJF", false);
        set.put("SchedStrategy_SRPT", false);
        set.put("SchedStrategy_SRPTPRIO", false);
        set.put("SchedStrategy_PSJF", false);
        set.put("SchedStrategy_FB", false);
        set.put("SchedStrategy_LRPT", false);
        set.put("SchedStrategy_SETF", false);
        set.put("SchedStrategy_FSP", false);
        set.put("SchedStrategy_PAS", false);
        set.put("SchedStrategy_OI", false);
        set.put("SchedStrategy_PS", false);
        set.put("SchedStrategy_DPS", false);
        set.put("SchedStrategy_GPS", false);
        set.put("SchedStrategy_PSPRIO", false);
        set.put("SchedStrategy_DPSPRIO", false);
        set.put("SchedStrategy_GPSPRIO", false);
        set.put("SchedStrategy_SIRO", false);
        set.put("SchedStrategy_HOL", false);
        set.put("SchedStrategy_EXT", false);
        set.put("SchedStrategy_POLLING", false);
        set.put("SchedStrategy_EDD", false);
        set.put("SchedStrategy_EDF", false);
        set.put("SchedStrategy_LPS", false);
        set.put("ReplacementStrategy_RR", false);
        set.put("ReplacementStrategy_FIFO", false);
        set.put("ReplacementStrategy_SFIFO", false);
        set.put("ReplacementStrategy_LRU", false);
        set.put("ReplacementStrategy_HLRU", false);
        set.put("ReplacementStrategy_CLIMB", false);
        set.put("ReplacementStrategy_QLRU", false);
        set.put("ClosedClass", false);
        set.put("OpenClass", false);
        set.put("SelfLoopingClass", false);
        set.put("OpenSignal", false);
        set.put("ClosedSignal", false);
        set.put("SignalType_NEGATIVE", false);
        set.put("SignalType_REPLY", false);
        set.put("SignalType_CATASTROPHE", false);
        set.put("SignalBatchRemoval", false);
        set.put("BatchArrival", false);
        set.put("SignalRemovalPolicy", false);
        set.put("LoadDependence", false);
        set.put("ClassDependence", false);
        set.put("JointDependence", false);
        set.put("GlobalDependence", false);
        set.put("SetupDelayOff", false);
        set.put("ServerParallelism", false);
        set.put("Retrial", false);
        set.put("Balking", false);
        set.put("Reneging", false);
        set.put("Breakdown", false);
        // LayeredNetwork (LQN) features
        set.put("Host", false);
        set.put("Processor", false);
        set.put("Task", false);
        set.put("Entry", false);
        set.put("Activity", false);
        set.put("SyncCall", false);
        set.put("AsyncCall", false);
        set.put("ActivityPrecedence_PRE_SEQ", false);
        set.put("ActivityPrecedence_POST_SEQ", false);
        set.put("ActivityPrecedence_PRE_AND", false);
        set.put("ActivityPrecedence_POST_AND", false);
        set.put("ActivityPrecedence_PRE_OR", false);
        set.put("ActivityPrecedence_POST_OR", false);
        set.put("SchedStrategy_REF", false);
        // Layered cache-queueing models: a CacheTask holds a segmented cache
        // whose items are ItemEntry entries, and a read is a call to one of them
        // whose bound activity branches on a POST_CACHE precedence into hit and
        // miss. SolverLDES.getLNFeatureSet DECLARES all three, and setTrue
        // line_errors on an unregistered name, so that whole feature set threw
        // before it could compare anything: supports(LayeredNetwork) could not
        // run at all, which is also why nobody noticed SchedStrategy_REF was
        // never recorded.
        set.put("CacheTask", false);
        set.put("ItemEntry", false);
        set.put("ActivityPrecedence_POST_CACHE", false);
        // Variable forking levels. A Fork emits tasksPerLink jobs on every
        // outgoing link; these three name the ways that degree stops being one
        // number. ForkFanoutVector: the count differs by destination or by class
        // (ForkNodeParam.fanOutLink). ForkFanoutRandom: the count is a draw from
        // a DiscreteSampler, redrawn per link and per forked job
        // (ForkNodeParam.fanOutDist). ForkBranchProbability: a branch fires only
        // with probability p, so the SIBLING COUNT is random even when each link
        // carries a fixed number (ForkNodeParam.fanOutProb). Appended at the
        // tail so every earlier index is unchanged.
        set.put("ForkFanoutVector", false);
        set.put("ForkFanoutRandom", false);
        set.put("ForkBranchProbability", false);
        // Two names the C++ port carried alone until 2026-08-22, for
        // capabilities this codebase can express but no gate could see.
        // HeteroServers: Queue.addServerType gives a station several server
        // POOLS with their own counts, class compatibilities and per-(type,
        // class) rates. Only SolverJMT and the LDES engine honour them; every
        // other solver reads sn.nservers and answers for a homogeneous station,
        // which is a different system. DepartureDiscipline:
        // Place.setDepartureDiscipline(class, FIFO) makes the depository
        // release a served token only after the earlier ones, which changes
        // which transitions are enabled. NO solver implements it in any
        // codebase, so declaring it nowhere is the point -- the model is
        // refused instead of being solved as if it were Normal.
        set.put("HeteroServers", false);
        set.put("DepartureDiscipline", false);
        // Ported from MATLAB SolverFeatureSet.m (2026-09-05), in this order, so
        // that every earlier index is unchanged. Both are "having something"
        // properties the registry could not name, so every rule about them
        // lived only in structural predicates and was invisible to the gate.
        // MultiServer: a finite-server station serving more than one job at
        // once, i.e. a Queue whose numberOfServers is finite and > 1 (a Delay /
        // INF station is not one). FiniteCapacity: a station or per-class
        // buffer that can BIND, exactly the condition Network.findBindingCapacity
        // tests (node-level cap / classCap below the population that can reach
        // it; an open class always binds; a Cache model is exempt), which is
        // also what NetworkSolver.bindingCapacityReason refuses on.
        set.put("MultiServer", false);
        set.put("FiniteCapacity", false);
    }

    /**
     * Checks if the used features are supported by the given solver
     *
     * @param supported - the features supported by the solver
     * @param used      - the used features
     * @return - true if the used features are supported, false otherwise
     */
    public static boolean supports(FeatureSet supported, FeatureSet used) {
        String reason = supportsReason(supported, used);
        if (!reason.isEmpty()) {
            line_warning(mfilename(new Object() {}), reason);
            return false;
        }
        return true;
    }

    /**
     * Same feature comparison as {@link #supports} but returns a human-readable
     * reason instead of a boolean (empty string when every used feature is
     * supported), and without the side-effect warning. Intended for the
     * method-aware gate and for feature-driven method selection, which probe
     * several candidate methods and must not warn on the non-covering ones.
     *
     * @param supported - the features supported by the (method of the) solver
     * @param used      - the used features
     * @return - empty string if supported, else the offending feature list
     */
    public static String supportsReason(FeatureSet supported, FeatureSet used) {
        List<String> unsupported = unsupportedFeatures(supported, used);
        if (!unsupported.isEmpty()) {
            return "Some features are not supported by the chosen solver (feature: "
                    + String.join(", ", unsupported) + ").";
        }
        return "";
    }

    /**
     * The names of the used features the given feature set does not cover, as a
     * list rather than a message. Dispatch decisions that depend on WHICH
     * features are missing (e.g. the MAP/MMPP random-environment fallback, which
     * fires only when the missing features are all non-renewal processes) read
     * this instead of parsing supportsReason.
     *
     * @param supported - the features supported by the (method of the) solver
     * @param used      - the used features
     * @return the offending feature names, empty when every used feature is supported
     */
    public static List<String> unsupportedFeatures(FeatureSet supported, FeatureSet used) {
        List<String> unsupported = new LinkedList<>();
        used.set.forEach((usedFeat, val) -> {
            if (!val || supported.set.getOrDefault(usedFeat, false)) {
                return;
            }
            // A specialization the solver did not name is covered by the general
            // capability when that one IS declared; see generalizationOf.
            String general = generalizationOf(usedFeat);
            if (general != null && supported.set.getOrDefault(general, false)) {
                return;
            }
            unsupported.add(usedFeat);
        });
        return unsupported;
    }

    /**
     * The registry name a specialization falls back to when it is not declared,
     * or null when the feature stands on its own.
     *
     * A few entries name a SPECIAL CASE of another entry rather than a
     * capability of their own: "Cox2" is a Coxian restricted to two phases and
     * "Trace" is a Replayer under another class name. Marking a model with only
     * the general name left the specific entry unreachable, which is dead
     * registry surface; marking it with the specific name alone would instead
     * REJECT the model at every solver that declares only the general one, i.e.
     * at every solver that accepts it today. So getFeatureName marks the most
     * specific name and the gate falls back here.
     *
     * The fallback runs one way only: a solver supporting just the special case
     * can still declare "Cox2" alone and keep refusing a five-phase Coxian.
     *
     * @param feature the used feature name
     * @return the more general feature name, or null if there is none
     */
    public static String generalizationOf(String feature) {
        if ("Cox2".equals(feature)) {
            return "Coxian";
        }
        if ("Trace".equals(feature)) {
            return "Replayer";
        }
        return null;
    }

    /**
     * Checks whether the given feature is used or not in the current feature set
     *
     * @param feature - the name of the given feature
     * @return - true if the given feature is used, false otherwise
     */
    public boolean inspectFeature(String feature) {
        return set.getOrDefault(feature, false);
    }

    /**
     * Returns the canonical registry of feature names (all features, enabled or
     * not). Used by cross-codebase parity tooling to enumerate and compare
     * per-method feature sets.
     *
     * @return the set of all feature names in this FeatureSet
     */
    public java.util.Set<String> featureNames() {
        return set.keySet();
    }

    /**
     * Sets multiple features to false in the feature set.
     * 
     * @param features array of feature names to disable
     * @throws RuntimeException if any feature name is not recognized
     */
    public void setFalse(String[] features) {
        for (String feature : features) {
            if (set.containsKey(feature))
                set.put(feature, false);
            else
                line_error(mfilename(new Object() {
                }), "Unrecognized feature to set to false in the feature set: " + feature);
        }
    }

    /**
     * Sets a single feature to false in the feature set.
     * 
     * @param feature the name of the feature to disable
     * @throws RuntimeException if the feature name is not recognized
     */
    public void setFalse(String feature) {
        if (set.containsKey(feature))
            set.put(feature, false);
        else
            line_error(mfilename(new Object() {
            }), "Unrecognized feature to set to true in the feature set: " + feature);
    }

    /**
     * Sets a single feature to true in the feature set.
     * 
     * @param feature the name of the feature to enable
     * @throws RuntimeException if the feature name is not recognized
     */
    public void setTrue(String feature) {
        if (set.containsKey(feature))
            set.put(feature, true);
        else
            line_error(mfilename(new Object() {
            }), "Unrecognized feature to set to true in the feature set: " + feature);
    }

    /**
     * Sets multiple features to true in the feature set.
     * 
     * @param features array of feature names to enable
     * @throws RuntimeException if any feature name is not recognized
     */
    public void setTrue(String[] features) {
        for (String feature : features) {
            if (set.containsKey(feature))
                set.put(feature, true);
            else
                line_error(mfilename(new Object() {
                }), "Unrecognized feature to set to true in the feature set: " + feature);
        }
    }
}
