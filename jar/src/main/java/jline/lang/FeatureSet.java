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
        set.put("InfiniteServer", false);
        set.put("Forker", false);
        set.put("Joiner", false);
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
        set.put("RoutingStrategy_RL", false);
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
        set.put("SetupDelayOff", false);
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
        List<String> unsupported = new LinkedList<>();
        used.set.forEach((usedFeat, val) -> {
            if (val && !supported.set.getOrDefault(usedFeat, false)) {
                unsupported.add(usedFeat);
            }
        });
        if (!unsupported.isEmpty()) {
            return "Some features are not supported by the chosen solver (feature: "
                    + String.join(", ", unsupported) + ").";
        }
        return "";
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
