/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import java.io.Serializable;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;

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
        set.put("Gamma", false);
        set.put("HyperExp", false);
        set.put("Immediate", false);
        set.put("Lognormal", false);
        set.put("MAP", false);
        set.put("MMAP", false);
        set.put("BMAP", false);
        set.put("MMPP2", false);
        set.put("Normal", false);
        set.put("Pareto", false);
        set.put("PH", false);
        set.put("ME", false);
        set.put("RAP", false);
        set.put("Replayer", false);
        set.put("Uniform", false);
        set.put("Weibull", false);
        set.put("StatelessClassSwitcher", false);
        set.put("CacheClassSwitcher", false);
        set.put("InfiniteServer", false);
        set.put("Forker", false);
        set.put("Joiner", false);
        set.put("LogTunnel", false);
        set.put("SharedServer", false);
        set.put("Buffer", false);
        set.put("Region", false);
        set.put("Linkage", false);
        set.put("Enabling", false);
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
        set.put("RoutingStrategy_KCHOICES", false);
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
        set.put("ClosedClass", false);
        set.put("OpenClass", false);
        set.put("SelfLoopingClass", false);
        set.put("LoadDependence", false);
    }

    /**
     * Checks if the used features are supported by the given solver
     *
     * @param supported - the features supported by the solver
     * @param used      - the used features
     * @return - true if the used features are supported, false otherwise
     */
    public static boolean supports(FeatureSet supported, FeatureSet used) {
        List<String> unsupported = new LinkedList<>();
        used.set.forEach((usedFeat, val) -> {
            if (val && !supported.set.get(usedFeat)) {
                unsupported.add(usedFeat);
            }
        });
        if (!unsupported.isEmpty()) {
            System.out.println("Some features are not supported by the chosen solver:");
            unsupported.forEach(System.out::println);
            return false;
        }
        return true;
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
