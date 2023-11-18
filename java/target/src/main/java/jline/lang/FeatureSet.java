package jline.lang;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

/**
 * Class representing the features of a particular solver
 */
public class FeatureSet {
    private final HashMap<String, Boolean> set;

    public FeatureSet() {
        this.set = new HashMap<>();
        set.put("ClassSwitch", false);
        set.put("Cache", false);
        set.put("Delay", false);
        set.put("Fork", false);
        set.put("Join", false);
        set.put("Logger", false);
        set.put("Place", false);
        set.put("Queue", false);
        set.put("Sink", false);
        set.put("Source", false);
        set.put("Router", false);
        set.put("Transition", false);
        set.put("Coxian", false);
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
        set.put("MMPP2", false);
        set.put("Normal", false);
        set.put("Pareto", false);
        set.put("PH", false);
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
        set.put("SchedStrategy_LCFS", false);
        set.put("SchedStrategy_LCFSPR", false);
        set.put("SchedStrategy_SEPT", false);
        set.put("SchedStrategy_LEPT", false);
        set.put("SchedStrategy_DPS", false);
        set.put("SchedStrategy_GPS", false);
        set.put("SchedStrategy_LJF", false);
        set.put("SchedStrategy_SJF", false);
        set.put("SchedStrategy_PS", false);
        set.put("SchedStrategy_SIRO", false);
        set.put("SchedStrategy_HOL", false);
        set.put("SchedStrategy_EXT", false);
        set.put("ClosedClass", false);
        set.put("OpenClass", false);
    }

    public void setTrue(String[] features){
        for(String feature : features){
            if(set.containsKey(feature))
                set.put(feature, true);
            else
                System.err.println("Unrecognized feature to set to true in the feature set: " + feature);
        }
    }

    public void setFalse(String[] features){
        for(String feature : features){
            if(set.containsKey(feature))
                set.put(feature, false);
            else
                System.err.println("Unrecognized feature to set to false in the feature set: " + feature);
        }
    }

    /**
     * Checks if the used features are supported by the given solver
     * @param supported - the features supported by the solver
     * @param used - the used features
     * @return - true if the used features are supported, false otherwise
     */
    public static boolean supports(FeatureSet supported, FeatureSet used){
        List<String> unsupported = new LinkedList<>();
        used.set.forEach((usedFeat, val) -> {
            if(val && !supported.set.get(usedFeat)){
                unsupported.add(usedFeat);
            }
        });
        if(!unsupported.isEmpty()){
            System.out.println("Some features are not supported by the chosen solver:");
            unsupported.forEach(System.out::println);
            return false;
        }
        return true;
    }

    /**
     * Checks whether the given feature is used or not in the current feature set
     * @param feature - the name of the given feature
     * @return - true if the given feature is used, false otherwise
     */
    public boolean inspectFeature(String feature){
        return set.getOrDefault(feature, false);
    }
}
