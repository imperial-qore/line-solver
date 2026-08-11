package jline.opt.objectives;

import jline.lang.nodes.Station;
import jline.opt.results.EvaluationResult;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Minimize infrastructure cost subject to service-level constraints. Cost is a
 * weighted sum of server costs ({@code name_servers}), rate costs (variable
 * names containing the station name and "rate"), and replica costs
 * ({@code name_replicas}). Mirrors native-Python {@code MinimizeCost}.
 */
public class MinimizeCost extends Objective {

    private final Map<String, Double> serverCost = new LinkedHashMap<String, Double>();
    private final Map<String, Double> rateCost = new LinkedHashMap<String, Double>();
    private final Map<String, Double> replicaCost = new LinkedHashMap<String, Double>();

    public MinimizeCost(Map<Station, Double> serverCost, Map<Station, Double> rateCost,
                        Map<Station, Double> replicaCost, List<Constraint> subjectTo) {
        if (serverCost != null) {
            for (Map.Entry<Station, Double> e : serverCost.entrySet()) {
                this.serverCost.put(e.getKey().getName() + "_servers", e.getValue());
            }
        }
        if (rateCost != null) {
            for (Map.Entry<Station, Double> e : rateCost.entrySet()) {
                this.rateCost.put(e.getKey().getName(), e.getValue());
            }
        }
        if (replicaCost != null) {
            for (Map.Entry<Station, Double> e : replicaCost.entrySet()) {
                this.replicaCost.put(e.getKey().getName() + "_replicas", e.getValue());
            }
        }
        if (subjectTo != null) {
            this.constraints = new ArrayList<Constraint>(subjectTo);
        }
    }

    /** Server-cost-only convenience constructor. */
    public static MinimizeCost serverCost(Map<Station, Double> serverCost, List<Constraint> subjectTo) {
        return new MinimizeCost(serverCost, null, null, subjectTo);
    }

    /** Rate-cost-only convenience constructor. */
    public static MinimizeCost rateCost(Map<Station, Double> rateCost, List<Constraint> subjectTo) {
        return new MinimizeCost(null, rateCost, null, subjectTo);
    }

    public boolean isMinimization() {
        return true;
    }

    public double evaluate(EvaluationResult result, Map<String, Object> variableValues) {
        double total = 0.0;

        for (Map.Entry<String, Double> e : serverCost.entrySet()) {
            Object value = variableValues.get(e.getKey());
            if (ObjectiveUtil.isScalar(value)) {
                total += e.getValue() * ((Number) value).doubleValue();
            }
        }
        for (Map.Entry<String, Double> e : rateCost.entrySet()) {
            String pattern = e.getKey();
            for (Map.Entry<String, Object> vv : variableValues.entrySet()) {
                String vname = vv.getKey();
                if (vname.contains(pattern) && vname.toLowerCase().contains("rate")
                        && ObjectiveUtil.isScalar(vv.getValue())) {
                    total += e.getValue() * ((Number) vv.getValue()).doubleValue();
                }
            }
        }
        for (Map.Entry<String, Double> e : replicaCost.entrySet()) {
            Object value = variableValues.get(e.getKey());
            if (ObjectiveUtil.isScalar(value)) {
                total += e.getValue() * ((Number) value).doubleValue();
            }
        }
        return total;
    }
}
