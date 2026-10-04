package jline.api.mapqn;

import java.util.Map;

/**
 * Solution returned by Mapqn LP solvers.
 */
public class Mapqn_solution {
    private final double objectiveValue;
    private final Map<String, Double> variables;

    public Mapqn_solution(double objectiveValue, Map<String, Double> variables) {
        this.objectiveValue = objectiveValue;
        this.variables = variables;
    }

    public double getObjectiveValue() { return objectiveValue; }
    public Map<String, Double> getVariables() { return variables; }

    public double getVariable(String name) {
        Double v = variables.get(name);
        return v == null ? 0.0 : v;
    }

    /**
     * Get utilization for queue i, phase k (1-based indices).
     * Reads the {@code U_i_k} variable populated by the linear-reduction
     * (LR) bound model, with an {@code e_i_k} fallback for solvers that
     * use that naming convention.
     */
    /**
     * Get utilization for queue i (1-based), aggregated over phases, as written
     * by the QR bound models. The phase-indexed overload below is the LR form.
     */
    public double getUtilization(int i) {
        return getVariable("U_" + i);
    }

    public double getUtilization(int i, int k) {
        double uVar = getVariable("U_" + i + "_" + k);
        if (uVar != 0.0) return uVar;
        return getVariable("e_" + i + "_" + k);
    }

    /**
     * Get queue length for queue i, phase k (1-based indices).
     * Reads the {@code Q_i_k} variable populated by the LR bound model.
     */
    public double getQueueLength(int i, int k) { return getVariable("Q_" + i + "_" + k); }

    /**
     * Get utilization from the MVA-version model, which uses the
     * {@code UN_i_k} naming convention (see {@link Mapqn_bnd_lr_mva}).
     */
    public double getUtilizationMVA(int i, int k) { return getVariable("UN_" + i + "_" + k); }
    public double getQueueLengthMVA(int i, int k) { return getVariable("QN_" + i + "_" + k); }
}
