/**
 * @file Result class for CTMC reward computation
 *
 * @since LINE 3.0
 */
package jline.solvers.ctmc.analyzers;

import java.util.List;
import java.util.Map;
import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result class for CTMC reward computation via value iteration.
 */
public final class RewardResult {
    private final Map<String, Matrix> valueFunction;
    private final double[] time;
    private final List<String> rewardNames;
    private final Matrix stateSpace;
    private final Map<String, Double> steadyState;
    private final double runtime;

    public RewardResult(Map<String, Matrix> valueFunction, double[] time, List<String> rewardNames,
                        Matrix stateSpace, Map<String, Double> steadyState, double runtime) {
        this.valueFunction = valueFunction;
        this.time = time;
        this.rewardNames = rewardNames;
        this.stateSpace = stateSpace;
        this.steadyState = steadyState;
        this.runtime = runtime;
    }

    public Map<String, Matrix> getValueFunction() { return valueFunction; }
    public double[] getTime() { return time; }
    public List<String> getRewardNames() { return rewardNames; }
    public Matrix getStateSpace() { return stateSpace; }
    public Map<String, Double> getSteadyState() { return steadyState; }
    public double getRuntime() { return runtime; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof RewardResult)) return false;
        RewardResult that = (RewardResult) o;
        return Double.compare(that.runtime, runtime) == 0
                && Objects.equals(valueFunction, that.valueFunction)
                && java.util.Arrays.equals(time, that.time)
                && Objects.equals(rewardNames, that.rewardNames)
                && Objects.equals(stateSpace, that.stateSpace)
                && Objects.equals(steadyState, that.steadyState);
    }

    @Override
    public int hashCode() {
        int result = Objects.hash(valueFunction, rewardNames, stateSpace, steadyState, runtime);
        result = 31 * result + java.util.Arrays.hashCode(time);
        return result;
    }

    @Override
    public String toString() {
        return "RewardResult(rewardNames=" + rewardNames + ", runtime=" + runtime + ")";
    }
}
