package odesolver;

/**
 * A per-accepted-step stopping test for {@link LSODA}.
 *
 * <p>WHY THIS EXISTS RATHER THAN A StepHandler. {@code LSODA} implements
 * {@code FirstOrderIntegrator}, so it carries {@code addStepHandler} and
 * {@code addEventHandler} -- and both are EMPTY STUBS that accept a handler and
 * discard it. Code registering one compiles, runs, and is silently never
 * called. This interface is the honest hook: {@code LSODA} consults it at the
 * point it records an accepted step, which is the same instant MATLAB's
 * OutputFcn and the native-Python step loop test at, so the four codebases stop
 * on the same condition at the same place.</p>
 *
 * <p>THE STATE IS 1-INDEXED. {@code LSODA} carries its solution in
 * {@code y[1..n]} with {@code y[0]} unused, the Fortran convention the
 * translation kept; an implementation must read from index 1, not 0.</p>
 *
 * @see LSODA#setStopCondition(StepStopCondition)
 */
public interface StepStopCondition {

    /**
     * @param t the time of the accepted step
     * @param y the state at {@code t}, 1-indexed: the coordinates are
     *          {@code y[1] .. y[n]} and {@code y[0]} is unused
     * @return true to end the integration here, reporting {@code t} as the
     *         window's end
     */
    boolean stop(double t, double[] y);
}
