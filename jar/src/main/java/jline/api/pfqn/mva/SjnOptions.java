/**
 * @file Options of the shortest-job-next response time equations
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

/** Options shared by Pfqn_mvasjn and Pfqn_amvasjn. */
public final class SjnOptions {
    /** Number of grid subdivisions of the job size axis, even. */
    public int ns = 32;
    /** Grid extent, in units of the largest mean service time at the station. */
    public double Lfactor = 8;
    /**
     * Priority levels, one per class, lower is higher priority. Null pools the classes and
     * compares every job by size; distinct levels select the priority reading, SJN within a
     * class under non-preemptive priority across classes.
     */
    public int[] prio = null;
    /** Convergence tolerance of the fixed point. */
    public double tol = 1e-8;
    /** Iteration cap of the fixed point. */
    public int iterMax = 1000;
    /** Utilization cap at an SJN station, strictly below one. */
    public double umax = 0.999;

    public SjnOptions() {}

    /** Reject the option combinations the equations cannot honour. */
    void validate(int R) {
        if (ns % 2 != 0) {
            throw new IllegalArgumentException("ns must be even, composite Simpson integrates over"
                    + " panels of two subdivisions");
        }
        if (umax <= 0 || umax >= 1) {
            throw new IllegalArgumentException("umax must lie strictly between zero and one");
        }
        if (prio != null) {
            if (prio.length != R) {
                throw new IllegalArgumentException("prio must have one priority level per class");
            }
            for (int i = 0; i < R; i++) {
                for (int j = i + 1; j < R; j++) {
                    if (prio[i] == prio[j]) {
                        throw new IllegalArgumentException("prio must assign distinct levels, ties"
                                + " across classes are not covered by the SJN priority equations");
                    }
                }
            }
        }
    }
}
