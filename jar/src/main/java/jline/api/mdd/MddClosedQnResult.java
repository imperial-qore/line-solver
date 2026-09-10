package jline.api.mdd;

import jline.util.matrix.Matrix;

/** Result of the MDD-stored exact closed-network solve of {@link Mdd_closedqn}. */
public class MddClosedQnResult {

    /** The MDD holding the reachable occupancy set. */
    public final MDD mdd;
    /** CTMC generator, rows aligned to {@link MDD#index} order. */
    public final Matrix Q;
    /** Stationary distribution over the reachable states, MDD index order. */
    public final Matrix pi;
    /** |S| x M occupancy states, in MDD index order. */
    public final int[][] states;
    /** Mean number of jobs per station. */
    public final double[] QLen;
    /** Utilisation: busy servers / servers, or mean busy jobs for a delay. */
    public final double[] U;
    /** Per-station throughput. */
    public final double[] X;
    /** MDD storage statistics of the reachable set. */
    public final MddStats stats;
    /** Reachable-set build time in seconds; 0 when the diagram was supplied. */
    public final double timeReach;
    /** Generator assembly time in seconds. */
    public final double timeGen;
    /** ctmc_solve time in seconds. */
    public final double timeSolve;
    /** Performance-measure time in seconds. */
    public final double timeMetrics;

    public MddClosedQnResult(MDD mdd, Matrix Q, Matrix pi, int[][] states, double[] QLen,
                             double[] U, double[] X, MddStats stats, double timeReach,
                             double timeGen, double timeSolve, double timeMetrics) {
        this.mdd = mdd;
        this.Q = Q;
        this.pi = pi;
        this.states = states;
        this.QLen = QLen;
        this.U = U;
        this.X = X;
        this.stats = stats;
        this.timeReach = timeReach;
        this.timeGen = timeGen;
        this.timeSolve = timeSolve;
        this.timeMetrics = timeMetrics;
    }
}
