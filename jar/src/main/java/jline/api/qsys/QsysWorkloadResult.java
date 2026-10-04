/**
 * @file Result of a workload distribution computation for a single station.
 *
 * @since LINE 3.0
 */
package jline.api.qsys;

/**
 * Holder for the stationary distribution of the quantity of work in a single
 * service stage, together with the grid it is reported on and the stationary
 * distribution of the number of requests in the system.
 */
public final class QsysWorkloadResult {

    /** F[j] = Pr{psi &lt;= t[j]}, the workload CDF. */
    public final double[] F;

    /** Grid at which F is reported. */
    public final double[] t;

    /** p[h] = Pr{x = h}, the stationary number in system, h = 0..N. */
    public final double[] p;

    public QsysWorkloadResult(double[] F, double[] t, double[] p) {
        this.F = F;
        this.t = t;
        this.p = p;
    }
}
