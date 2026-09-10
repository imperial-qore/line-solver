/**
 * @file Time-dependent index-of-dispersion function handle
 *
 * Functional interface for a vector-valued index of dispersion for counts (IDC)
 * evaluated at a scalar time argument, used by the Robust Queueing Network
 * Analyzer traffic equations.
 *
 * @since LINE 3.0
 */
package jline.api.npfqn;

public interface IdcFunction {
    /**
     * Evaluate the IDC vector at time t.
     *
     * @param t time argument (t &gt; 0)
     * @return array of per-queue IDC values
     */
    double[] eval(double t);
}
