/**
 * @file Robust Queueing Theory (RQT) worst-case system time of a G/G/1 queue
 *
 * Single-server case of {@link jline.api.qsys.Qsys_gigk_rqt}, per C. Bandi,
 * D. Bertsimas, N. Youssef (2015), "Robust Queueing Theory", Operations
 * Research 63(3), 676-700, Theorem 2 and eq. (12).
 *
 * @since LINE 3.0
 */
package jline.api.qsys;

public final class Qsys_gig1_rqt {
    private Qsys_gig1_rqt() {}

    /**
     * Robust Queueing Theory worst-case system time of a G/G/1 FCFS queue.
     *
     * @param lambda  arrival rate
     * @param mu      service rate
     * @param Gamma_a variability parameter of the arrival uncertainty set
     * @param Gamma_s variability parameter of the service uncertainty set
     * @param alpha_a arrival tail coefficient in (1,2]
     * @param alpha_s service tail coefficient in (1,2]
     * @return array {W, rhohat, Sworst}: closed-form bound on the system time
     *         (Theorem 2), modified utilization, and exact worst-case system
     *         time over the uncertainty sets (eq. 12)
     */
    public static double[] qsys_gig1_rqt(double lambda, double mu, double Gamma_a, double Gamma_s,
                                         double alpha_a, double alpha_s) {
        return Qsys_gigk_rqt.qsys_gigk_rqt(lambda, mu, Gamma_a, Gamma_s, 1, alpha_a, alpha_s);
    }

    /**
     * Worst-case system time with the finite-variance tail coefficients
     * alpha_a = alpha_s = 2.
     *
     * @param lambda  arrival rate
     * @param mu      service rate
     * @param Gamma_a variability parameter of the arrival uncertainty set
     * @param Gamma_s variability parameter of the service uncertainty set
     * @return array {W, rhohat, Sworst}
     */
    public static double[] qsys_gig1_rqt(double lambda, double mu, double Gamma_a, double Gamma_s) {
        return Qsys_gigk_rqt.qsys_gigk_rqt(lambda, mu, Gamma_a, Gamma_s, 1, 2.0, 2.0);
    }
}
