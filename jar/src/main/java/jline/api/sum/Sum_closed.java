package jline.api.sum;

import jline.util.matrix.Matrix;

/**
 * Summation method (SUM) for closed queueing networks, including the
 * extended SUM (ESUM) node functions for non-product-form networks with
 * generally distributed service times.
 *
 * The method expresses the mean queue length of each station as a function
 * of its throughput, Ki = fi(lambdai), and solves the population constraint
 * sum_i Ki = K. Single-class models are solved by bisection on the system
 * throughput (Bolch et al., Sec. 9.2.1); multiclass models by Gauss-Seidel
 * sweeps of per-class bisections on the population constraints, a robust
 * alternative to the successive substitution of Sec. 9.2.2.
 *
 * Node functions:
 * - Product-form stations (scv=1, or insensitive disciplines PS/LCFS-PR,
 *   for which the caller must pass scv=1): Eq. (9.15)/(9.19).
 * - FCFS stations with general service (scv!=1): ESUM corrections,
 *   Eq. (10.88) for -/G/1 and Eq. (10.89) for -/G/m, with
 *   ai=(1+scv_i)/2 and Erlang-C waiting probability P_mi.
 * - Infinite-server stations (mi=Inf) and think times Z: Ki = lambdai*Li.
 *
 * Reference: G. Bolch, S. Greiner, H. de Meer, K.S. Trivedi, Queueing
 * Networks and Markov Chains, 2nd ed., Wiley, 2006, Secs. 9.2 and 10.1.4.4.
 */
public class Sum_closed {

    /** Result of the summation method. */
    public static final class Result {
        /** 1xR class throughputs */
        public final Matrix XN;
        /** MxR mean queue lengths */
        public final Matrix QN;
        /** MxR utilizations (per-server for queueing stations, X*L for IS) */
        public final Matrix UN;
        /** MxR residence times, RN=QN/XN */
        public final Matrix RN;
        /** number of iterations */
        public final int it;

        Result(Matrix XN, Matrix QN, Matrix UN, Matrix RN, int it) {
            this.XN = XN;
            this.QN = QN;
            this.UN = UN;
            this.RN = RN;
            this.it = it;
        }
    }

    /**
     * Summation method for closed queueing networks.
     *
     * @param L       MxR service demand matrix, L(i,r) = e(i,r)/mu(i,r)
     * @param N       1xR population vector
     * @param Z       1xR think times (aggregated as a delay term)
     * @param mi      Mx1 number of servers (Double.POSITIVE_INFINITY for IS)
     * @param scv     MxR squared coefficient of variation of service times
     * @param tol     convergence tolerance (e.g. 1e-6)
     * @param maxiter maximum number of iterations (e.g. 10000)
     * @return throughputs, queue lengths, utilizations, residence times
     */
    public static Result sum_closed(Matrix L, Matrix N, Matrix Z, Matrix mi,
                                    Matrix scv, double tol, int maxiter) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        double K = 0;
        for (int r = 0; r < R; r++) {
            if (Double.isFinite(N.get(r))) {
                K += N.get(r);
            }
        }
        Matrix XN = new Matrix(1, R);
        Matrix QN = new Matrix(M, R);
        Matrix UN = new Matrix(M, R);
        Matrix RN = new Matrix(M, R);
        int it = 0;

        if (K == 0) {
            return new Result(XN, QN, UN, RN, it);
        }

        if (R == 1) {
            // single class: bisection on the throughput (Sec. 9.2.1)
            double lambda_l = 0;
            double lambda_u = Double.POSITIVE_INFINITY;
            for (int i = 0; i < M; i++) {
                if (L.get(i, 0) > 0) {
                    if (Double.isInfinite(mi.get(i))) {
                        lambda_u = Math.min(lambda_u, K / L.get(i, 0));
                    } else {
                        lambda_u = Math.min(lambda_u, mi.get(i) / L.get(i, 0));
                    }
                }
            }
            if (Z.get(0) > 0) {
                lambda_u = Math.min(lambda_u, K / Z.get(0));
            }
            if (Double.isInfinite(lambda_u)) {
                throw new RuntimeException("sum_closed: all service demands are zero.");
            }
            double lambda = lambda_u;
            for (it = 1; it <= maxiter; it++) {
                lambda = (lambda_l + lambda_u) / 2;
                XN.set(0, 0, lambda);
                Matrix Qir = nodeQlen(L, XN, mi, scv, K);
                double g = lambda * Z.get(0) + Qir.sumCols(0);
                if (Math.abs(g - K) <= tol || (lambda_u - lambda_l) <= tol * lambda_u) {
                    break;
                }
                if (g > K) {
                    lambda_u = lambda;
                } else {
                    lambda_l = lambda;
                }
            }
            XN.set(0, 0, lambda);
        } else {
            // multiclass: Gauss-Seidel sweeps of per-class bisections
            for (it = 1; it <= maxiter; it++) {
                double delta = 0;
                for (int r = 0; r < R; r++) {
                    if (N.get(r) == 0) {
                        continue;
                    }
                    double ub = Double.POSITIVE_INFINITY;
                    for (int i = 0; i < M; i++) {
                        if (L.get(i, r) > 0) {
                            if (Double.isInfinite(mi.get(i))) {
                                ub = Math.min(ub, K / L.get(i, r));
                            } else {
                                double rowLoad = 0;
                                for (int q = 0; q < R; q++) {
                                    rowLoad += XN.get(q) * L.get(i, q);
                                }
                                double rem = mi.get(i) - (rowLoad - XN.get(r) * L.get(i, r));
                                ub = Math.min(ub, Math.max(rem, 0) / L.get(i, r));
                            }
                        }
                    }
                    if (Z.get(r) > 0) {
                        ub = Math.min(ub, N.get(r) / Z.get(r));
                    }
                    if (Double.isInfinite(ub)) {
                        throw new RuntimeException("sum_closed: all service demands are zero.");
                    }
                    double lambda_old = XN.get(r);
                    double lambda_l = 0;
                    double lambda_u = ub;
                    while ((lambda_u - lambda_l) > tol * Math.max(ub, 1) / 1e3) {
                        double lambda = (lambda_l + lambda_u) / 2;
                        XN.set(0, r, lambda);
                        Matrix Qir = nodeQlen(L, XN, mi, scv, K);
                        double g = lambda * Z.get(r) + Qir.sumCols(r);
                        if (g > N.get(r)) {
                            lambda_u = lambda;
                        } else {
                            lambda_l = lambda;
                        }
                    }
                    XN.set(0, r, (lambda_l + lambda_u) / 2);
                    delta = Math.max(delta, Math.abs(XN.get(r) - lambda_old));
                }
                if (delta <= tol) {
                    break;
                }
            }
        }

        QN = nodeQlen(L, XN, mi, scv, K);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                if (Double.isInfinite(mi.get(i))) {
                    UN.set(i, r, XN.get(r) * L.get(i, r));
                } else {
                    UN.set(i, r, XN.get(r) * L.get(i, r) / mi.get(i));
                }
                if (XN.get(r) > 0) {
                    RN.set(i, r, QN.get(i, r) / XN.get(r));
                }
            }
        }
        return new Result(XN, QN, UN, RN, it);
    }

    /** per-station per-class mean queue lengths Ki_r = fir(lambda_r) */
    static Matrix nodeQlen(Matrix L, Matrix XN, Matrix mi, Matrix scv, double K) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        Matrix Qir = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            if (Double.isInfinite(mi.get(i))) {
                for (int r = 0; r < R; r++) {
                    Qir.set(i, r, XN.get(r) * L.get(i, r)); // Type 3, Eq. (9.15)
                }
                continue;
            }
            double m = mi.get(i);
            double[] Uir = new double[R];
            double Ui = 0;
            double ci2num = 0;
            for (int r = 0; r < R; r++) {
                Uir[r] = XN.get(r) * L.get(i, r);
                Ui += Uir[r];
                ci2num += Uir[r] * scv.get(i, r);
            }
            if (Ui == 0) {
                continue;
            }
            // per-server utilization; the correction factors keep the node
            // functions finite at rho=1 (Ki(1)<=K)
            double rho = Math.min(Ui / m, 1);
            double ci2 = ci2num / Ui; // demand-weighted node service SCV
            double ai = (1 + ci2) / 2;
            if (K <= m) {
                // never more than m jobs at a m-server node: no queueing
                for (int r = 0; r < R; r++) {
                    Qir.set(i, r, Uir[r]);
                }
                continue;
            }
            if (m == 1) {
                if (ci2 == 1 || K <= 1) {
                    // Type 1,2,4 with mi=1, Eq. (9.15)/(9.19)
                    double den = 1 - (K - 1) / K * rho;
                    for (int r = 0; r < R; r++) {
                        Qir.set(i, r, Uir[r] / den);
                    }
                } else {
                    // -/G/1 FCFS, Eq. (10.88)
                    double den = 1 - (K - 1 - ai) / (K - 1) * rho;
                    for (int r = 0; r < R; r++) {
                        Qir.set(i, r, Uir[r] * (1 + rho * ai / den));
                    }
                }
            } else {
                double Pm = erlangC((int) m, rho);
                if (ci2 == 1) {
                    // Type 1 with mi>1, Eq. (9.15)/(9.19)
                    double den = 1 - (K - m - 1) / (K - m) * rho;
                    for (int r = 0; r < R; r++) {
                        Qir.set(i, r, Uir[r] + (Uir[r] / m) * Pm / den);
                    }
                } else {
                    // -/G/m FCFS, Eq. (10.89)
                    double den = 1 - (K - m - ai) / (K - m) * rho;
                    for (int r = 0; r < R; r++) {
                        Qir.set(i, r, Uir[r] + (Uir[r] / m) * ai * Pm / den);
                    }
                }
            }
        }
        return Qir;
    }

    /** Erlang-C probability of waiting for an M/M/m queue (Eq. 6.28) */
    static double erlangC(int m, double rho) {
        if (rho >= 1) {
            return 1;
        }
        double a = m * rho;
        double s = 0;
        double term = 1; // a^k/k!
        for (int k = 0; k < m; k++) {
            if (k > 0) {
                term *= a / k;
            }
            s += term;
        }
        double last = term * a / m / (1 - rho); // a^m/(m!(1-rho))
        return last / (s + last);
    }
}
