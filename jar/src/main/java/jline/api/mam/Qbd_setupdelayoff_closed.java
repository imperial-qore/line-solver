/**
 * @file Finite-population queue with setup delay and server switch-off
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.api.mc.Ctmc_solve;
import jline.GlobalConstants;
import jline.lang.processes.Coxian;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Qbd_setupdelayoff_closed {
    private Qbd_setupdelayoff_closed() {}

    /** Mean queue length and throughput of the closed setup/delay-off queue. */
    public static final class Result {
        /** Mean number of jobs at the station. */
        public final double QN;
        /** Throughput of the station. */
        public final double XN;

        Result(double QN, double XN) {
            this.QN = QN;
            this.XN = XN;
        }
    }

    /**
     * Canonical PH sub-generator of a phase given its rate, entered at phase 1
     * for every SCV. The same helper as {@code Qbd_setupdelayoff}: an
     * exponential phase is built from the rate directly, which avoids the
     * mean-based FineTol cutoff in the fitters, and the Coxian form is what
     * guarantees the entry vector is [1 0 ... 0] so an arrival to an off server
     * enters the setup at phase 1.
     *
     * @param rate the phase rate
     * @param scv the squared coefficient of variation of the phase
     * @return the sub-generator D0 of the phase
     */
    private static Matrix coxianSubgen(double rate, double scv) {
        MatrixCell proc = scv == 1.0 ? new Exp(rate).getProcess()
                                     : Coxian.fitMeanAndSCV(1.0 / rate, scv).getProcess();
        return proc.get(0);
    }

    /** Completion rate of each phase: t(i) = -sum_j D0(i,j). */
    private static double[] completionRates(Matrix D0) {
        int n = D0.getNumRows();
        double[] t = new double[n];
        for (int i = 0; i < n; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < n; j++) {
                rowSum += D0.get(i, j);
            }
            t[i] = -rowSum;
        }
        return t;
    }

    /**
     * Mean queue length and throughput of a FINITE-POPULATION queue with setup
     * delay and delay-off.
     *
     * <p>The closed twin of {@code Qbd_setupdelayoff}. The population N is finite
     * and Z is the complementary delay, the mean time a customer spends away
     * from this station, so the arrival rate is state dependent,
     * lambda(n) = (N - n)/Z, and the level index is bounded by N. That makes the
     * chain a LEVEL-DEPENDENT QBD over finitely many levels, i.e. a finite CTMC,
     * and it is solved exactly rather than by a matrix-geometric tail.
     *
     * <p>THE SEMANTICS ARE THE SIMULATOR'S, not the mean-value shortcut's. When
     * the queue empties the server begins a delay-off period; an arrival DURING
     * it finds the server still warm and resumes without setup (Solver_ssj's
     * cancelDelayoff), and only an arrival after the delay-off has expired pays
     * the setup. That is an M/M/1 with setup time AND close-down time. The
     * per-instance cold-start race {@code p_cold*E[setup] + S} this replaces
     * raced the delay-off against the per-instance idle time and carried NO
     * queueing term, so it described a serverless instance pool rather than a
     * single-server vacation queue and left the reported response time
     * byte-identical across a tenfold change in the setup mean.
     *
     * <p>The phase index is overloaded by level exactly as in the open twin: at
     * level 0 phase 1 is the OFF server and the rest are the delay-off; above
     * level 0 the phases are the setup and the last one is the busy server. Only
     * the REACHABLE states are enumerated, because a finite chain cannot carry an
     * unreachable row: it would be absorbing and the stationary solve singular.
     *
     * @param N population of the closed chain
     * @param Z complementary delay, the mean time a customer spends away
     * @param mu service rate of the station
     * @param alpharate rate of the setup phase
     * @param alphascv squared coefficient of variation of the setup phase
     * @param betarate rate of the delay-off phase
     * @param betascv squared coefficient of variation of the delay-off phase
     * @return the mean queue length and the throughput
     */
    public static Result qbd_setupdelayoff_closed(double N, double Z, double mu, double alpharate,
                                                  double alphascv, double betarate,
                                                  double betascv) {
        int pop = (int) Math.round(N);
        if (pop <= 0 || mu <= 0) {
            return new Result(0.0, 0.0);
        }
        double Zc = Math.max(Z, GlobalConstants.FineTol);

        Matrix Ta = coxianSubgen(alpharate, alphascv);
        int na = Ta.getNumRows();
        double[] ta = completionRates(Ta);
        Matrix Tb = coxianSubgen(betarate, betascv);
        int nb = Tb.getNumRows();
        double[] tb = completionRates(Tb);

        final int off = 0;             // level 0, server off
        final int base = 1 + nb;       // level 0 delay-off occupies 1..nb
        int m = base + pop * (na + 1);
        Matrix Q = Matrix.zeros(m, m);

        double lam0 = pop > 0 ? pop / Zc : 0.0;
        if (lam0 > 0) {
            Q.set(off, setupIdx(base, na, 1, 0), Q.get(off, setupIdx(base, na, 1, 0)) + lam0);
        }
        for (int j = 0; j < nb; j++) {
            int rj = 1 + j;
            for (int j2 = 0; j2 < nb; j2++) {
                if (j2 != j) {
                    Q.set(rj, 1 + j2, Q.get(rj, 1 + j2) + Tb.get(j, j2));
                }
            }
            Q.set(rj, off, Q.get(rj, off) + tb[j]);
            // an arrival during the delay-off cancels it and resumes WITHOUT setup
            if (lam0 > 0) {
                int b1 = busyIdx(base, na, 1);
                Q.set(rj, b1, Q.get(rj, b1) + lam0);
            }
        }
        for (int n = 1; n <= pop; n++) {
            double lam = n < pop ? (pop - n) / Zc : 0.0;
            for (int i = 0; i < na; i++) {
                int ri = setupIdx(base, na, n, i);
                for (int i2 = 0; i2 < na; i2++) {
                    if (i2 != i) {
                        int ci = setupIdx(base, na, n, i2);
                        Q.set(ri, ci, Q.get(ri, ci) + Ta.get(i, i2));
                    }
                }
                int bn = busyIdx(base, na, n);
                Q.set(ri, bn, Q.get(ri, bn) + ta[i]);
                // an arrival during the setup joins the queue and the setup carries
                // on in the SAME phase: the level rises, the phase does not move
                if (lam > 0) {
                    int up = setupIdx(base, na, n + 1, i);
                    Q.set(ri, up, Q.get(ri, up) + lam);
                }
            }
            int bn = busyIdx(base, na, n);
            if (lam > 0) {
                int up = busyIdx(base, na, n + 1);
                Q.set(bn, up, Q.get(bn, up) + lam);
            }
            // a completion that empties the queue starts the delay-off at its phase 1
            int down = n - 1 >= 1 ? busyIdx(base, na, n - 1) : 1;
            Q.set(bn, down, Q.get(bn, down) + mu);
        }
        for (int i = 0; i < m; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < m; j++) {
                if (j != i) {
                    rowSum += Q.get(i, j);
                }
            }
            Q.set(i, i, -rowSum);
        }

        Matrix pi = Ctmc_solve.ctmc_solve(Q);
        double total = 0.0;
        double[] p = new double[m];
        for (int i = 0; i < m; i++) {
            p[i] = Math.max(0.0, pi.get(i));
            total += p[i];
        }
        if (total <= 0) {
            return new Result(0.0, 0.0);
        }
        double QN = 0.0;
        double pbusy = 0.0;
        for (int n = 1; n <= pop; n++) {
            double level = p[busyIdx(base, na, n)] / total;
            for (int i = 0; i < na; i++) {
                level += p[setupIdx(base, na, n, i)] / total;
            }
            QN += n * level;
            pbusy += p[busyIdx(base, na, n)] / total;
        }
        return new Result(QN, mu * pbusy);
    }

    private static int setupIdx(int base, int na, int n, int i) {
        return base + (n - 1) * (na + 1) + i;
    }

    private static int busyIdx(int base, int na, int n) {
        return base + (n - 1) * (na + 1) + na;
    }
}
