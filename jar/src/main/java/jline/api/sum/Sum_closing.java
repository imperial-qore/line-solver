package jline.api.sum;

import jline.util.matrix.Matrix;

/**
 * Closing method for open and mixed non-product-form queueing networks
 * (Bolch et al., Sec. 10.1.5), solved with the summation method.
 *
 * The external world of each open class is replaced by an additional
 * -/G/1 station with service rate mu_inf,r = Ropen*lambda0(r), where Ropen
 * is the number of open classes, service SCV equal to the interarrival
 * time SCV of the open class, and unit visit ratio. The resulting closed
 * network is then solved by Sum_closed with a large population Kclosed
 * for the open classes (5000 is recommended for the summation method).
 * Closed classes are passed through unchanged, which makes the method
 * applicable to mixed networks.
 *
 * Reference: G. Bolch, S. Greiner, H. de Meer, K.S. Trivedi, Queueing
 * Networks and Markov Chains, 2nd ed., Wiley, 2006, Sec. 10.1.5.
 */
public class Sum_closing {

    /** Result of the closing method. */
    public static final class Result {
        /** 1xR class throughputs (open classes approach lambda0 from below) */
        public final Matrix XN;
        /** MxR mean queue lengths at the original stations */
        public final Matrix QN;
        /** MxR utilizations at the original stations */
        public final Matrix UN;
        /** MxR residence times at the original stations */
        public final Matrix RN;
        /** 1xR mean response time in the original network, TN=sum(QN)/XN */
        public final Matrix TN;
        /** number of iterations */
        public final int it;

        Result(Matrix XN, Matrix QN, Matrix UN, Matrix RN, Matrix TN, int it) {
            this.XN = XN;
            this.QN = QN;
            this.UN = UN;
            this.RN = RN;
            this.TN = TN;
            this.it = it;
        }
    }

    /**
     * Closing method for open and mixed networks solved with SUM.
     *
     * @param lambda0 1xR external arrival rates (0 for closed classes)
     * @param scva    1xR interarrival time SCVs of the open classes (1 if Poisson)
     * @param L       MxR service demand matrix of the original network, with
     *                visit ratios of open classes normalized per external arrival
     * @param mi      Mx1 number of servers (Double.POSITIVE_INFINITY for IS)
     * @param scv     MxR service time SCVs (pass 1 for insensitive stations)
     * @param N       1xR populations: Double.POSITIVE_INFINITY for open classes,
     *                finite integers for closed classes
     * @param Z       1xR think times
     * @param Kclosed closing population for the open classes (e.g. 5000)
     * @param tol     convergence tolerance (e.g. 1e-6)
     * @param maxiter maximum number of iterations (e.g. 10000)
     * @return throughputs and original-station metrics
     */
    public static Result sum_closing(Matrix lambda0, Matrix scva, Matrix L, Matrix mi,
                                     Matrix scv, Matrix N, Matrix Z, double Kclosed,
                                     double tol, int maxiter) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        int Ropen = 0;
        for (int r = 0; r < R; r++) {
            if (lambda0.get(r) > 0) {
                Ropen++;
            }
        }
        if (Ropen == 0) {
            throw new RuntimeException("sum_closing: no open class, use sum_closed for closed networks.");
        }

        // augment with the closing -/G/1 station, visited by open classes only
        Matrix Laug = new Matrix(M + 1, R);
        Matrix scvaug = new Matrix(M + 1, R);
        Matrix Naug = new Matrix(1, R);
        Matrix miaug = new Matrix(M + 1, 1);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                Laug.set(i, r, L.get(i, r));
                scvaug.set(i, r, scv.get(i, r));
            }
            miaug.set(i, 0, mi.get(i));
        }
        miaug.set(M, 0, 1);
        for (int r = 0; r < R; r++) {
            Naug.set(0, r, N.get(r));
            scvaug.set(M, r, 1);
            if (lambda0.get(r) > 0) {
                Laug.set(M, r, 1 / (Ropen * lambda0.get(r)));
                scvaug.set(M, r, scva.get(r));
                Naug.set(0, r, Kclosed);
            }
        }

        Sum_closed.Result res = Sum_closed.sum_closed(Laug, Naug, Z, miaug, scvaug, tol, maxiter);

        Matrix XN = res.XN;
        Matrix QN = new Matrix(M, R);
        Matrix UN = new Matrix(M, R);
        Matrix RN = new Matrix(M, R);
        Matrix TN = new Matrix(1, R);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                QN.set(i, r, res.QN.get(i, r));
                UN.set(i, r, res.UN.get(i, r));
                RN.set(i, r, res.RN.get(i, r));
            }
        }
        for (int r = 0; r < R; r++) {
            if (XN.get(r) > 0) {
                TN.set(0, r, QN.sumCols(r) / XN.get(r));
            }
        }
        return new Result(XN, QN, UN, RN, TN, res.it);
    }
}
