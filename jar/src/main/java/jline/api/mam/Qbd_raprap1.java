/**
 * @file Quasi-Birth-Death process RAP/RAP/1 queue analysis
 *
 * <p>The two RAPs are INDEPENDENT of each other, so the QBD phase space is the
 * product of the two phase spaces and the blocks factor as Kronecker products.
 * This class is the thin product-space wrapper; the analysis itself is the
 * block-level core {@link Qbd_rap}, which solves an arbitrary QBD with RAP
 * components. A model whose arrival process and sequence of service times share
 * a phase space, and are therefore cross-correlated, has no such product
 * structure and must call {@link Qbd_rap} directly.</p>
 *
 * <p>References: N. G. Bean and B. F. Nielsen, "Quasi-Birth-and-Death Processes
 * with Rational Arrival Process Components", Stochastic Models, 26(3), 2010,
 * pp. 309-334. The analysis rests on the prediction-process interpretation of a
 * RAP due to Asmussen and Bladt, which is what allows a QBD argument to be
 * carried over to matrices that are not nonnegative. The same prediction
 * process underlies the conditional-vector RAP sampler in {@link Rap_sample}.</p>
 *
 * <p>Phase ordering: the QBD phase is the pair (arrival phase, service phase)
 * laid out with the ARRIVAL phase major and the service phase minor, i.e. the
 * phase index is a*ns + s. That is the ordering produced by
 * {@code RAPa.kron(I_ns)} and {@code I_na.kron(RAPs)}, and pqueue is indexed by
 * it downstream, so the two Kronecker factors must not be swapped.</p>
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Qbd_raprap1 {
    private Qbd_raprap1() {}

    /**
     * Cap on the number of reported level vectors. The level series is cut at
     * the same point in MATLAB, the JAR and Python, otherwise the mean queue
     * lengths disagree in the truncated tail.
     */
    private static final int MAX_NUM_COMP = 100;

    /** Mass threshold at which the level series is truncated. */
    private static final double MASS_TOL = 1e-10;

    public static QbdRapRap1Result qbd_raprap1(MatrixCell RAPa, MatrixCell RAPs, Double util) {
        int na = RAPa.get(0).getNumRows();
        int ns = RAPs.get(0).getNumRows();

        MatrixCell scaledRAPs = RAPs;
        if (util != null) {
            double lambdaScale = Map_lambda.map_lambda(RAPa);
            scaledRAPs = Map_scale.map_scale(RAPs, util.doubleValue() / lambdaScale);
        }

        double lambdaA = Map_lambda.map_lambda(RAPa);

        // QBD blocks of the RAP/RAP/1 queue: a level is the number in system, a
        // phase is the (arrival,service) RAP phase pair, arrival phase major.
        // The level rises on an arrival (F), falls on a service completion (B),
        // and the two RAPs evolve independently between level changes (L, the
        // Kronecker sum of the two hidden generators). At level 0 the queue is
        // empty, so no service completion can occur and only the arrival RAP
        // evolves: the boundary local block is B1 = kron(Ca,I) and the boundary
        // up-block is F.
        Matrix IA = Matrix.eye(na);
        Matrix IS = Matrix.eye(ns);
        Matrix F = RAPa.get(1).kron(IS);
        Matrix L = RAPa.get(0).kron(IS).add(IA.kron(scaledRAPs.get(0)));
        Matrix B = IA.kron(scaledRAPs.get(1));
        Matrix B1 = RAPa.get(0).kron(IS);

        // Theorem 7 of Bean and Nielsen (2010): G from the quadratic matrix
        // equation, U = L + F*G, R = F*inv(-U), and pi0 from the boundary
        // equation pi0*(B1 + R*B) = 0 normalised by pi0*inv(I-R)*e = 1.
        QbdRapResult core = Qbd_rap.qbd_rap(F, L, B, F, B1, 0);
        Matrix R = core.getR();
        Matrix G = core.getG();
        Matrix pi0 = core.getPi0();

        // Level series, truncated when the accumulated mass reaches 1-1e-10 or
        // at MAX_NUM_COMP level vectors. Qbd_rap returns the exact mean queue
        // length in closed form, but this method keeps reporting the truncated
        // series because pqueue is the documented return value and the three
        // codebases must cut the tail at exactly the same point.
        int n = na * ns;
        java.util.List<Matrix> levels = new java.util.ArrayList<Matrix>();
        levels.add(pi0);
        double sumpi = rowSum(pi0);
        int numit = 1;
        while (sumpi < 1.0 - MASS_TOL && numit < 1 + MAX_NUM_COMP) {
            Matrix next = levels.get(numit - 1).mult(R);
            levels.add(next);
            numit++;
            sumpi += rowSum(next);
        }

        int numLevels = levels.size();
        Matrix pqueue = new Matrix(numLevels, n);
        for (int i = 0; i < numLevels; i++) {
            Matrix lvl = levels.get(i);
            for (int j = 0; j < n; j++) {
                pqueue.set(i, j, lvl.get(0, j));
            }
        }

        double etaVal = core.getSpr();
        Matrix eta = new Matrix(1, 1);
        eta.set(0, 0, etaVal);

        double XN = lambdaA;
        double UN;
        double QN;

        if (na == 1 && ns == 1) {
            UN = 1.0 - pqueue.get(0, 0);
        } else {
            double sum0 = 0.0;
            for (int j = 0; j < n; j++) {
                sum0 += pqueue.get(0, j);
            }
            UN = 1.0 - sum0;
        }

        QN = 0.0;
        for (int i = 0; i < numLevels; i++) {
            double levelProb = 0.0;
            for (int j = 0; j < n; j++) {
                levelProb += pqueue.get(i, j);
            }
            QN += (double) i * levelProb;
        }

        return new QbdRapRap1Result(XN, QN, UN, pqueue, R, eta, G, B, L, F);
    }

    public static QbdRapRap1Result qbd_raprap1(MatrixCell RAPa, MatrixCell RAPs) {
        return qbd_raprap1(RAPa, RAPs, null);
    }

    private static double rowSum(Matrix row) {
        double s = 0.0;
        for (int j = 0; j < row.getNumCols(); j++) {
            s += row.get(0, j);
        }
        return s;
    }
}
