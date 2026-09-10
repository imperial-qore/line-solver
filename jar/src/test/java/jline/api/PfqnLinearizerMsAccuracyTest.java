package jline.api;

import jline.api.pfqn.ld.Pfqn_mvald;
import jline.api.pfqn.mva.Pfqn_linearizerms;
import jline.io.Ret;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Accuracy of the multiserver Linearizer against exact load-dependent MVA.
 *
 * The partially-idle-server term used to carry the raw demand summed over all
 * classes instead of the arriving class's own demand divided by the number of
 * servers, which inflated the residence times by roughly a factor m. The single
 * job case pins the correction: with one job in the network no queueing is
 * possible, so the residence time must equal the service demand exactly.
 */
public class PfqnLinearizerMsAccuracyTest {

    /** Exact reference: LD-MVA with rate min(n, m) at every station. */
    private static Matrix exactResidenceTimes(Matrix L, Matrix N, Matrix Z, Matrix nservers) {
        int M = L.getNumRows();
        int Ntot = (int) N.elementSum();
        Matrix mu = new Matrix(M, Ntot);
        for (int i = 0; i < M; i++) {
            for (int n = 1; n <= Ntot; n++) {
                mu.set(i, n - 1, Math.min(n, nservers.get(i)));
            }
        }
        Ret.pfqnMVALD ret = Pfqn_mvald.pfqn_mvald(L, N, Z, mu);
        double X = ret.X.get(0);
        Matrix W = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            W.set(i, 0, ret.Q.get(i, 0) / X);
        }
        return W;
    }

    private static Matrix col(double... v) {
        Matrix m = new Matrix(v.length, 1);
        for (int i = 0; i < v.length; i++) m.set(i, 0, v[i]);
        return m;
    }

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) m.set(0, i, v[i]);
        return m;
    }

    @Test
    public void singleJobResidenceTimeEqualsDemand() {
        Matrix L = col(1.0, 0.75);
        Matrix Z = row(0.0);
        for (double[] ms : new double[][]{{1, 2}, {1, 3}, {2, 4}, {1, 5}, {3, 3}}) {
            Matrix nservers = col(ms[0], ms[1]);
            Ret.pfqnAMVAMS ret = Pfqn_linearizerms.pfqn_linearizerms(L, row(1.0), Z, nservers);
            for (int i = 0; i < 2; i++) {
                assertEquals(L.get(i, 0), ret.R.get(i, 0), 5e-3,
                        "single-job residence time at station " + i + " with m=" + ms[i]);
            }
        }
    }

    @Test
    public void residenceTimesTrackExactLoadDependentMva() {
        Matrix L = col(1.0, 0.75);
        Matrix Z = row(0.0);
        double worst = 0.0;
        for (double[] ms : new double[][]{{1, 2}, {1, 3}, {2, 4}, {1, 5}, {3, 3}}) {
            Matrix nservers = col(ms[0], ms[1]);
            for (int N = 2; N <= 8; N += 2) {
                Matrix Nv = row((double) N);
                Ret.pfqnAMVAMS amva = Pfqn_linearizerms.pfqn_linearizerms(L, Nv, Z, nservers);
                Matrix exact = exactResidenceTimes(L, Nv, Z, nservers);
                for (int i = 0; i < 2; i++) {
                    worst = Math.max(worst, Math.abs(amva.R.get(i, 0) - exact.get(i, 0)) / exact.get(i, 0));
                }
            }
        }
        // Before the fix the same sweep peaked above 200%.
        assertTrue(worst < 0.05, "worst relative residence-time error was " + (100 * worst) + "%");
    }
}
