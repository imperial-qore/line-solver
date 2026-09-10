/**
 * @file Probabilistic class-oriented method of moments for two-station load-dependent models
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Pfqn_procomom2 {
    private Pfqn_procomom2() {}

    /**
     * Compute marginal state probabilities for the queue in a model consisting of a
     * queueing station and a delay station only.
     */
    public static Ret.pfqnProcomom2 pfqn_procomom2(Matrix L, Matrix N, Matrix Z, Matrix mu, Integer m) {
        int Nsum = (int) N.elementSum();
        int mServers = (m == null) ? 1 : m.intValue();

        Matrix muMatrix;
        if (mu == null || mu.isEmpty()) {
            muMatrix = Matrix.ones(mServers, Nsum + 1);
        } else {
            Matrix muFlat = new Matrix(1, mu.length() + 1);
            muFlat.set(0, 0, 1.0);
            for (int i = 0; i < mu.length(); i++) {
                muFlat.set(0, i + 1, mu.get(i));
            }
            muMatrix = muFlat;
        }

        Matrix p0 = new Matrix(Nsum + 1, 1);
        p0.set(Nsum, 0, 1.0);

        int R = L.getNumCols();
        List<Matrix> T = new ArrayList<Matrix>();
        for (int r = 0; r < R; r++) {
            Matrix Tr = new Matrix(1 + Nsum, 1 + Nsum);
            for (int n = Nsum; n >= 1; n--) {
                int row = Nsum - n;
                Tr.set(row, row, Z.get(r));
                Tr.set(row, row + 1, (n + mServers - 1) * L.get(r) / muMatrix.get(0, n));
            }
            Tr.set(Nsum, Nsum, Z.get(r));
            T.add(Tr);
        }

        Matrix F = Matrix.eye(Nsum + 1);
        Matrix B = Matrix.eye(Nsum + 1);

        for (int r = 0; r < R; r++) {
            int Nr = (int) N.get(r);
            Matrix Tpower = Matrix.pow(T.get(r), Nr);
            F = F.mult(Tpower).scale(1.0 / CombinatoricsUtils.factorial(Nr));
            B = B.mult(T.get(r));
        }

        Matrix pk = F.mult(p0).transpose();
        double G = pk.elementSum();

        double lG;
        if (!Double.isFinite(pk.get(0, 0)) || !Double.isFinite(G)) {
            Matrix logValues = new Matrix(1, pk.getNumCols());
            for (int i = 0; i < pk.getNumCols(); i++) {
                logValues.set(0, i, FastMath.log(pk.get(0, i)));
            }
            lG = Matrix.logsumexp(logValues);
        } else {
            lG = FastMath.log(G);
        }

        pk = pk.scale(1.0 / G);

        Matrix pkReversed = new Matrix(pk.getNumRows(), pk.getNumCols());
        for (int i = 0; i < pk.getNumCols(); i++) {
            pkReversed.set(0, i, pk.get(0, pk.getNumCols() - 1 - i));
        }

        return new Ret.pfqnProcomom2(pkReversed, lG, G, T, F, B);
    }

    public static Ret.pfqnProcomom2 pfqn_procomom2(Matrix L, Matrix N, Matrix Z, Matrix mu) {
        return pfqn_procomom2(L, N, Z, mu, null);
    }

    public static Ret.pfqnProcomom2 pfqn_procomom2(Matrix L, Matrix N, Matrix Z) {
        return pfqn_procomom2(L, N, Z, null, null);
    }
}
