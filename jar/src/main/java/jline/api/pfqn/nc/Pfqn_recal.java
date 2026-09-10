/**
 * @file RECAL normalizing constant method for closed queueing networks
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import org.apache.commons.math3.util.FastMath;

import jline.api.pfqn.Pfqn_replicas;
import jline.io.Ret;
import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Pfqn_recal {
    private Pfqn_recal() {}

    /**
     * RECAL method to compute the normalizing constant of a load-independent closed queueing network model.
     */
    public static Ret.pfqnNc pfqn_recal(Matrix L, Matrix N) {
        int R = L.getNumCols();
        Matrix Z = new Matrix(1, R);
        return pfqn_recal(L, N, Z);
    }

    public static Ret.pfqnNc pfqn_recal(Matrix L, Matrix N, Matrix Z) {
        int M = L.getNumRows();
        Matrix m0 = new Matrix(1, M);
        m0.fill(1.0);
        return pfqn_recal(L, N, Z, m0);
    }

    public static Ret.pfqnNc pfqn_recal(Matrix L, Matrix N, Matrix Z, Matrix m0) {
        int R = L.getNumCols();

        Matrix m0_local = m0;
        if (m0_local.getNumRows() > 1) {
            m0_local = m0_local.transpose();
        }

        jline.api.pfqn.PfqnUniqueResult uniqueResult = Pfqn_replicas.pfqn_unique(L);
        Matrix L_reduced = uniqueResult.getL_unique();
        int[] mapping = uniqueResult.getMapping();
        int M = L_reduced.getNumRows();

        m0_local = Pfqn_replicas.pfqn_combine_mi(m0_local, mapping, M);

        int Ntot = (int) N.elementSum();

        int G_1_size = (int) (Maths.nCk((double) (Ntot + (M + 1) - 1), (double) Ntot) + 0.5);

        Matrix G_1 = new Matrix(1, G_1_size);
        G_1.fill(1.0);
        Matrix G = G_1.copy();

        Matrix I_1 = Maths.multichoose((double) (M + 1), (double) Ntot);
        int n = 0;

        for (int r = 0; r < R; r++) {
            for (int nr = 1; nr <= (int) N.get(0, r); nr++) {
                n++;
                Matrix I = Maths.multichoose((double) (M + 1), (double) (Ntot + 1 - (n + 1)));

                G = new Matrix(1, I.getNumRows());

                for (int i = 0; i < I.getNumRows(); i++) {
                    Matrix m = I.getRow(i);
                    Matrix mZ = Matrix.extractColumns(m, 0, M);
                    Matrix I_1_subset = Matrix.extractColumns(I_1, 0, M);
                    int mzIndex = Matrix.matchrow(I_1_subset, mZ);

                    G.set(i, Z.get(0, r) * G_1.get(mzIndex) / nr);

                    for (int j = 0; j < M; j++) {
                        m.set(j, m.get(j) + 1);
                        int fullIndex = Matrix.matchrow(I_1, m);
                        G.set(i, G.get(i) + (m.get(j) + m0_local.get(0, j) - 1) * L_reduced.get(j, r) * G_1.get(fullIndex) / nr);
                        m.set(j, m.get(j) - 1);
                    }
                }
                I_1 = I;
                G_1 = G;
            }
        }

        return new Ret.pfqnNc(G.get(0), FastMath.log(G.get(0)));
    }
}
