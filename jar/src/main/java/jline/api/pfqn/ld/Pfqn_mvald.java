/**
 * Load-Dependent Mean Value Analysis
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.util.PopulationLattice;
import jline.util.matrix.Matrix;

public final class Pfqn_mvald {
    private Pfqn_mvald() {}

    public static Ret.pfqnMVALD pfqn_mvald(Matrix L, Matrix N, Matrix Z, Matrix mu) {
        return pfqn_mvald(L, N, Z, mu, true);
    }

    public static Ret.pfqnMVALD pfqn_mvald(Matrix L, Matrix N, Matrix Z, Matrix mu, boolean stabilize) {
        boolean warn = true;
        boolean isNumStable = true;
        int M = L.getNumRows();
        int R = L.getNumCols();
        double prodN = 1.0;
        for (int i = 0; i < N.getNumRows(); i++) {
            for (int j = 0; j < N.getNumCols(); j++) {
                prodN *= (1 + N.get(i, j));
            }
        }
        Matrix Xs = new Matrix(R, (int) prodN);
        Matrix[] pi = new Matrix[M];
        for (int ist = 0; ist < M; ist++) {
            pi[ist] = Matrix.ones((int) N.elementSum() + 1, (int) prodN);
        }
        Matrix WN = new Matrix(M, R);
        Matrix n = PopulationLattice.pprod(N);
        List<Double> lGN = new ArrayList<Double>();
        lGN.add(Double.valueOf(0.0));
        while (!(n.getNumRows() == 1 && n.getNumCols() == 1 && n.value() == -1.0)) {
            WN = new Matrix(WN.getNumRows(), WN.getNumCols());
            for (int s = 0; s < R; s++) {
                if (n.get(s) > 0) {
                    for (int ist = 0; ist < M; ist++) {
                        WN.set(ist, s, 0);
                        int k = 0;
                        while (k < n.elementSum()) {
                            WN.set(ist, s, WN.get(ist, s)
                                    + ((L.get(ist, s) / mu.get(ist, k)) * (k + 1)
                                    * pi[ist].get(k, PopulationLattice.hashpop(
                                            Matrix.oner(n, new ArrayList<Integer>(Arrays.asList(Integer.valueOf(s)))), N))));
                            k++;
                        }
                    }
                    Xs.set(s, PopulationLattice.hashpop(n, N),
                            n.get(s) / (Z.get(s) + Matrix.extractColumn(WN, s, null).elementSum()));
                }
            }
            int k = 0;
            while (k < n.elementSum()) {
                for (int ist = 0; ist < M; ist++) {
                    pi[ist].set(k + 1, PopulationLattice.hashpop(n, N), 0);
                }
                for (int s = 0; s < R; s++) {
                    if (n.get(s) > 0) {
                        for (int ist = 0; ist < M; ist++) {
                            pi[ist].set(k + 1, PopulationLattice.hashpop(n, N),
                                    pi[ist].get(k + 1, PopulationLattice.hashpop(n, N))
                                            + ((L.get(ist, s) / mu.get(ist, k))
                                            * Xs.get(s, PopulationLattice.hashpop(n, N))
                                            * pi[ist].get(k, PopulationLattice.hashpop(
                                                    Matrix.oner(n, new ArrayList<Integer>(Arrays.asList(Integer.valueOf(s)))), N))));
                        }
                    }
                }
                k++;
            }
            for (int ist = 0; ist < M; ist++) {
                double sumpi = 0.0;
                int kk = 0;
                while (kk < n.elementSum()) {
                    sumpi += pi[ist].get(kk + 1, PopulationLattice.hashpop(n, N));
                    kk++;
                }
                double p0 = 1 - sumpi;
                if (p0 < 0) {
                    if (warn) {
                        System.err.println("pfqn_mvald: MVA-LD is numerically unstable on this model, "
                                + "LINE will force all probabilities to be non-negative.");
                        warn = false;
                        isNumStable = false;
                    }
                    if (stabilize) {
                        pi[ist].set(0, PopulationLattice.hashpop(n, N), FastMath.ulp(1.0));
                    } else {
                        pi[ist].set(0, PopulationLattice.hashpop(n, N), p0);
                    }
                } else {
                    pi[ist].set(0, PopulationLattice.hashpop(n, N), p0);
                }
            }
            Matrix nfind = n.find();
            if (!nfind.isEmpty()) {
                int idx = nfind.getNumRows() - 1;
                while (idx >= 0 && n.get((int) nfind.get(idx)) < 0) {
                    idx--;
                }
                int last_nnz = (int) nfind.get(idx);
                double sumn = 0.0;
                double sumN = 0.0;
                double sumnp = 0.0;
                for (int i = 0; i < R; i++) {
                    if (i < last_nnz) {
                        sumn += n.get(i);
                        sumN += N.get(i);
                    } else if (i > last_nnz) {
                        sumnp += n.get(i);
                    }
                }
                if (sumn == sumN && sumnp == 0.0) {
                    double logX = FastMath.log(Xs.get(last_nnz, PopulationLattice.hashpop(n, N)));
                    lGN.add(Double.valueOf(lGN.get(lGN.size() - 1).doubleValue() - logX));
                }
            }
            n = PopulationLattice.pprod(n, N);
        }
        Matrix XN = Matrix.extractColumn(Xs, PopulationLattice.hashpop(N, N), null).transpose();
        Matrix newpi = new Matrix(M, (int) N.elementSum() + 1);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < (int) N.elementSum() + 1; j++) {
                newpi.set(i, j, pi[i].get(j, PopulationLattice.hashpop(N, N)));
            }
        }
        Matrix QN;
        if (WN.isEmpty()) {
            QN = new Matrix(0, 0);
        } else {
            QN = WN.elementMult(XN.repmat(M, 1), null);
        }
        Matrix UN = Matrix.ones(newpi.getNumRows(), 1).sub(Matrix.extractColumn(newpi, 0, null));
        // see _kb/03-api-layer.md for rationale
        Matrix CN = new Matrix(N.getNumRows(), N.getNumCols());
        for (int i = 0; i < CN.getNumRows(); i++) {
            for (int j = 0; j < CN.getNumCols(); j++) {
                if (N.get(i, j) <= 0) {
                    continue;
                }
                CN.set(i, j, N.get(i, j) / XN.get(i, j) - Z.get(i, j));
            }
        }
        return new Ret.pfqnMVALD(XN, QN, UN, CN, lGN, isNumStable, newpi);
    }
}
