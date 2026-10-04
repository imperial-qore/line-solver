/**
 * @file Auxiliary EC terms computation for load-dependent mixed MVA
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Pfqn_ldmx_ec {
    private Pfqn_ldmx_ec() {}

    /**
     * Auxiliary function used by pfqn_mvaldmx to compute the EC terms.
     */
    public static Ret.pfqnLDMXEC pfqn_ldmx_ec(Matrix lambda, Matrix D, Matrix mu) {
        int M = mu.getNumRows();
        Matrix Lo = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            Lo.set(i, lambda.mult(Matrix.extractRows(D, i, i + 1, null).transpose()).get(0));
        }
        Matrix b = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            int idx = 0;
            while (idx < mu.getNumCols() && mu.get(i, idx) != mu.get(i, mu.getNumCols() - 1)) {
                idx++;
            }
            b.set(i, (double) idx);
        }
        int Nt = mu.getNumCols();
        int oldEnd = mu.getNumCols() - 1;
        mu.expandMatrix(mu.getNumRows(),
                mu.getNumCols() + 2 + (int) b.elementMax(),
                mu.getNonZeroLength() + (2 + (int) b.elementMax())
                        * Matrix.extractColumn(mu, oldEnd, null).getNonZeroLength());
        for (int i = 0; i < mu.getNumRows(); i++) {
            for (int j = oldEnd + 1; j < mu.getNumCols(); j++) {
                mu.set(i, j, mu.get(i, oldEnd));
            }
        }
        Matrix C = new Matrix(mu.getNumRows(), mu.getNumCols());
        mu.divide(1.0, C, false);
        Matrix EC = new Matrix(M, Nt);
        Matrix E = new Matrix(M, 1 + Nt);
        Matrix Eprime = new Matrix(M, 1 + Nt);
        for (int i = 0; i < M; i++) {
            Matrix E1 = new Matrix(1 + Nt, 1 + Nt);
            Matrix E2 = new Matrix(1 + Nt, 1 + Nt);
            Matrix E3 = new Matrix(1 + Nt, 1 + Nt);
            Matrix F2 = new Matrix(1 + Nt, 2 + (int) b.get(i) - 2);
            Matrix F3 = new Matrix(1 + Nt, 2 + (int) b.get(i) - 2);

            Matrix E2prime = new Matrix(1 + Nt, 1 + Nt);
            Matrix F2prime = new Matrix(1 + Nt, 2 + (int) b.get(i) - 2);
            for (int n = 0; n <= Nt; n++) {
                if (n >= b.get(i) + 1) {
                    E.set(i, n, 1.0 / FastMath.pow(1 - Lo.get(i) * C.get(i, (int) b.get(i)), n + 1));
                    Eprime.set(i, n, C.get(i, (int) b.get(i)) * E.get(i, n));
                } else {
                    if (n == 0) {
                        E1.set(n, 1 / (1 - Lo.get(i) * C.get(i, (int) b.get(i))));
                        int j = 0;
                        while (j < b.get(i) - 1 + 1) {
                            E1.set(n, E1.get(n) * C.get(i, j) / C.get(i, (int) b.get(i)));
                            j++;
                        }
                    } else {
                        E1.set(n, 1 / (1 - Lo.get(i) * C.get(i, (int) b.get(i)))
                                * C.get(i, (int) b.get(i)) / C.get(i, n - 1) * E1.get(n - 1));
                    }

                    {
                        int n0 = 0;
                        while (n0 <= b.get(i) - 2 + 1) {
                            if (n0 == 0) {
                                F2.set(n, n0, 1);
                            } else {
                                F2.set(n, n0, ((double) (n + n0) / n0 * Lo.get(i)
                                        * C.get(i, n + n0 - 1) * F2.get(n, n0 - 1)));
                            }
                            n0++;
                        }
                    }

                    double sumf2 = 0.0;
                    {
                        int k = -1;
                        while (k < b.get(i) - 2 + 1) {
                            sumf2 += F2.get(n, 1 + k);
                            k++;
                        }
                    }
                    E2.set(n, sumf2);

                    {
                        int n0 = 0;
                        while (n0 <= b.get(i) - 2 + 1) {
                            if (n == 0 && n0 == 0) {
                                F3.set(n, n0, 1);
                                int j = 0;
                                while (j < b.get(i) - 1 + 1) {
                                    F3.set(n, n0, F3.get(n, n0) * C.get(i, j) / C.get(i, (int) b.get(i)));
                                    j++;
                                }
                            } else if (n > 0 && n0 == 0) {
                                F3.set(n, n0, C.get(i, (int) b.get(i)) / C.get(i, n - 1) * F3.get(n - 1, 0));
                            } else {
                                F3.set(n, n0, (double) (n + n0) / n0 * Lo.get(i)
                                        * C.get(i, (int) b.get(i)) * F3.get(n, n0 - 1));
                            }
                            n0++;
                        }
                    }

                    double sumf3 = 0.0;
                    {
                        int k = -1;
                        while (k < b.get(i) - 2 + 1) {
                            sumf3 += F3.get(n, 1 + k);
                            k++;
                        }
                    }
                    E3.set(n, sumf3);

                    int n0 = 0;
                    while (n0 <= b.get(i) - 2 + 1) {
                        if (n0 == 0) {
                            F2prime.set(n, n0, C.get(i, n));
                        } else {
                            F2prime.set(n, n0, (double) (n + n0) / n0 * Lo.get(i)
                                    * C.get(i, n + n0) * F2prime.get(n, n0 - 1));
                        }
                        n0++;
                    }

                    double sumf2p = 0.0;
                    int k = -1;
                    while (k < b.get(i) - 2 + 1) {
                        sumf2p += F2prime.get(n, 1 + k);
                        k++;
                    }
                    E2prime.set(n, sumf2p);

                    E.set(i, n, E1.get(n) + E2.get(n) - E3.get(n));
                    if (n < b.get(i) - 1 + 1) {
                        Eprime.set(i, n, C.get(i, (int) b.get(i)) * E1.get(n)
                                + E2prime.get(n) - C.get(i, (int) b.get(i)) * E3.get(n));
                    } else {
                        Eprime.set(i, n, C.get(i, (int) b.get(i)) * E.get(i, n));
                    }
                }
            }
            for (int n = 0; n < Nt; n++) {
                EC.set(i, n, C.get(i, n) * E.get(i, n + 1) / E.get(i, n));
            }
        }
        return new Ret.pfqnLDMXEC(EC, E, Eprime, Lo);
    }
}
