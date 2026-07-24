/**
 * Maximum Entropy algorithm for Mixed Queueing Networks.
 *
 * Extension of the Kouvatsos (1994) Maximum Entropy Method to mixed
 * open/closed multiclass networks. The 1994 survey notes that the closed
 * two-stage treatment carries over to mixed networks (Section 2.3) but
 * gives no algorithm; this implementation composes the open (Section 3.2)
 * and closed (Section 3.3) ME algorithms by product-form-style
 * conditioning.
 *
 * @since LINE 3.0
 */
package jline.api.nc;

import jline.util.matrix.Matrix;

public final class Me_mqn {
    private Me_mqn() {}

    /**
     * Maximum Entropy algorithm for Mixed Queueing Networks.
     *
     * Composes the open (Section 3.2) and closed (Section 3.3) ME
     * algorithms: (1) the open classes are solved by the open GE-type
     * fixed point on the station set, ignoring the closed classes; (2) the
     * closed classes are solved by the two-stage pseudo-open plus
     * convolution algorithm on servers whose capacity is reduced by the
     * open-class utilization; (3) the open mean queue lengths are inflated
     * by the closed occupancy at single-server stations. Steps 2-3 are
     * exact in the BCMP product-form limit. Only single-server and
     * infinite-server stations are supported.
     *
     * @param M           number of queues (stations)
     * @param R           number of job classes
     * @param openClasses openClasses[r] true when class r is open
     * @param lambda0     external arrival rates [M x R], zero for closed
     * @param Ca0         external arrival scvs [M x R]
     * @param N           class populations [1 x R], infinite for open
     * @param mu          service rates [M x R]
     * @param Cs          service scvs [M x R]
     * @param P           routing probabilities, P[j][i].get(r,0) = p_ji,r
     * @param c           servers per queue [M x 1]; infinity marks IS
     * @param refstat     reference station per class [R], 0-based
     * @param options     algorithm options
     * @return closed-form composition result (X holds the arrival rates
     *         for open classes and the reference-station throughputs for
     *         closed classes)
     */
    public static MeCqnResult me_mqn(int M, int R, boolean[] openClasses,
                                      Matrix lambda0, Matrix Ca0, Matrix N,
                                      Matrix mu, Matrix Cs, Matrix[][] P,
                                      Matrix c, int[] refstat, MeOqnOptions options) {
        return me_mqn(M, R, openClasses, lambda0, Ca0, N, mu, Cs, P, c, refstat, null, options);
    }

    /**
     * Maximum Entropy algorithm for Mixed Queueing Networks with
     * discipline-aware building blocks (insens[i] true for PS/LCFS-PR).
     */
    public static MeCqnResult me_mqn(int M, int R, boolean[] openClasses,
                                      Matrix lambda0, Matrix Ca0, Matrix N,
                                      Matrix mu, Matrix Cs, Matrix[][] P,
                                      Matrix c, int[] refstat, boolean[] insens,
                                      MeOqnOptions options) {
        if (insens == null) {
            insens = new boolean[M];
        }
        int Ro = 0;
        int Rc = 0;
        for (int r = 0; r < R; r++) {
            if (openClasses[r]) {
                Ro++;
            } else {
                Rc++;
            }
        }
        int[] oc = new int[Ro];
        int[] cc = new int[Rc];
        int po = 0;
        int pc = 0;
        for (int r = 0; r < R; r++) {
            if (openClasses[r]) {
                oc[po] = r;
                po++;
            } else {
                cc[pc] = r;
                pc++;
            }
        }

        Matrix L = new Matrix(M, R);
        Matrix W = new Matrix(M, R);
        Matrix Ca = Matrix.ones(M, R);
        Matrix Cd = Matrix.ones(M, R);
        Matrix lambda = new Matrix(M, R);
        Matrix rho = new Matrix(M, R);
        Matrix X = new Matrix(1, R);
        int iters = 0;

        // Step 1: open classes by the Section 3.2 GE-type fixed point
        double[] rho_o = new double[M]; // per-station aggregate open utilization
        if (Ro > 0) {
            Matrix lambda0o = subCols(lambda0, oc);
            Matrix Ca0o = subCols(Ca0, oc);
            Matrix muo = subCols(mu, oc);
            Matrix Cso = subCols(Cs, oc);
            Matrix[][] Po = subRouting(P, oc, M);
            MeOqnResult open = Me_oqn.me_oqn(M, Ro, lambda0o, Ca0o, muo, Cso, Po, c, insens, options);
            iters += open.getIter();
            for (int i = 0; i < M; i++) {
                for (int k = 0; k < Ro; k++) {
                    L.set(i, oc[k], open.getL().get(i, k));
                    Ca.set(i, oc[k], open.getCa().get(i, k));
                    Cd.set(i, oc[k], open.getCd().get(i, k));
                    lambda.set(i, oc[k], open.getLambda().get(i, k));
                    rho.set(i, oc[k], open.getRho().get(i, k));
                }
                if (!Double.isInfinite(c.get(i, 0))) {
                    for (int k = 0; k < Ro; k++) {
                        rho_o[i] += open.getRho().get(i, k);
                    }
                }
            }
            for (int k = 0; k < Ro; k++) {
                double ext = 0.0;
                for (int i = 0; i < M; i++) {
                    ext += lambda0.get(i, oc[k]);
                }
                X.set(0, oc[k], ext);
            }
        }

        // Step 2: closed classes by the Section 3.3 algorithm on servers
        // with capacity reduced by the open-class utilization
        if (Rc > 0) {
            Matrix Nc = new Matrix(1, Rc);
            int[] refc = new int[Rc];
            for (int k = 0; k < Rc; k++) {
                Nc.set(0, k, N.get(0, cc[k]));
                refc[k] = refstat[cc[k]];
            }
            Matrix muc = subCols(mu, cc);
            for (int i = 0; i < M; i++) {
                if (!Double.isInfinite(c.get(i, 0))) {
                    double fac = Math.max(1.0 - rho_o[i], 0.0);
                    for (int k = 0; k < Rc; k++) {
                        muc.set(i, k, muc.get(i, k) * fac);
                    }
                }
            }
            Matrix Csc = subCols(Cs, cc);
            Matrix[][] Pc = subRouting(P, cc, M);
            MeCqnResult closed = Me_cqn.me_cqn(M, Rc, Nc, muc, Csc, Pc, c, refc, insens, options);
            iters += closed.getIter();
            for (int i = 0; i < M; i++) {
                double fac = Double.isInfinite(c.get(i, 0)) ? 1.0 : Math.max(1.0 - rho_o[i], 0.0);
                for (int k = 0; k < Rc; k++) {
                    L.set(i, cc[k], closed.getL().get(i, k));
                    W.set(i, cc[k], closed.getW().get(i, k));
                    Ca.set(i, cc[k], closed.getCa().get(i, k));
                    Cd.set(i, cc[k], closed.getCd().get(i, k));
                    lambda.set(i, cc[k], closed.getLambda().get(i, k));
                    // Closed utilizations are relative to the reduced
                    // capacity; rescale to the physical busy fraction
                    if (Double.isInfinite(c.get(i, 0))) {
                        rho.set(i, cc[k], closed.getRho().get(i, k));
                    } else {
                        rho.set(i, cc[k], closed.getRho().get(i, k) * fac);
                    }
                }
            }
            for (int k = 0; k < Rc; k++) {
                X.set(0, cc[k], closed.getX().get(0, k));
            }
        }

        // see _kb/03-api-layer.md for rationale
        if (Ro > 0 && Rc > 0) {
            for (int i = 0; i < M; i++) {
                if (!Double.isInfinite(c.get(i, 0))) {
                    double Lc_i = 0.0;
                    for (int k = 0; k < Rc; k++) {
                        Lc_i += L.get(i, cc[k]);
                    }
                    for (int k = 0; k < Ro; k++) {
                        L.set(i, oc[k], L.get(i, oc[k]) * (1.0 + Lc_i));
                    }
                }
            }
        }
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < Ro; k++) {
                if (lambda.get(i, oc[k]) > 0) {
                    W.set(i, oc[k], L.get(i, oc[k]) / lambda.get(i, oc[k]));
                }
            }
        }

        return new MeCqnResult(L, W, Ca, Cd, lambda, rho, X, iters);
    }

    /** Extracts the given columns of a matrix. */
    private static Matrix subCols(Matrix A, int[] cols) {
        Matrix B = new Matrix(A.getNumRows(), cols.length);
        for (int i = 0; i < A.getNumRows(); i++) {
            for (int k = 0; k < cols.length; k++) {
                B.set(i, k, A.get(i, cols[k]));
            }
        }
        return B;
    }

    /** Extracts the class slices of a routing array. */
    private static Matrix[][] subRouting(Matrix[][] P, int[] cols, int M) {
        Matrix[][] Q = new Matrix[M][M];
        for (int j = 0; j < M; j++) {
            for (int i = 0; i < M; i++) {
                Q[j][i] = new Matrix(cols.length, 1);
                for (int k = 0; k < cols.length; k++) {
                    Q[j][i].set(k, 0, P[j][i].get(cols[k], 0));
                }
            }
        }
        return Q;
    }
}
