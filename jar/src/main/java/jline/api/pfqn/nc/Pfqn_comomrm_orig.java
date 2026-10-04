/**
 * @file Original CoMoM implementation for finite repairman model.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import jline.GlobalConstants;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

public final class Pfqn_comomrm_orig {
    private Pfqn_comomrm_orig() {}

    public static double pfqn_comomrm_orig(Matrix L, Matrix N, Matrix Z) {
        return pfqn_comomrm_orig(L, N, Z, GlobalConstants.FineTol);
    }

    public static double pfqn_comomrm_orig(Matrix L, Matrix N, Matrix Z, double atol) {
        if (L.getNumRows() != 1) {
            throw new RuntimeException("pfqn_comomrm_orig: The solver accepts at most a single queueing station.");
        }
        int m = 1;
        Matrix lambda = new Matrix(1, N.getNumCols());
        lambda.fill(0.0);
        jline.io.Ret.pfqnNcSanitize ret = Pfqn_nc_sanitize.pfqn_nc_sanitize(lambda, L, N, Z, atol);
        Matrix L_new = ret.L;
        Matrix N_new = ret.N;
        Matrix Z_new = ret.Z;
        double lG0 = ret.lGremaind;

        int M = L_new.getNumRows();
        int R = L_new.getNumCols();

        Matrix Lmax = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            double maxVal = L_new.get(0, r);
            if (maxVal < atol) maxVal = Z_new.get(r);
            Lmax.set(r, maxVal);
        }
        for (int r = 0; r < R; r++) {
            if (Lmax.get(r) > 0) {
                L_new.set(0, r, L_new.get(0, r) / Lmax.get(r));
                Z_new.set(r, Z_new.get(r) / Lmax.get(r));
            }
        }

        int[] rsort = new int[R];
        for (int i = 0; i < R; i++) rsort[i] = i;
        for (int i = 1; i < R; i++) {
            int key = rsort[i];
            double keyVal = Z_new.get(key);
            int j = i - 1;
            while (j >= 0 && Z_new.get(rsort[j]) > keyVal) {
                rsort[j + 1] = rsort[j];
                j--;
            }
            rsort[j + 1] = key;
        }
        Matrix L_sorted = new Matrix(1, R);
        Matrix Z_sorted = new Matrix(1, R);
        Matrix N_sorted = new Matrix(1, R);
        Matrix Lmax_sorted = new Matrix(1, R);
        for (int i = 0; i < R; i++) {
            L_sorted.set(0, i, L_new.get(0, rsort[i]));
            Z_sorted.set(i, Z_new.get(rsort[i]));
            N_sorted.set(i, N_new.get(rsort[i]));
            Lmax_sorted.set(i, Lmax.get(rsort[i]));
        }
        L_new = L_sorted;
        Z_new = Z_sorted;
        N_new = N_sorted;

        Matrix nvec = new Matrix(1, R);
        nvec.fill(0.0);
        Matrix h = new Matrix(2, 1);
        h.set(0, 1.0);
        h.set(1, 1.0);
        Matrix lh = new Matrix(2, 1);
        for (int i = 0; i < 2; i++) {
            lh.set(i, FastMath.log(h.get(i)) + Maths.factln(nvec.elementSum() + M - 1) - Matrix.factln(nvec).elementSum());
        }
        for (int i = 0; i < 2; i++) h.set(i, FastMath.exp(lh.get(i)));

        Matrix scale = new Matrix(1, (int) N_new.elementSum());
        scale.fill(0.0);

        Matrix h_1 = h.copy();
        int nt = 0;

        for (int r = 0; r < R; r++) {
            // see _kb/03-api-layer.md for rationale
            Matrix F1r = null;
            Matrix F2r = null;
            for (int Nr = 1; Nr <= (int) N_new.get(r); Nr++) {
                nvec.set(r, nvec.get(r) + 1);

                if (Nr == 1) {
                    if (r > 0) {
                        int r1 = r;
                        Matrix P = new Matrix(2 * r1, 2 * (r + 1));
                        P.fill(0.0);
                        for (int i = 0; i < r1; i++) P.set(i, i, 1.0);
                        for (int i = 0; i < r1; i++) P.set(r1 + i, r + 1 + i, 1.0);
                        Matrix h1 = new Matrix(2 * (r + 1), 1);
                        h1.fill(0.0);
                        for (int i = 0; i < 2 * (r + 1); i++) {
                            double sum = 0.0;
                            for (int j = 0; j < 2 * r1; j++) sum += P.get(j, i) * h.get(j);
                            h1.set(i, sum);
                        }
                        double nvecSum = nvec.elementSum();
                        if (nvecSum > 1 && nt > 0) {
                            h1.set(r, h_1.get(0) * nvec.get(r1 - 1) / (nvecSum - 1) / scale.get(nt - 1));
                            h1.set(h1.length() - 1, h_1.get(r1) * nvec.get(r1 - 1) / (nvecSum - 1) / scale.get(nt - 1));
                        }
                        h = h1;
                    }
                    int sz = 2 * (r + 1);
                    Matrix A = new Matrix(sz, sz);
                    A.fill(0.0);
                    Matrix B = new Matrix(sz, sz);
                    B.fill(0.0);

                    A.set(0, 0, 1.0);
                    for (int s = 0; s < r; s++) A.set(0, 1 + s, -L_new.get(0, s));
                    A.set(0, r + 1, -1.0);
                    B.set(0, 0, L_new.get(0, r));

                    for (int s = 0; s < r; s++) {
                        A.set(1 + s, r + 1, N_new.get(s));
                        A.set(1 + s, r + 1 + 1 + s, -Z_new.get(s));
                        A.set(1 + s, 1 + s, -m * L_new.get(0, s));
                    }

                    for (int i = 0; i <= r; i++) {
                        A.set(r + 1 + i, r + 1 + i, (double) Nr);
                        B.set(r + 1 + i, i, m * L_new.get(0, r));
                        B.set(r + 1 + i, r + 1 + i, Z_new.get(r));
                    }

                    int rp1 = r + 1;
                    Matrix C = new Matrix(rp1, rp1);
                    Matrix A12 = new Matrix(rp1, rp1);
                    for (int i = 0; i < rp1; i++) {
                        for (int j = 0; j < rp1; j++) {
                            C.set(i, j, A.get(i, j));
                            A12.set(i, j, A.get(i, rp1 + j));
                        }
                    }

                    Matrix B1r = new Matrix(rp1, sz);
                    Matrix B2r = new Matrix(rp1, sz);
                    for (int i = 0; i < rp1; i++) {
                        for (int j = 0; j < sz; j++) {
                            B1r.set(i, j, B.get(i, j));
                            B2r.set(i, j, B.get(rp1 + i, j));
                        }
                    }

                    Matrix iC = C.inv();
                    Matrix iCB1r = iC.mult(B1r);
                    F1r = new Matrix(sz, sz);
                    F1r.fill(0.0);
                    for (int i = 0; i < rp1; i++) {
                        for (int j = 0; j < sz; j++) F1r.set(i, j, iCB1r.get(i, j));
                    }

                    Matrix iCA12 = iC.mult(A12);
                    Matrix iCA12B2r = iCA12.mult(B2r);
                    F2r = new Matrix(sz, sz);
                    F2r.fill(0.0);
                    for (int i = 0; i < rp1; i++) {
                        for (int j = 0; j < sz; j++) F2r.set(i, j, -iCA12B2r.get(i, j));
                    }
                    for (int i = 0; i < rp1; i++) {
                        for (int j = 0; j < sz; j++) F2r.set(rp1 + i, j, B2r.get(i, j));
                    }

                    h_1 = h.copy();
                    Matrix combined = new Matrix(sz, sz);
                    for (int i = 0; i < sz; i++) {
                        for (int j = 0; j < sz; j++) {
                            combined.set(i, j, nvec.get(r) * F1r.get(i, j) + F2r.get(i, j));
                        }
                    }
                    h = combined.mult(h_1);
                    double denom = nvec.elementSum() + M - 1;
                    for (int i = 0; i < h.length(); i++) h.set(i, h.get(i) / denom);
                } else {
                    // see _kb/03-api-layer.md for rationale
                    h_1 = h.copy();
                    int sz = 2 * (r + 1);
                    Matrix combined = new Matrix(sz, sz);
                    for (int i = 0; i < sz; i++) {
                        for (int j = 0; j < sz; j++) {
                            combined.set(i, j, nvec.get(r) * F1r.get(i, j) + F2r.get(i, j));
                        }
                    }
                    h = combined.mult(h_1);
                    double denom = nvec.elementSum() + M - 1;
                    for (int i = 0; i < h.length(); i++) h.set(i, h.get(i) / denom);
                }
                nt = (int) nvec.elementSum();
                scale.set(nt - 1, FastMath.abs(h.sort().elementSum()));
                h.absEq();
                h.scaleEq(1.0 / scale.get(nt - 1));
            }
        }

        double logScale = 0.0;
        for (int i = 0; i < scale.length(); i++) logScale += FastMath.log(scale.get(i));

        double NlogLmax = 0.0;
        for (int i = 0; i < R; i++) NlogLmax += N_sorted.get(i) * FastMath.log(Lmax_sorted.get(i));

        double lG = lG0 + FastMath.log(h.get(h.length() - 1 - (R - 1)))
                + Maths.factln(N_new.elementSum() + M - 1) - Matrix.factln(N_new).elementSum()
                + NlogLmax + logScale;
        return lG;
    }
}
