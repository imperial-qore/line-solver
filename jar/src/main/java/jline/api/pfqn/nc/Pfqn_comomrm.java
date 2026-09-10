/**
 * @file Class-Oriented Method of Moments for Repairman models (COMOM-RM)
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Pfqn_comomrm {
    private Pfqn_comomrm() {}

    /**
     * Compute the normalizing constant of a repairmen model using COMOM
     */
    public static Ret.pfqnComomrm pfqn_comomrm(Matrix L, Matrix N, Matrix Z, Integer m, double atol) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        if (M != 1) {
            throw new RuntimeException("pfqn_comomrm: The solver accepts at most a single queueing station.");
        }
        if (m == null) {
            m = 1;
        }
        Matrix lambda = N.copy();
        lambda.fill(0.0);
        Ret.pfqnNcSanitize ret = Pfqn_nc_sanitize.pfqn_nc_sanitize(lambda, L, N, Z, atol);
        Matrix L_new = ret.L;
        Matrix N_new = ret.N;
        Matrix Z_new = ret.Z;
        double lG0 = ret.lGremaind;
        // see _kb/03-api-layer.md for rationale
        R = L_new.getNumCols();
        List<Integer> zerothinktimes = new ArrayList<Integer>();

        for (int i = 0; i < Z_new.length(); i++) {
            if (Z_new.get(i) < GlobalConstants.FineTol) {
                zerothinktimes.add(i);
            }
        }
        Matrix nvec = new Matrix(1, R);
        nvec.fill(0.0);

        Matrix lh;

        if (!zerothinktimes.isEmpty()) {
            for (int i = 0; i < zerothinktimes.size(); i++) {
                nvec.set(zerothinktimes.get(i), N_new.get(zerothinktimes.get(i)));
            }
            lh = new Matrix(0, 1);
            Matrix tmp = new Matrix(1, 1);
            tmp.set(0, Maths.factln(nvec.elementSum() + m) - Matrix.factln(nvec).elementSum());
            lh = Matrix.concatRows(lh, tmp, null);

            // see _kb/03-api-layer.md for rationale
            for (int i = 0; i < zerothinktimes.size(); i++) {
                int s = zerothinktimes.get(i);
                Matrix nvec_s = nvec.copy();
                nvec_s.set(s, nvec_s.get(s) - 1);
                tmp.set(0, Maths.factln(nvec_s.elementSum() + m) - Matrix.factln(nvec_s).elementSum());
                lh = Matrix.concatRows(lh, tmp, null);
            }
            tmp.set(0, Maths.factln(nvec.elementSum() + m - 1) - Matrix.factln(nvec).elementSum());
            lh = Matrix.concatRows(lh, tmp, null);
            // Same oner index as above (pfqn_comomrm.m:66-67).
            for (int i = 0; i < zerothinktimes.size(); i++) {
                int s = zerothinktimes.get(i);
                Matrix nvec_s = nvec.copy();
                nvec_s.set(s, nvec_s.get(s) - 1);
                tmp.set(0, Maths.factln(nvec_s.elementSum() + m - 1) - Matrix.factln(nvec_s).elementSum());
                lh = Matrix.concatRows(lh, tmp, null);
            }
        } else {
            lh = new Matrix(2, 1);
            lh.fill(0.0);
        }
        Matrix h = lh.copy();
        for (int i = 0; i < h.length(); i++) {
            h.set(i, FastMath.exp(h.get(i)));
        }

        double lG;
        Matrix lGbasis;

        if (zerothinktimes.size() == R) {
            lGbasis = h.copy();
            for (int i = 0; i < lGbasis.length(); i++) {
                lGbasis.set(i, FastMath.log(lGbasis.get(i)));
            }
            lG = lG0 + FastMath.log(h.get(h.length() - 1 - R));
        } else {
            Matrix scale = new Matrix(1, (int) N_new.elementSum());
            scale.fill(1.0);
            double nt = nvec.elementSum();
            Matrix h_1 = h.copy();
            for (int r = zerothinktimes.size() + 1; r <= R; r++) {
                Matrix F1r = null;
                Matrix F2r = null;
                int Nr = 1;
                while (Nr <= N_new.get(r - 1)) {
                    nvec.set(r - 1, nvec.get(r - 1) + 1);
                    if (Nr == 1) {
                        if (r > zerothinktimes.size() + 1) {
                            Matrix hr = new Matrix(2 * r, 1);
                            hr.fill(0.0);
                            for (int i = 0; i < r - 1; i++) {
                                hr.set(i, h.get(i));
                            }
                            for (int i = r; i < 2 * r - 1; i++) {
                                hr.set(i, h.get(i - 1));
                            }
                            h = hr;
                            if (nt > 0) {
                                h.set(r - 1, h_1.get(0) / scale.get((int) nt - 1));
                                h.set(h.length() - 1, h_1.get(r - 1) / scale.get((int) nt - 1));
                            }
                        }

                        Matrix A12 = new Matrix(r, r);
                        A12.fill(0.0);
                        A12.set(0, 0, -1);
                        for (int s = 1; s < r; s++) {
                            A12.set(s, 0, N_new.get(s - 1));
                            A12.set(s, s, -Z_new.get(s - 1));
                        }

                        Matrix B2r = Matrix.eye(r);
                        Matrix B2r_tmp = Matrix.eye(r);
                        for (int i = 0; i < r; i++) {
                            B2r.set(i, i, m * L_new.get(0, r - 1));
                            B2r_tmp.set(i, i, Z_new.get(r - 1));
                        }
                        B2r = Matrix.concatColumns(B2r, B2r_tmp, null);

                        Matrix iC = Matrix.eye(r);
                        for (int i = 0; i < r; i++) {
                            iC.set(i, i, 1.0 / m);
                            iC.set(0, i, 1.0 / m);
                        }
                        iC.set(0, 0, -1);

                        F1r = new Matrix(2 * r, 2 * r);
                        F1r.fill(0.0);
                        F1r.set(0, 0, 1);

                        F2r = iC.mult(A12).mult(B2r);
                        F2r = Matrix.concatRows(F2r, B2r, null);
                    }
                    h_1 = h;
                    Matrix tmp_mat = F1r.copy();
                    for (int i = 0; i < tmp_mat.getNumRows(); i++) {
                        for (int j = 0; j < tmp_mat.getNumCols(); j++) {
                            tmp_mat.set(i, j, F1r.get(i, j) + F2r.get(i, j) / nvec.get(r - 1));
                        }
                    }
                    h = tmp_mat.mult(h_1);
                    nt = nvec.elementSum();
                    scale.set((int) nt - 1, FastMath.abs(h.elementSum()));
                    for (int i = 0; i < h.length(); i++) {
                        h.set(i, FastMath.abs(h.get(i)) / scale.get((int) nt - 1));
                    }
                    Nr++;
                }
            }
            Matrix log_scale = scale.copy();
            for (int i = 0; i < log_scale.length(); i++) {
                log_scale.set(i, FastMath.log(log_scale.get(i)));
            }
            lG = lG0 + FastMath.log(h.get(h.length() - 1 - (R - 1))) + log_scale.elementSum();
            lGbasis = h.copy();
            for (int i = 0; i < lGbasis.length(); i++) {
                lGbasis.set(i, FastMath.log(h.get(i)) + log_scale.elementSum());
            }
        }
        return new Ret.pfqnComomrm(lG, lGbasis);
    }
}
