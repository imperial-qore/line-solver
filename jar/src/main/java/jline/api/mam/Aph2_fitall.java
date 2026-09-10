/**
 * @file Absorbing Phase-type distribution comprehensive fitting
 *
 * Fits multiple APH(2) distributions to match given moments with exhaustive parameter search.
 * Provides comprehensive fitting solutions for phase-type distribution approximation.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.util.FastMath;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Aph2_fitall {
    private Aph2_fitall() {}

    /**
     * Fits a set of acyclic phase-type (APH) distributions with two phases (APH(2)) to match the given moments.
     */
    public static Map<Integer, MatrixCell> aph2_fitall(double M1, double M2, double M3) {
        Map<Integer, MatrixCell> APHS = new HashMap<Integer, MatrixCell>();
        double degentol = 1e-8;
        double SCV = (M2 - FastMath.pow(M1, 2)) / FastMath.pow(M1, 2);
        double M3lb = 3 * FastMath.pow(M1, 3) * (3 * SCV - 1 + FastMath.sqrt(2.0) * FastMath.pow(1 - SCV, 1.5));
        double tmp0;
        if (SCV <= 1 && FastMath.abs(M3 - M3lb) < degentol) {
            tmp0 = 0.0;
        } else {
            tmp0 = FastMath.pow(M3, 2) / 9.0 + ((8 * FastMath.pow(M1, 3)) / 3.0 - 2 * M2 * M1) * M3
                    - 3 * FastMath.pow(M1, 2) * FastMath.pow(M2, 2) + 2 * FastMath.pow(M2, 3);
            if (tmp0 < 0) {
                APHS.put(0, Aph_fit.aph_fit(M1, M2, M3, 2));
                return APHS;
            }
        }

        double tmp1 = 3 * FastMath.sqrt(tmp0);
        double tmp2 = M3 - 3 * M1 * M2;
        double tmp3 = (6 * M2 - 12 * FastMath.pow(M1, 2));
        int n = (tmp0 == 0.0) ? 1 : 2;

        Matrix h1v = new Matrix(n, 1, n);
        Matrix h2v = new Matrix(n, 1, n);
        Matrix r1v = new Matrix(n, 1, n);

        if (n == 1) {
            h2v.set(0, 0, tmp2 / tmp3);
            h1v.set(0, 0, tmp2 / tmp3);
        } else {
            h2v.set(0, 0, (tmp2 + tmp1) / tmp3);
            h2v.set(1, 0, (tmp2 - tmp1) / tmp3);
            h1v.set(1, 0, h2v.get(0, 0));
            h1v.set(0, 0, h2v.get(1, 0));
        }

        for (int j = 0; j < n; j++) {
            double h1 = h1v.get(j);
            double h2 = h2v.get(j);
            r1v.set(j, (M1 - h1) / h2);
        }

        int idx = 0;
        for (int j = 0; j < n; j++) {
            double h1 = h1v.get(j);
            double h2 = h2v.get(j);
            double r1 = r1v.get(j);
            if (h1 > 0 && h2 > 0 && r1 >= -degentol && r1 <= (1 + degentol)) {
                r1 = FastMath.max(Math.min(r1, 1.0), 0.0);
                APHS.put(idx, Aph2_assemble.aph2_assemble(h1, h2, r1));
                idx++;
            }
        }

        if (APHS.isEmpty()) {
            APHS.put(0, Aph_fit.aph_fit(M1, M2, M3, 2));
        }

        return APHS;
    }
}
