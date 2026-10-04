/**
 * @file M3PP theoretical superposition fitting.
 *
 * Superposes k M3PP to fit the counting-process characteristics of a MMAP[k].
 * Faithful port of MATLAB {@code m3a/m3pp/m3pp_superpos_fitc_theoretical.m}.
 *
 * @since LINE 3.0
 */
package jline.api.mam.m3pp;

import java.util.List;

import jline.api.mam.Mmap_count_idc;
import jline.api.mam.Mmap_count_mean;
import jline.api.mam.Mmap_count_moment;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class M3pp_superpos_fitc_theoretical {
    private M3pp_superpos_fitc_theoretical() {}

    /**
     * Superposes k M3PP to fit the characteristics of {@code MMAP}, using default time
     * scales t = 1, tinf = 1000.
     */
    public static Pair<MatrixCell, List<MatrixCell>> m3pp_superpos_fitc_theoretical(MatrixCell MMAP) {
        return m3pp_superpos_fitc_theoretical(MMAP, 1.0, 1000.0);
    }

    /**
     * Superposes k M3PP to fit the characteristics of a MMAP[k].
     *
     * @param MMAP process to fit
     * @param t    finite time scale
     * @param tinf near-infinite time scale
     * @return Pair of (superposed M3PP[m], list of per-class fitted M3PP components)
     */
    public static Pair<MatrixCell, List<MatrixCell>> m3pp_superpos_fitc_theoretical(
            MatrixCell MMAP, double t, double tinf) {

        // number of classes
        int m = MMAP.size() - 2;

        // per-class rates
        double[] av = Mmap_count_mean.mmap_count_mean(MMAP, 1.0).toArray1D();

        // per-class IDC(t), IDC(inf)
        double[] btv = Mmap_count_idc.mmap_count_idc(MMAP, t).toArray1D();
        double[] binfv = Mmap_count_idc.mmap_count_idc(MMAP, tinf).toArray1D();

        // per-class third central moment of counts at t
        Matrix mtv = Mmap_count_moment.mmap_count_moment(MMAP, t, new int[]{1, 2, 3});
        double[] m3tv = new double[m];
        for (int i = 0; i < m; i++) {
            double m1 = mtv.get(0, i);
            double m2 = mtv.get(1, i);
            double m3 = mtv.get(2, i);
            m3tv[i] = m3 - 3 * m2 * m1 + 2 * m1 * m1 * m1;
        }

        // fit superposition
        return M3pp_superpos_fitc.m3pp_superpos_fitc(av, btv, binfv, m3tv, t, tinf);
    }
}
