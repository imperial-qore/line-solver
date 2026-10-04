/**
 * @file M3PP superposition fitting for multiple processes
 *
 * Fits k second-order M3PP processes and superposes them into a higher-order M3PP.
 * Faithful port of MATLAB {@code m3a/m3pp/m3pp_superpos_fitc.m}.
 *
 * @since LINE 3.0
 */
package jline.api.mam.m3pp;

import java.util.ArrayList;
import java.util.List;

import jline.api.mam.Mmap_super;
import jline.api.mam.Mmpp2_fitc;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class M3pp_superpos_fitc {
    private M3pp_superpos_fitc() {}

    /**
     * Fits k second-order M3PP[m_j] and superposes them into a M3PP[m] of order k+1,
     * with m = sum_j m_j.
     *
     * <p>Port of MATLAB {@code m3pp_superpos_fitc.m}: one MMPP(2) is fitted per class via
     * {@link Mmpp2_fitc#mmpp2_fitc}, each is wrapped as a single-class MMAP {D0, D1, D1},
     * and the per-class MMAPs are superposed with {@link Mmap_super#mmap_super}.</p>
     *
     * @param av    per-class rates (length m)
     * @param btv   per-class IDC(t)   (length m)
     * @param binfv per-class IDC(inf) (length m)
     * @param m3tv  per-class third central moment of counts at t (length m)
     * @param t     finite time scale
     * @param tinf  near-infinite time scale
     * @return Pair of (superposed M3PP[m] of order k+1, list of per-class fitted M3PP[m_j])
     */
    public static Pair<MatrixCell, List<MatrixCell>> m3pp_superpos_fitc(double[] av, double[] btv,
                                                                        double[] binfv, double[] m3tv,
                                                                        double t, double tinf) {
        int m = av.length;
        if (btv.length != m || binfv.length != m || m3tv.length != m) {
            throw new IllegalArgumentException("av, btv, binfv, m3tv must all have length m");
        }

        // fit m3pp[2] processes: one MMPP(2) per class, wrapped as single-class MMAP
        List<MatrixCell> m3pps = new ArrayList<MatrixCell>();
        for (int i = 0; i < m; i++) {
            Matrix[] mmpp = Mmpp2_fitc.mmpp2_fitc(av[i], btv[i], btv[i], binfv[i], m3tv[i], t, tinf);
            MatrixCell comp = new MatrixCell(3);
            comp.set(0, mmpp[0]);
            comp.set(1, mmpp[1]);
            comp.set(2, mmpp[1]);
            m3pps.add(comp);
        }

        // perform superposition (fold pairwise, matching mmap_superpos over the cell array)
        MatrixCell fit = m3pps.get(0);
        for (int i = 1; i < m; i++) {
            fit = Mmap_super.mmap_super(fit, m3pps.get(i), "default");
        }

        return new Pair<MatrixCell, List<MatrixCell>>(fit, m3pps);
    }
}
