package jline.api;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;

import org.junit.jupiter.api.Test;

import jline.api.pfqn.mva.Pfqn_bs;
import jline.api.pfqn.mva.Pfqn_chow;
import jline.api.pfqn.mva.Pfqn_clust;
import jline.api.pfqn.mva.Pfqn_dmlin;
import jline.api.pfqn.mva.Pfqn_lcp;
import jline.api.pfqn.mva.Pfqn_linearizer;
import jline.api.pfqn.mva.Pfqn_looping;
import jline.api.pfqn.mva.Pfqn_mva;
import jline.api.pfqn.mva.Pfqn_pam;
import jline.io.Ret;
import jline.lang.constant.SchedStrategy;
import jline.util.matrix.Matrix;

/**
 * The six approximate MVA algorithms of Chapter 2 of H. Wang, "Approximate MVA
 * Algorithms for Solving Queueing Network Models", M.Sc. thesis, University of
 * Toronto, 1997, that LINE lacked until 2026-08-03: Bard LCP, Chow SA, the
 * Hsieh-Lam PAM family, Eager Looping, dSeS-Lavenberg-Muntz Clustering and the
 * dSeS-Muntz Improved Linearizer.
 *
 * <p>The reference values are the MATLAB, Python and C++ ports on the same
 * input; all four agree to the printed digits. Two of the checks pin PROPERTIES
 * rather than numbers, and those are the ones that catch a bad transcription:
 * dmlin MUST equal linearizer (IL is a cost reduction, not an approximation;
 * survey eq. (2.50) drops the Delta^(j) correction and breaks it), and looping
 * MUST bracket the exact solution.</p>
 */
public class PfqnAmvaSurveyTest {

    private static final double TOL = 1e-6;

    private static Matrix demands() {
        Matrix L = new Matrix(3, 2);
        L.set(0, 0, 0.10); L.set(0, 1, 0.30);
        L.set(1, 0, 0.20); L.set(1, 1, 0.05);
        L.set(2, 0, 0.40); L.set(2, 1, 0.15);
        return L;
    }

    private static Matrix pop() {
        Matrix N = new Matrix(1, 2);
        N.set(0, 0, 5);
        N.set(0, 1, 4);
        return N;
    }

    private static Matrix think() {
        Matrix Z = new Matrix(1, 2);
        Z.set(0, 0, 1.0);
        Z.set(0, 1, 2.0);
        return Z;
    }

    /** Max relative throughput error against exact MVA. */
    private static double xerr(Matrix X) {
        Ret.pfqnMVA ex = Pfqn_mva.pfqn_mva(demands(), pop(), think());
        double e = 0.0;
        for (int c = 0; c < X.getNumCols(); c++) {
            e = Math.max(e, Math.abs(X.get(c) - ex.X.get(c)) / ex.X.get(c));
        }
        return e;
    }

    @Test
    public void lcpMatchesTheOtherPorts() {
        Ret.pfqnAMVA r = Pfqn_lcp.pfqn_lcp(demands(), pop(), think());
        assertEquals(1.502598, r.X.get(0), TOL);
        assertEquals(1.188360, r.X.get(1), TOL);
    }

    @Test
    public void lcpIsLessAccurateThanBardSchweitzer() {
        // Table 2.1 ranks LCP below the PE algorithm it seeded: dropping the
        // (N-1)/N factor can only over-count the queue the arrival sees.
        double eLcp = xerr(Pfqn_lcp.pfqn_lcp(demands(), pop(), think()).X);
        double eBs = xerr(Pfqn_bs.pfqn_bs(demands(), pop(), think()).X);
        assertTrue(eLcp > eBs, "LCP (" + eLcp + ") should be worse than BS (" + eBs + ")");
    }

    @Test
    public void chowMatchesTheOtherPortsAndBeatsBardSchweitzer() {
        Ret.pfqnAMVA r = Pfqn_chow.pfqn_chow(demands(), pop(), think());
        assertEquals(1.721975, r.X.get(0), TOL);
        assertEquals(1.250086, r.X.get(1), TOL);
        // Rank 3 against PE's rank 4; a theta-correction of zero would show here.
        assertTrue(xerr(r.X) < xerr(Pfqn_bs.pfqn_bs(demands(), pop(), think()).X));
    }

    @Test
    public void chowBackwardEstimatorIsADifferentAnswer() {
        SchedStrategy[] type = new SchedStrategy[3];
        Arrays.fill(type, SchedStrategy.PS);
        Ret.pfqnAMVA f = Pfqn_chow.pfqn_chow(demands(), pop(), think(), 1e-6, 1000, null, type,
                Pfqn_chow.FORWARD);
        Ret.pfqnAMVA b = Pfqn_chow.pfqn_chow(demands(), pop(), think(), 1e-6, 1000, null, type,
                Pfqn_chow.BACKWARD);
        assertNotEquals(f.X.get(0), b.X.get(0), 1e-9);
    }

    @Test
    public void pamVariants() {
        Ret.pfqnAMVA b = Pfqn_pam.pfqn_pam(demands(), pop(), think(), Pfqn_pam.PAMB);
        assertEquals(1.351351, b.X.get(0), TOL);
        assertEquals(1.024515, b.X.get(1), TOL);
        assertEquals(1, b.totiter);   // noniterative by construction

        // No centre is overloaded here, so the PAMI capping is inert and the
        // two must coincide; a capping applied unconditionally would not.
        Ret.pfqnAMVA i = Pfqn_pam.pfqn_pam(demands(), pop(), think(), Pfqn_pam.PAMI);
        assertEquals(b.X.get(0), i.X.get(0), 1e-12);
        assertEquals(b.X.get(1), i.X.get(1), 1e-12);

        // PAMT unrolls one MVA step more, which moves the answer.
        Ret.pfqnAMVA t = Pfqn_pam.pfqn_pam(demands(), pop(), think(), Pfqn_pam.PAMT);
        assertEquals(1.672982, t.X.get(0), TOL);
        assertEquals(1.197762, t.X.get(1), TOL);
    }

    @Test
    public void dmlinIsLinearizer() {
        // The defining claim of de Souza e Silva and Muntz (1990): same fixed
        // point, lower cost. Transcribing survey eq. (2.50) literally breaks it.
        SchedStrategy[] type = new SchedStrategy[3];
        Arrays.fill(type, SchedStrategy.PS);
        Ret.pfqnAMVA d = Pfqn_dmlin.pfqn_dmlin(demands(), pop(), think());
        Ret.pfqnAMVA l = Pfqn_linearizer.pfqn_linearizer(demands(), pop(), think(), type);
        for (int c = 0; c < d.X.getNumCols(); c++) {
            assertEquals(l.X.get(c), d.X.get(c), 1e-9);
        }
        for (int i = 0; i < d.Q.getNumRows(); i++) {
            for (int c = 0; c < d.Q.getNumCols(); c++) {
                assertEquals(l.Q.get(i, c), d.Q.get(i, c), 1e-9);
            }
        }
        assertEquals(1.790568, d.X.get(0), TOL);
    }

    @Test
    public void clustMatchesTheOtherPorts() {
        Ret.pfqnAMVA r = Pfqn_clust.pfqn_clust(demands(), pop(), think());
        assertEquals(1.799264, r.X.get(0), TOL);
        assertEquals(1.241524, r.X.get(1), TOL);
    }

    @Test
    public void clustWithPeInsideIsADifferentAlgorithm() {
        // Linearizer inside is the only setting that buys anything: CA with the
        // PE algorithm inside every subnetwork IS global PE, per the paper.
        Ret.pfqnAMVA lin = Pfqn_clust.pfqn_clust(demands(), pop(), think(), null, null, "lin",
                1e-6, 1000);
        Ret.pfqnAMVA pe = Pfqn_clust.pfqn_clust(demands(), pop(), think(), null, null, "bs",
                1e-6, 1000);
        assertNotEquals(lin.X.get(0), pe.X.get(0), 1e-9);
    }

    @Test
    public void loopingMatchesTheOtherPorts() {
        Pfqn_looping.Result b = Pfqn_looping.pfqn_looping(demands(), pop(), think());
        assertEquals(1.341502, b.Xlo.get(0), TOL);
        assertEquals(0.987254, b.Xlo.get(1), TOL);
        assertEquals(2.129780, b.Xup.get(0), TOL);
        assertEquals(1.455485, b.Xup.get(1), TOL);
    }

    @Test
    public void loopingBracketsTheExactSolution() {
        // What makes it a bound rather than an estimate.
        Pfqn_looping.Result b = Pfqn_looping.pfqn_looping(demands(), pop(), think());
        Ret.pfqnMVA ex = Pfqn_mva.pfqn_mva(demands(), pop(), think());
        for (int c = 0; c < b.Xlo.getNumCols(); c++) {
            assertTrue(b.Xlo.get(c) <= ex.X.get(c) + 1e-9);
            assertTrue(b.Xup.get(c) >= ex.X.get(c) - 1e-9);
        }
    }
}
