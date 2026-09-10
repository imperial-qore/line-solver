package jline.api;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;

import org.junit.jupiter.api.Test;

import jline.api.pfqn.mva.Pfqn_bs;
import jline.api.pfqn.mva.Pfqn_cntol;
import jline.api.pfqn.mva.Pfqn_linearizer;
import jline.io.Ret;
import jline.lang.constant.SchedStrategy;
import jline.util.matrix.Matrix;

/**
 * Chandy-Neuse population-scaled termination test (Pfqn_cntol).
 *
 * The cutoff 1/(4000 + 16*sum(N)) and the metric max_{i,r}|dQ(i,r)|/N_r it is compared
 * against are published in K. M. Chandy, D. Neuse, Commun. ACM 25(2):126-134, 1982, p.129
 * and appendix; LQNS runs the same expression from SchweitzerCommon. It is opt-in here:
 * the sentinel is a NaN tolerance, and the default tolerances are untouched, so the
 * figures below double as a guard that the sentinel path is the only thing it changes.
 *
 * The asserted values are MATLAB's, printed by pfqn_linearizer/pfqn_bs with tol = 'cn'
 * on the same model.
 */
public class PfqnCntolTest {

    private static Matrix demands() {
        Matrix L = new Matrix(3, 2);
        L.set(0, 0, 1.0);  L.set(0, 1, 2.0);
        L.set(1, 0, 3.0);  L.set(1, 1, 1.0);
        L.set(2, 0, 0.5);  L.set(2, 1, 0.5);
        return L;
    }

    private static Matrix population() {
        Matrix N = new Matrix(1, 2);
        N.set(0, 8.0);
        N.set(1, 1.0);
        return N;
    }

    @Test
    public void cutoffMatchesThePublishedExpression() {
        assertEquals(1.0 / (4000.0 + 160.0), Pfqn_cntol.pfqn_cntol(10.0), 1e-15);
        assertEquals(1.0 / (4000.0 + 16.0 * 9.0), Pfqn_cntol.pfqn_cntol(population()), 1e-15);
        // Below 0.00025 even at a population of one, as the paper's appendix states.
        assertTrue(Pfqn_cntol.pfqn_cntol(1.0) < 0.00025);
        // Decreasing in the population, which is the point of the scaling.
        assertTrue(Pfqn_cntol.pfqn_cntol(100.0) < Pfqn_cntol.pfqn_cntol(10.0));
    }

    @Test
    public void nanIsTheSentinelAndNothingElseIs() {
        assertTrue(Pfqn_cntol.isCntol(Double.NaN));
        assertFalse(Pfqn_cntol.isCntol(1e-8));
        assertFalse(Pfqn_cntol.isCntol(0.0));
    }

    @Test
    public void linearizerUnderTheChandyNeuseTestMatchesMatlab() {
        Matrix L = demands();
        Matrix N = population();
        Matrix Z = new Matrix(1, 2);
        SchedStrategy[] type = new SchedStrategy[3];
        Arrays.fill(type, SchedStrategy.PS);

        Ret.pfqnAMVA cn = Pfqn_linearizer.pfqn_linearizer(L, N, Z, type, Double.NaN, 1000);
        assertEquals(0.304718171165, cn.X.get(0), 1e-9);
        assertEquals(0.084038879640, cn.X.get(1), 1e-9);
        assertEquals(0.555262542605, cn.Q.get(0, 0), 1e-9);
        assertEquals(7.255702544745, cn.Q.get(1, 0), 1e-9);
        assertEquals(0.189034912650, cn.Q.get(2, 0), 1e-9);
        assertEquals(0.251901005767, cn.Q.get(0, 1), 1e-9);
        assertEquals(0.697720815306, cn.Q.get(1, 1), 1e-9);
        assertEquals(0.050378178927, cn.Q.get(2, 1), 1e-9);

        // The looser published cutoff must stop earlier than the 1e-8 default, and the two
        // answers must differ: an identical result would mean the sentinel did nothing.
        Ret.pfqnAMVA def = Pfqn_linearizer.pfqn_linearizer(L, N, Z, type, 1e-8, 1000);
        assertTrue(cn.totiter < def.totiter);
        assertTrue(Math.abs(cn.X.get(0) - def.X.get(0)) > 1e-12);
        assertEquals(def.X.get(0), cn.X.get(0), 1e-4);
    }

    @Test
    public void bardSchweitzerUnderTheChandyNeuseTestMatchesMatlab() {
        Matrix L = demands();
        Matrix N = population();
        Matrix Z = new Matrix(1, 2);

        Ret.pfqnAMVA cn = Pfqn_bs.pfqn_bs(L, N, Z, Double.NaN, 1000);
        assertEquals(0.301066244185, cn.X.get(0), 1e-9);
        assertEquals(0.083878140075, cn.X.get(1), 1e-9);

        Ret.pfqnAMVA def = Pfqn_bs.pfqn_bs(L, N, Z, 1e-6, 1000);
        assertTrue(cn.totiter < def.totiter);
        assertEquals(def.X.get(0), cn.X.get(0), 1e-4);
    }
}
