package jline.api.fj;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

/**
 * The k-of-n (quorum) fork-join tail.
 *
 * <p>A quorum join fires on the k-th of n siblings, so the request response time is the k-th
 * ORDER STATISTIC of the branch times and not their maximum. Reading the maximum there returns
 * the AND-join tail under a quorum's name -- the same number for every k, which is what the
 * forktail route did before this class existed.</p>
 *
 * <p>The reference values are agreed across MATLAB, native python and the C++ port, and were
 * checked against a 400k-sample Monte Carlo of the SAME fitted GE branches (so the check is of
 * the order-statistic inversion, not of the GE fit): every k agreed to within 1.2% at the 99th
 * percentile, which is the sampling error there.</p>
 */
public class FJTailOrdstatTest {

    private static final double REL = 1e-5;
    private static final double ET_HOM = 2.0;
    private static final double VT_HOM = 6.0;
    private static final int N_HOM = 4;
    private static final double[] ET_HET = {1.0, 2.0, 4.0};
    private static final double[] VT_HET = {1.0, 8.0, 40.0};

    @Test
    public void homogeneousReference() {
        double[][] ref = {{0.871483, 2.179109}, {2.146981, 4.220359},
                          {4.189345, 7.427549}, {8.718015, 15.050364}};
        for (int k = 1; k <= N_HOM; k++) {
            assertEquals(ref[k - 1][0],
                    FJ_tail_ordstat.fj_tail_ordstat(ET_HOM, VT_HOM, N_HOM, 0.90, k),
                    ref[k - 1][0] * REL, "k=" + k + " p=0.90");
            assertEquals(ref[k - 1][1],
                    FJ_tail_ordstat.fj_tail_ordstat(ET_HOM, VT_HOM, N_HOM, 0.99, k),
                    ref[k - 1][1] * REL, "k=" + k + " p=0.99");
        }
    }

    @Test
    public void heterogeneousReference() {
        double[][] ref = {{0.980116, 2.397990}, {3.006171, 7.388476}, {12.421746, 29.881069}};
        for (int k = 1; k <= 3; k++) {
            assertEquals(ref[k - 1][0], FJ_tail_ordstat.fj_tail_ordstat(ET_HET, VT_HET, 0.90, k),
                    ref[k - 1][0] * REL, "k=" + k + " p=0.90");
            assertEquals(ref[k - 1][1], FJ_tail_ordstat.fj_tail_ordstat(ET_HET, VT_HET, 0.99, k),
                    ref[k - 1][1] * REL, "k=" + k + " p=0.99");
        }
    }

    @Test
    public void fullJoinIsForkTailExactly() {
        // No existing result may move: a full join must return the ForkTail root itself.
        double[] ps = {0.5, 0.9, 0.99, 0.999};
        for (double p : ps) {
            assertEquals(FJ_tail_forktail.fj_tail_forktail(ET_HOM, VT_HOM, N_HOM, p),
                    FJ_tail_ordstat.fj_tail_ordstat(ET_HOM, VT_HOM, N_HOM, p, N_HOM), 0.0,
                    "homogeneous p=" + p);
            assertEquals(FJ_tail_forktail.fj_tail_forktail(ET_HET, VT_HET, p),
                    FJ_tail_ordstat.fj_tail_ordstat(ET_HET, VT_HET, p, 3), 0.0,
                    "heterogeneous p=" + p);
        }
    }

    @Test
    public void percentileGrowsWithTheQuorum() {
        // Waiting for more siblings can only take longer.
        double[] ps = {0.5, 0.9, 0.99};
        for (double p : ps) {
            for (int k = 1; k < N_HOM; k++) {
                assertTrue(FJ_tail_ordstat.fj_tail_ordstat(ET_HOM, VT_HOM, N_HOM, p, k)
                        < FJ_tail_ordstat.fj_tail_ordstat(ET_HOM, VT_HOM, N_HOM, p, k + 1),
                        "homogeneous p=" + p + " k=" + k);
            }
            for (int k = 1; k < 3; k++) {
                assertTrue(FJ_tail_ordstat.fj_tail_ordstat(ET_HET, VT_HET, p, k)
                        < FJ_tail_ordstat.fj_tail_ordstat(ET_HET, VT_HET, p, k + 1),
                        "heterogeneous p=" + p + " k=" + k);
            }
        }
    }

    @Test
    public void outOfRangeQuorumOrPercentileIsRefused() {
        assertThrows(IllegalArgumentException.class,
                () -> FJ_tail_ordstat.fj_tail_ordstat(ET_HET, VT_HET, 0.99, 0));
        assertThrows(IllegalArgumentException.class,
                () -> FJ_tail_ordstat.fj_tail_ordstat(ET_HET, VT_HET, 0.99, 4));
        assertThrows(IllegalArgumentException.class,
                () -> FJ_tail_ordstat.fj_tail_ordstat(ET_HET, VT_HET, 0.0, 2));
        assertThrows(IllegalArgumentException.class,
                () -> FJ_tail_ordstat.fj_tail_ordstat(ET_HET, VT_HET, 1.0, 2));
    }
}
