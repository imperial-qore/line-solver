/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.advanced;

import jline.examples.java.advanced.PassAndSwapExample;
import jline.solvers.NetworkAvgTable;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * CTMC tests for the pass-and-swap (PAS) / order-independent queue examples.
 *
 * <p>The PAS station has the order-independent product-form stationary
 * distribution (Dorsman & Gardner 2024). Golden QLen/Util/Tput values for the
 * PASQueue station rows were produced by the MATLAB reference implementation and
 * agree across MATLAB, Python-native, Java and Kotlin to CTMC precision. The
 * mirrored {@code PassAndSwapExamplesTest.kt} asserts the same values.
 *
 * @see PassAndSwapExample
 */
public class PassAndSwapExamplesTest {

    private static final double TOL = 1e-6;

    /** Assert the PASQueue station rows (the last {@code nclasses} rows). */
    private static void assertQueueRows(NetworkAvgTable t, int nclasses,
                                        double[] qlen, double[] util, double[] tput) {
        int from = t.getQLen().size() - nclasses;
        assertVec("QLen", t.getQLen(), from, qlen);
        assertVec("Util", t.getUtil(), from, util);
        assertVec("Tput", t.getTput(), from, tput);
    }

    private static void assertVec(String metric, List<Double> actual, int from, double[] expected) {
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], actual.get(from + i), TOL,
                    "PAS CTMC " + metric + "[" + i + "] mismatch");
        }
    }

    @Test
    @DisplayName("PAS CTMC: M/M/K order-independent queue")
    public void testMmk() {
        assertQueueRows(PassAndSwapExample.mmk(), 2,
                new double[]{0.8032786885245902, 0.5737704918032788},
                new double[]{0.3248781568049634, 0.2320558262931842},
                new double[]{0.6497563136099268, 0.4641116525863684});
    }

    @Test
    @DisplayName("PAS CTMC: five-class compatibility (paper Figs 1-2)")
    public void testCompatibility() {
        assertQueueRows(PassAndSwapExample.compatibility(), 5,
                new double[]{0.561490291, 0.4214071462, 0.3053202126, 0.1123506815, 0.1017734042},
                new double[]{0.1303302987, 0.104264239, 0.07819817923, 0.03366135172, 0.02606605974},
                new double[]{0.3909908962, 0.3127927169, 0.2345945377, 0.1563963585, 0.07819817923});
    }

    @Test
    @DisplayName("PAS CTMC: self-loop swapping graph")
    public void testSelfloop() {
        assertQueueRows(PassAndSwapExample.selfloop(), 3,
                new double[]{0.7216685979, 0.4199304751, 0.2940903824},
                new double[]{0.1599073001, 0.1066048667, 0.07995365006},
                new double[]{0.4797219003, 0.3198146002, 0.2398609502});
    }
}
