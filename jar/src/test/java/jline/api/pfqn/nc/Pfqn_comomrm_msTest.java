/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn.nc;

import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Multiserver repairman CoMoM.
 *
 * The zero-population case is the one that matters: pfqn_nc_sanitize drops such
 * a class, and both codebases used to loop over the pre-sanitize class count,
 * overrunning the shortened arrays. Fixtures are MATLAB values.
 */
public class Pfqn_comomrm_msTest {

    private static final double TOL = 1e-9;

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) m.set(0, i, v[i]);
        return m;
    }

    /** S=1 degenerates to the single-server repairman shared by the comom family. */
    @Test
    public void singleServerMatchesComomFamily() {
        assertEquals(0.610851937833,
                Pfqn_comomrm_ms.pfqn_comomrm_ms(row(0.6, 0.4), row(2, 1), row(1, 0.5), 1, 1).lG,
                TOL);
    }

    /** S=2 is the same model pfqn_ncld solves with mu = min(1:n, 2). */
    @Test
    public void twoServersMatchLoadDependentNc() {
        assertEquals(0.17227122094,
                Pfqn_comomrm_ms.pfqn_comomrm_ms(row(0.6, 0.4), row(2, 1), row(1, 0.5), 1, 2).lG,
                TOL);
    }

    /** The four-argument form supplies S=1, so it must agree with the explicit call. */
    @Test
    public void fourArgumentFormDefaultsServersToOne() {
        assertEquals(0.610851937833,
                Pfqn_comomrm_ms.pfqn_comomrm_ms(row(0.6, 0.4), row(2, 1), row(1, 0.5), 1).lG,
                TOL);
    }

    /** A zero-population class is dropped by sanitize; the loop must not overrun. */
    @Test
    public void zeroPopulationClassIsDropped() {
        assertEquals(0.603222473032,
                Pfqn_comomrm_ms.pfqn_comomrm_ms(row(0.6, 0.4, 0.5), row(2, 0, 1),
                        row(1, 0.5, 0.3), 1, 1).lG,
                TOL);
    }
}
