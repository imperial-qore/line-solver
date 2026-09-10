/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn.nc;

import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Validation of the Class-Oriented Method of Moments against exact convolution.
 *
 * COMOM is exact for a product-form network, so its log normalizing constant
 * must equal Buzen's to numerical precision. The single-station restriction
 * (M == 1) is inherent to this implementation, so every case below has one
 * queueing station and an arbitrary number of closed classes with think times.
 *
 * These cases exist because pfqn_comom had no coverage at all, and two index
 * defects survived because of it: the CoMoM hash was resolved against the
 * target population rather than the running population vector, which made
 * every basis lookup miss and returned lG = -Inf; and the CE/PC gate summed
 * one class too few, emitting rows MATLAB does not emit.
 */
public class Pfqn_comomTest {

    private static final double TOL = 1e-9;

    /**
     * Reference value cross-checked in MATLAB (pfqn_comom, pfqn_comomrm and
     * pfqn_ca all return this) and against a brute-force sum over the state
     * space: G = 1.842, log G = 0.6108519378.
     */
    @Test
    public void twoClassSingleStationMatchesMatlabReference() {
        Matrix L = new Matrix(1, 2);
        L.set(0, 0, 0.6);
        L.set(0, 1, 0.4);
        Matrix N = new Matrix(1, 2);
        N.set(0, 0, 2);
        N.set(0, 1, 1);
        Matrix Z = new Matrix(1, 2);
        Z.set(0, 0, 1.0);
        Z.set(0, 1, 0.5);

        assertEquals(0.6108519378, Pfqn_comom.pfqn_comom(L, N, Z, 1e-6), 1e-9);
    }

    /** COMOM equals convolution on randomly generated two- to four-class models. */
    @Test
    public void matchesConvolutionOnRandomModels() {
        Random rng = new Random(42);
        for (int trial = 0; trial < 12; trial++) {
            int R = 2 + rng.nextInt(3);
            Matrix L = new Matrix(1, R);
            Matrix N = new Matrix(1, R);
            Matrix Z = new Matrix(1, R);
            for (int r = 0; r < R; r++) {
                L.set(0, r, 0.1 + rng.nextDouble());
                N.set(0, r, 1 + rng.nextInt(3));
                Z.set(0, r, rng.nextDouble());
            }
            assertEquals(Pfqn_ca.pfqn_ca(L, N, Z).lG,
                    Pfqn_comom.pfqn_comom(L, N, Z, 1e-6), TOL,
                    "COMOM disagrees with convolution on trial " + trial);
        }
    }
}
