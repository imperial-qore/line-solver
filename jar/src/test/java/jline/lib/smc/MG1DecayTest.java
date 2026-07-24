/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.lib.smc;

import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static jline.lib.smc.MG1_Decay.mg1_decay;
import static jline.lib.smc.MG1_Shifts.mg1_shifts;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class MG1DecayTest {

    @Test
    public void scalarCase() {
        // A = [a0 a1 a2] with a0+a1+a2=1, drift = a1+2*a2 < 1.
        // For m=1, eta is the non-unit root of a2*z^2 + (a1-1)*z + a0 = 0.
        // By Vieta's (z=1 is one root): eta = a0/a2.
        Matrix A = new Matrix(1, 3);
        A.set(0, 0, 0.5);
        A.set(0, 1, 0.4);
        A.set(0, 2, 0.1);

        MG1DecayResult result = mg1_decay(A, true);
        assertEquals(5.0, result.getEta(), 1e-12, "scalar eta should be a0/a2 = 5.0");
        assertNotNull(result.getUT());
        assertEquals(1, result.getUT().getNumRows());
        assertEquals(1, result.getUT().getNumCols());
    }

    @Test
    public void matrix2x2Case() {
        // 2x2 case verified against MATLAB MG1_Decay: eta = 2.0
        // A0 = [0.3 0.1; 0.2 0.2]
        // A1 = [0.2 0.2; 0.1 0.3]
        // A2 = [0.1 0.1; 0.1 0.1]
        Matrix A = new Matrix(2, 6);
        A.set(0, 0, 0.3); A.set(0, 1, 0.1);
        A.set(1, 0, 0.2); A.set(1, 1, 0.2);
        A.set(0, 2, 0.2); A.set(0, 3, 0.2);
        A.set(1, 2, 0.1); A.set(1, 3, 0.3);
        A.set(0, 4, 0.1); A.set(0, 5, 0.1);
        A.set(1, 4, 0.1); A.set(1, 5, 0.1);

        MG1DecayResult result = mg1_decay(A, true);
        assertEquals(2.0, result.getEta(), 1e-10, "2x2 eta should match MATLAB MG1_Decay = 2.0");
        assertNotNull(result.getUT());
        assertEquals(1, result.getUT().getNumRows());
        assertEquals(2, result.getUT().getNumCols());
    }

    @Test
    public void shiftsTauTypeRunsCleanly() {
        // Smoke test: mg1_shifts with shiftType="tau" and "dbl" (drift<1 branch)
        // should now run without throwing. Use the 2x2 case.
        Matrix A = new Matrix(2, 6);
        A.set(0, 0, 0.3); A.set(0, 1, 0.1);
        A.set(1, 0, 0.2); A.set(1, 1, 0.2);
        A.set(0, 2, 0.2); A.set(0, 3, 0.2);
        A.set(1, 2, 0.1); A.set(1, 3, 0.3);
        A.set(0, 4, 0.1); A.set(0, 5, 0.1);
        A.set(1, 4, 0.1); A.set(1, 5, 0.1);

        Quadruple<Matrix, Double, Double, Matrix> tauResult = mg1_shifts(A, "tau");
        assertEquals(2.0, tauResult.getThird(), 1e-10, "tau in mg1_shifts should equal mg1_decay eta");
        assertTrue(tauResult.getSecond() < 1.0, "drift should be < 1 for this test matrix");

        Quadruple<Matrix, Double, Double, Matrix> dblResult = mg1_shifts(A, "dbl");
        assertEquals(2.0, dblResult.getThird(), 1e-10, "dbl tau in mg1_shifts should also equal mg1_decay eta");
    }
}
