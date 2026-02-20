/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lib.lti;

import org.junit.jupiter.api.Test;
import java.math.BigDecimal;
import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for LTI (Laplace Transform Inversion) library functions.
 */
public class LTITest {

    @Test
    public void testGetomegaHandlesOddN() {
        // Test that the algorithm handles odd n by converting to even n
        BigDecimal[] result_odd = gaverstehfest.INSTANCE.getomega(11);  // Odd n
        BigDecimal[] result_even = gaverstehfest.INSTANCE.getomega(10); // Even n
        
        // Both should have the same length (10) since odd is decremented
        assertEquals(result_odd.length, result_even.length);
        assertEquals(10, result_odd.length);
    }

    @Test
    public void testGetalphaNeedstEvenN() {
        // Test that getalpha also handles odd n correctly
        BigDecimal[] result_odd = gaverstehfest.INSTANCE.getalpha(11);  // Odd n
        BigDecimal[] result_even = gaverstehfest.INSTANCE.getalpha(10); // Even n
        
        // Both should have the same length (10) since odd is decremented
        assertEquals(result_odd.length, result_even.length);
        assertEquals(10, result_odd.length);
    }
    
    @Test
    public void testAlgorithmConsistency() {
        // Test that both functions handle parity consistently
        BigDecimal[] omega_12 = gaverstehfest.INSTANCE.getomega(12);
        BigDecimal[] alpha_12 = gaverstehfest.INSTANCE.getalpha(12);
        
        BigDecimal[] omega_11 = gaverstehfest.INSTANCE.getomega(11); // Should behave like 10
        BigDecimal[] alpha_11 = gaverstehfest.INSTANCE.getalpha(11); // Should behave like 10
        
        BigDecimal[] omega_10 = gaverstehfest.INSTANCE.getomega(10);
        BigDecimal[] alpha_10 = gaverstehfest.INSTANCE.getalpha(10);
        
        // Verify even n works correctly
        assertEquals(12, omega_12.length);
        assertEquals(12, alpha_12.length);
        
        // Verify odd n is converted to even
        assertEquals(10, omega_11.length);
        assertEquals(10, alpha_11.length);
        assertEquals(10, omega_10.length);
        assertEquals(10, alpha_10.length);
        
        // Verify odd n produces same results as corresponding even n
        assertEquals(omega_11.length, omega_10.length);
        assertEquals(alpha_11.length, alpha_10.length);
    }
}