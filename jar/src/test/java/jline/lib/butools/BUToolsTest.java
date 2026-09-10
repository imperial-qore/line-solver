package jline.lib.butools;

import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.ZERO_TOL;
import static jline.TestTools.FINE_TOL;

/**
 * Tests for BUTools library functions.
 */
public class BUToolsTest {


    @Test
    public void testCheckMoments() {
        // Test valid moment sequence for exponential distribution with rate 1
        // Moments are 1!, 2!, 3!, ... = 1, 2, 6, 24, 120
        Matrix validMoments = new Matrix(new double[]{1.0, 2.0, 6.0, 24.0, 120.0});
        assertTrue(CheckMoments.checkMoments(validMoments, 1e-14));
        
        // Test invalid moment sequence
        Matrix invalidMoments = new Matrix(new double[]{1.0, 0.5, 2.0});
        assertFalse(CheckMoments.checkMoments(invalidMoments, 1e-14));
    }

    @Test
    public void testFactorialMomsFromMoms() {
        // Test with simple moments
        Matrix rawMoments = new Matrix(new double[]{1.0, 2.0, 6.0});
        Matrix factorialMoments = FactorialMomsFromMoms.factorialMomsFromMoms(rawMoments);
        
        // For exponential distribution: factorial moments equal raw moments
        assertEquals(3, factorialMoments.length());
        assertEquals(1.0, factorialMoments.get(0), ZERO_TOL);
        assertEquals(1.0, factorialMoments.get(1), ZERO_TOL);
        assertEquals(2.0, factorialMoments.get(2), ZERO_TOL);
    }

    @Test
    public void testReducedMomsConversion() {
        // Test conversion between raw and reduced moments
        Matrix rawMoments = new Matrix(new double[]{1.0, 2.0, 6.0, 24.0});
        
        // Convert to reduced moments
        Matrix reducedMoments = ReducedMomsFromMoms.reducedMomsFromMoms(rawMoments);
        
        // Convert back to raw moments
        Matrix reconstructedMoments = MomsFromReducedMoms.momsFromReducedMoms(reducedMoments);
        
        // Should be identical
        assertEquals(rawMoments.length(), reconstructedMoments.length());
        for (int i = 0; i < rawMoments.length(); i++) {
            assertEquals(rawMoments.get(i), reconstructedMoments.get(i), ZERO_TOL);
        }
    }

    @Test
    public void testHankelMomsConversion() {
        // Test with exponential distribution moments
        Matrix rawMoments = new Matrix(new double[]{1.0, 2.0, 6.0});
        
        // Convert to Hankel moments
        Matrix hankelMoments = HankelMomsFromMoms.hankelMomsFromMoms(rawMoments);
        
        // Convert back to raw moments
        Matrix reconstructedMoments = MomsFromHankelMoms.momsFromHankelMoms(hankelMoments);
        
        // Should be approximately identical (allowing for numerical errors)
        assertEquals(rawMoments.length(), reconstructedMoments.length());
        for (int i = 0; i < rawMoments.length(); i++) {
            assertEquals(rawMoments.get(i), reconstructedMoments.get(i), FINE_TOL);
        }
    }

    @Test
    public void testJointFactorialMomsConversion() {
        // Test with simple 2x2 joint moments matrix
        Matrix jointMoments = new Matrix(new double[][]{
            {1.0, 2.0},
            {2.0, 6.0}
        });
        
        // Convert to joint factorial moments
        Matrix jointFactorialMoments = JFactorialMomsFromJMoms.jFactorialMomsFromJMoms(jointMoments);
        
        // Convert back to joint raw moments
        Matrix reconstructedMoments = JMomsFromJFactorialMoms.jMomsFromJFactorialMoms(jointFactorialMoments);
        
        // Should be approximately identical
        assertEquals(jointMoments.getNumRows(), reconstructedMoments.getNumRows());
        assertEquals(jointMoments.getNumCols(), reconstructedMoments.getNumCols());
        
        for (int i = 0; i < jointMoments.getNumRows(); i++) {
            for (int j = 0; j < jointMoments.getNumCols(); j++) {
                assertEquals(jointMoments.get(i, j), reconstructedMoments.get(i, j), FINE_TOL);
            }
        }
    }
}