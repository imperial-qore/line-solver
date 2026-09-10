package jline.util;

import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

class MathsTest {

    @Test
    void maxpos_shouldReturnIndexOfMaximumValue() {
        Matrix v = new Matrix(1, 5);
        v.set(0, 0, 1.0);
        v.set(0, 1, 5.0);
        v.set(0, 2, 3.0);
        v.set(0, 3, 7.0);
        v.set(0, 4, 2.0);
        
        int maxIndex = Maths.maxpos(v);
        assertEquals(3, maxIndex);
    }

    @Test
    void maxpos_shouldReturnFirstIndexWhenMultipleMaxValues() {
        Matrix v = new Matrix(1, 4);
        v.set(0, 0, 3.0);
        v.set(0, 1, 7.0);
        v.set(0, 2, 7.0);
        v.set(0, 3, 5.0);
        
        int maxIndex = Maths.maxpos(v);
        assertEquals(1, maxIndex);
    }

    @Test
    void maxpos_shouldHandleSingleElement() {
        Matrix v = new Matrix(1, 1);
        v.set(0, 0, 42.0);
        
        int maxIndex = Maths.maxpos(v);
        assertEquals(0, maxIndex);
    }

    @Test
    void maxpos_shouldHandleNegativeValues() {
        Matrix v = new Matrix(1, 3);
        v.set(0, 0, -5.0);
        v.set(0, 1, -2.0);
        v.set(0, 2, -10.0);
        
        int maxIndex = Maths.maxpos(v);
        assertEquals(1, maxIndex);
    }

    @Test
    void maxpos_shouldHandleColumnVector() {
        Matrix v = new Matrix(4, 1);
        v.set(0, 0, 1.0);
        v.set(1, 0, 8.0);
        v.set(2, 0, 3.0);
        v.set(3, 0, 5.0);
        
        int maxIndex = Maths.maxpos(v);
        assertEquals(1, maxIndex);
    }

    @Test
    void maxpos_shouldThrowExceptionForNullMatrix() {
        assertThrows(IllegalArgumentException.class, () -> Maths.maxpos(null));
    }

    @Test
    void maxpos_shouldThrowExceptionForEmptyMatrix() {
        Matrix v = new Matrix(0, 0);
        assertThrows(IllegalArgumentException.class, () -> Maths.maxpos(v));
    }

    @Test
    void maxpos_withN_shouldReturnTopNIndices() {
        Matrix v = new Matrix(1, 5);
        v.set(0, 0, 1.0);
        v.set(0, 1, 5.0);
        v.set(0, 2, 3.0);
        v.set(0, 3, 7.0);
        v.set(0, 4, 2.0);
        
        int[] topIndices = Maths.maxpos(v, 3);
        assertEquals(3, topIndices.length);
        assertEquals(3, topIndices[0]); // 7.0
        assertEquals(1, topIndices[1]); // 5.0
        assertEquals(2, topIndices[2]); // 3.0
    }

    @Test
    void maxpos_withN_shouldReturnAllIndicesWhenNEqualsLength() {
        Matrix v = new Matrix(1, 3);
        v.set(0, 0, 1.0);
        v.set(0, 1, 3.0);
        v.set(0, 2, 2.0);
        
        int[] topIndices = Maths.maxpos(v, 3);
        assertEquals(3, topIndices.length);
        assertEquals(1, topIndices[0]); // 3.0
        assertEquals(2, topIndices[1]); // 2.0
        assertEquals(0, topIndices[2]); // 1.0
    }

    @Test
    void maxpos_withN_shouldThrowExceptionWhenNIsZero() {
        Matrix v = new Matrix(1, 3);
        v.set(0, 0, 1.0);
        v.set(0, 1, 2.0);
        v.set(0, 2, 3.0);
        
        assertThrows(IllegalArgumentException.class, () -> Maths.maxpos(v, 0));
    }

    @Test
    void maxpos_withN_shouldThrowExceptionWhenNIsNegative() {
        Matrix v = new Matrix(1, 3);
        v.set(0, 0, 1.0);
        v.set(0, 1, 2.0);
        v.set(0, 2, 3.0);
        
        assertThrows(IllegalArgumentException.class, () -> Maths.maxpos(v, -1));
    }

    @Test
    void maxpos_withN_shouldThrowExceptionWhenNIsGreaterThanLength() {
        Matrix v = new Matrix(1, 3);
        v.set(0, 0, 1.0);
        v.set(0, 1, 2.0);
        v.set(0, 2, 3.0);
        
        assertThrows(IllegalArgumentException.class, () -> Maths.maxpos(v, 4));
    }

    @Test
    void probchoose_shouldReturnIndexBasedOnProbability() {
        // Set deterministic seed for reproducible test
        Maths.setMatlabRandomSeed(12345);
        
        List<Double> probs = Arrays.asList(0.1, 0.3, 0.4, 0.2);
        
        // Test multiple times to ensure it's working correctly
        int[] counts = new int[4];
        int numTests = 1000;
        
        for (int i = 0; i < numTests; i++) {
            int choice = Maths.probchoose(probs);
            assertTrue(choice >= 0 && choice < 4);
            counts[choice]++;
        }
        
        // Check that we got some selections for each category
        // (with 1000 tests, even low probability should get some hits)
        for (int i = 0; i < 4; i++) {
            assertTrue(counts[i] > 0, "Category " + i + " should have been selected at least once");
        }
    }

    @Test
    void probchoose_shouldReturnZeroForSingleElement() {
        List<Double> probs = Arrays.asList(1.0);
        
        int choice = Maths.probchoose(probs);
        assertEquals(0, choice);
    }

    @Test
    void probchoose_shouldReturnLastElementWhenAllProbabilitiesAreZero() {
        List<Double> probs = Arrays.asList(0.0, 0.0, 0.0, 0.0);
        
        int choice = Maths.probchoose(probs);
        assertEquals(3, choice);
    }

    @Test
    void probchoose_shouldHandleExtremelySmallProbabilities() {
        List<Double> probs = Arrays.asList(1e-10, 1e-10, 1.0 - 2e-10);
        
        int choice = Maths.probchoose(probs);
        assertTrue(choice >= 0 && choice < 3);
    }

    @Test
    void probchoose_shouldThrowExceptionForNullList() {
        assertThrows(IllegalArgumentException.class, () -> Maths.probchoose((List<Double>) null));
    }

    @Test
    void probchoose_shouldThrowExceptionForEmptyList() {
        List<Double> probs = new ArrayList<Double>();
        assertThrows(IllegalArgumentException.class, () -> Maths.probchoose(probs));
    }

    @Test
    void probchoose_matrix_shouldReturnIndexBasedOnProbability() {
        // Set deterministic seed for reproducible test
        Maths.setMatlabRandomSeed(12345);
        
        Matrix probs = new Matrix(1, 4);
        probs.set(0, 0, 0.1);
        probs.set(0, 1, 0.3);
        probs.set(0, 2, 0.4);
        probs.set(0, 3, 0.2);
        
        // Test multiple times to ensure it's working correctly
        int[] counts = new int[4];
        int numTests = 1000;
        
        for (int i = 0; i < numTests; i++) {
            int choice = Maths.probchoose(probs);
            assertTrue(choice >= 0 && choice < 4);
            counts[choice]++;
        }
        
        // Check that we got some selections for each category
        for (int i = 0; i < 4; i++) {
            assertTrue(counts[i] > 0, "Category " + i + " should have been selected at least once");
        }
    }

    @Test
    void probchoose_matrix_shouldHandleColumnVector() {
        Matrix probs = new Matrix(3, 1);
        probs.set(0, 0, 0.3);
        probs.set(1, 0, 0.4);
        probs.set(2, 0, 0.3);
        
        int choice = Maths.probchoose(probs);
        assertTrue(choice >= 0 && choice < 3);
    }

    @Test
    void probchoose_matrix_shouldThrowExceptionForNullMatrix() {
        assertThrows(IllegalArgumentException.class, () -> Maths.probchoose((Matrix) null));
    }

    @Test
    void probchoose_matrix_shouldThrowExceptionForEmptyMatrix() {
        Matrix probs = new Matrix(0, 0);
        assertThrows(IllegalArgumentException.class, () -> Maths.probchoose(probs));
    }

    @Test
    void probchoose_shouldHandleEdgeCaseWithProbabilitiesThatDontSumToOne() {
        // This tests the robustness of the implementation
        List<Double> probs = Arrays.asList(0.1, 0.2, 0.3); // sum = 0.6

        int choice = Maths.probchoose(probs);
        assertTrue(choice >= 0 && choice < 3);
    }

    @Test
    void factln_shouldNotOverflowForLargeArguments() {
        // log(Gamma.gamma(1+n)) overflowed to Inf for n >= 170; the
        // log-gamma form must stay finite (MATLAB gammaln parity).
        assertEquals(863.2319871924054, Maths.factln(200), 1e-8);
        assertEquals(5912.128178488163, Maths.factln(1000.0), 1e-6);
        assertEquals(Math.log(120.0), Maths.factln(5), 1e-12);
        assertEquals(0.0, Maths.factln(0), 1e-12);
    }
}