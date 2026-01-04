package jline.api;

import jline.util.matrix.MatrixCell;
import org.junit.jupiter.api.Test;

import static jline.TestTools.FINE_TOL;
import static jline.api.mam.Map_exponentialKt.map_exponential;
import static jline.api.mam.Map_erlangKt.map_erlang;
import static jline.api.mam.Map_hyperexpKt.map_hyperexp;
import static jline.api.polling.Polling_qsys_1limitedKt.polling_qsys_1limited;
import static jline.api.polling.Polling_qsys_gatedKt.polling_qsys_gated;
import static jline.api.polling.Polling_qsys_exhaustiveKt.polling_qsys_exhaustive;
import static org.junit.jupiter.api.Assertions.*;

/**
 * Test class for polling system API functions.
 *
 * <p>This class tests the three main polling algorithms:
 * <ul>
 *   <li>polling_qsys_1limited - 1-limited service discipline</li>
 *   <li>polling_qsys_gated - Gated service discipline</li>
 *   <li>polling_qsys_exhaustive - Exhaustive service discipline</li>
 * </ul>
 *
 * <p>Each test uses the same input configuration:
 * <ul>
 *   <li>A = {map_exponential(1), map_erlang(1,2), map_hyperexp(1,10)}</li>
 *   <li>S = {map_erlang(0.1,2), map_erlang(0.1,3), map_hyperexp(0.1,2)}</li>
 *   <li>T = {map_erlang(0.2,2), map_erlang(0.2,2), map_erlang(0.2,2)}</li>
 * </ul>
 */
class PollingAPITest {

    /**
     * Creates test data arrays for arrival, service, and switching time MAPs.
     *
     * @return array containing A, S, T matrices
     */
    private MatrixCell[][] createTestData() {
        MatrixCell[] A = new MatrixCell[3];
        MatrixCell[] S = new MatrixCell[3];
        MatrixCell[] T = new MatrixCell[3];

        // Arrival MAPs: A={map_exponential(1), map_erlang(1,2), map_hyperexp(1,10)}
        A[0] = map_exponential(1.0);
        A[1] = map_erlang(1.0, 2);
        A[2] = map_hyperexp(1.0, 10.0, 0.99);

        // Service MAPs: S={map_erlang(0.1,2), map_erlang(0.1,3), map_hyperexp(0.1,2)}
        S[0] = map_erlang(0.1, 2);
        S[1] = map_erlang(0.1, 3);
        S[2] = map_hyperexp(0.1, 2.0, 0.99);

        // Switching time MAPs: T={map_erlang(0.2,2), map_erlang(0.2,2), map_erlang(0.2,2)}
        T[0] = map_erlang(0.2, 2);
        T[1] = map_erlang(0.2, 2);
        T[2] = map_erlang(0.2, 2);

        return new MatrixCell[][]{A, S, T};
    }

    /**
     * Test polling_qsys_1limited function with expected values:
     * [1.741666666666670, 1.741666666666670, 1.741666666666668]
     */
    @Test
    void test_polling_qsys_1limited() {
        MatrixCell[][] testData = createTestData();
        MatrixCell[] A = testData[0];
        MatrixCell[] S = testData[1];
        MatrixCell[] T = testData[2];

        double[] result = polling_qsys_1limited(A, S, T);

        assertNotNull(result);
        assertEquals(3, result.length);

        assertEquals(1.741666666666670, result[0], FINE_TOL);
        assertEquals(1.741666666666670, result[1], FINE_TOL);
        assertEquals(1.741666666666668, result[2], FINE_TOL);
    }

    /**
     * Test polling_qsys_gated function with expected values:
     * [0.561601459703463, 0.562779821430099, 0.564904433152153]
     */
    @Test
    void test_polling_qsys_gated() {
        MatrixCell[][] testData = createTestData();
        MatrixCell[] A = testData[0];
        MatrixCell[] S = testData[1];
        MatrixCell[] T = testData[2];

        double[] result = polling_qsys_gated(A, S, T);

        assertNotNull(result);
        assertEquals(3, result.length);

        assertEquals(0.561601459703463, result[0], FINE_TOL);
        assertEquals(0.562779821430099, result[1], FINE_TOL);
        assertEquals(0.564904433152153, result[2], FINE_TOL);
    }

    /**
     * Test polling_qsys_exhaustive function with expected values:
     * [0.476962383126767, 0.479587953903023, 0.475592520113068]
     */
    @Test
    void test_polling_qsys_exhaustive() {
        MatrixCell[][] testData = createTestData();
        MatrixCell[] A = testData[0];
        MatrixCell[] S = testData[1];
        MatrixCell[] T = testData[2];

        double[] result = polling_qsys_exhaustive(A, S, T);

        assertNotNull(result);
        assertEquals(3, result.length);

        assertEquals(0.476962383126767, result[0], FINE_TOL);
        assertEquals(0.479587953903023, result[1], FINE_TOL);
        assertEquals(0.475592520113068, result[2], FINE_TOL);
    }
}
