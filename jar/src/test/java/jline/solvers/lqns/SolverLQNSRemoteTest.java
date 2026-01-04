package jline.solvers.lqns;

import jline.VerboseLevel;
import jline.examples.java.basic.LayeredModel;
import jline.lang.layered.LayeredNetwork;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.util.Maths;
import org.junit.jupiter.api.*;

import java.net.HttpURLConnection;
import java.net.URL;

import static jline.TestTools.withSuppressedOutput;
import static org.junit.jupiter.api.Assertions.*;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

/**
 * Test remote LQNS execution via the lqns-rest Docker container.
 *
 * Prerequisites:
 * - Docker image imperialqore/lqns-rest:latest must be available
 * - Container must be running on port 8082 (or update remote_url)
 *
 * To start the container:
 * docker run -p 8082:8080 imperialqore/lqns-rest:latest
 */
@Tag("remote")
public class SolverLQNSRemoteTest {

    private static final String REMOTE_URL = "http://localhost:8082";

    private static final String DOCKER_UNAVAILABLE_WARNING =
        "WARNING: LQNS Docker container not available at " + REMOTE_URL + ". " +
        "To run these tests, start the container with: docker run -p 8082:8080 imperialqore/lqns-rest:latest";

    @BeforeAll
    public static void setUp() {
        Maths.setRandomNumbersMatlab(true);
    }

    /**
     * Check if the LQNS Docker container is available.
     * Prints a warning to stderr if not available.
     * @return true if the container is reachable, false otherwise
     */
    private static boolean isDockerAvailable() {
        try {
            URL url = new URL(REMOTE_URL + "/api/v1/health");
            HttpURLConnection connection = (HttpURLConnection) url.openConnection();
            connection.setRequestMethod("GET");
            connection.setConnectTimeout(3000);
            connection.setReadTimeout(3000);
            int responseCode = connection.getResponseCode();
            connection.disconnect();
            return responseCode >= 200 && responseCode < 400;
        } catch (Exception e) {
            // Try root endpoint as fallback
            try {
                URL url = new URL(REMOTE_URL);
                HttpURLConnection connection = (HttpURLConnection) url.openConnection();
                connection.setRequestMethod("GET");
                connection.setConnectTimeout(3000);
                connection.setReadTimeout(3000);
                int responseCode = connection.getResponseCode();
                connection.disconnect();
                return responseCode >= 200 && responseCode < 500;
            } catch (Exception e2) {
                return false;
            }
        }
    }

    /**
     * Check Docker availability and skip test with warning if not available.
     */
    private void assumeDockerAvailable() {
        boolean available = isDockerAvailable();
        if (!available) {
            System.err.println(DOCKER_UNAVAILABLE_WARNING);
        }
        assumeTrue(available, DOCKER_UNAVAILABLE_WARNING);
    }

    @Test
    @Timeout(60)
    @DisplayName("Remote LQNS execution - simple serial model")
    public void testRemoteLQNSSerial() throws Exception {
        assumeDockerAvailable();

        // Create a simple layered network model
        LayeredNetwork model = LayeredModel.lqn_serial();

        // Configure solver for remote execution
        SolverOptions options = new SolverOptions();
        options.config.remote = true;
        options.config.remote_url = REMOTE_URL;
        options.verbose = VerboseLevel.DEBUG;
        options.method = "lqns";

        // Run solver remotely
        final LayeredNetworkAvgTable[] resultHolder = new LayeredNetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverLQNS solver = new SolverLQNS(model, options);
            resultHolder[0] = solver.getAvgTable();
        });

        LayeredNetworkAvgTable result = resultHolder[0];

        // Verify results are computed
        assertNotNull(result, "Result should not be null");
        assertNotNull(result.getUtil(), "Utilization metrics should be computed");
        assertNotNull(result.getTput(), "Throughput metrics should be computed");

        // Verify table has expected number of elements (2 processors, 2 tasks, 2 entries, 4 activities = 10 total)
        assertEquals(10, result.getNodeNames().size(),
            "Expected 10 nodes in the model");

        // Verify some utilization values are non-zero
        boolean hasNonZeroUtil = result.getUtil().stream()
            .anyMatch(u -> u != null && !Double.isNaN(u) && u > 0.0);
        assertTrue(hasNonZeroUtil, "At least one node should have non-zero utilization");

        // Verify some throughput values are non-zero
        boolean hasNonZeroTput = result.getTput().stream()
            .anyMatch(t -> t != null && !Double.isNaN(t) && t > 0.0);
        assertTrue(hasNonZeroTput, "At least one node should have non-zero throughput");
    }

    @Test
    @Timeout(60)
    @DisplayName("Remote LQSIM execution - simple serial model")
    public void testRemoteLQSIMSerial() throws Exception {
        assumeDockerAvailable();

        // Create a simple layered network model
        LayeredNetwork model = LayeredModel.lqn_serial();

        // Configure solver for remote LQSIM execution
        SolverOptions options = new SolverOptions();
        options.config.remote = true;
        options.config.remote_url = REMOTE_URL;
        options.verbose = VerboseLevel.DEBUG;
        options.method = "lqsim";
        options.samples = 10000;  // Short simulation for testing
        options.seed = 12345;

        // Run solver remotely
        final LayeredNetworkAvgTable[] resultHolder = new LayeredNetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverLQNS solver = new SolverLQNS(model, options);
            resultHolder[0] = solver.getAvgTable();
        });

        LayeredNetworkAvgTable result = resultHolder[0];

        // Verify results are computed
        assertNotNull(result, "Result should not be null");
        assertNotNull(result.getUtil(), "Utilization metrics should be computed");
        assertNotNull(result.getTput(), "Throughput metrics should be computed");

        // Verify table has expected number of elements
        assertEquals(10, result.getNodeNames().size(),
            "Expected 10 nodes in the model");

        // Verify some utilization values are non-zero
        boolean hasNonZeroUtil = result.getUtil().stream()
            .anyMatch(u -> u != null && !Double.isNaN(u) && u > 0.0);
        assertTrue(hasNonZeroUtil, "At least one node should have non-zero utilization");
    }

    @Test
    @Timeout(60)
    @DisplayName("Remote LQNS execution with exactmva method")
    public void testRemoteLQNSExactMVA() throws Exception {
        assumeDockerAvailable();

        // Create a model
        LayeredNetwork model = LayeredModel.lqn_multi_solvers();

        // Configure solver for remote execution with exactmva method
        SolverOptions options = new SolverOptions();
        options.config.remote = true;
        options.config.remote_url = REMOTE_URL;
        options.verbose = VerboseLevel.DEBUG;
        options.method = "exactmva";

        // Run solver remotely
        final LayeredNetworkAvgTable[] resultHolder = new LayeredNetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverLQNS solver = new SolverLQNS(model, options);
            resultHolder[0] = solver.getAvgTable();
        });

        LayeredNetworkAvgTable result = resultHolder[0];

        // Verify results are computed
        assertNotNull(result, "Result should not be null");
        assertTrue(result.getNodeNames().size() > 0, "Should have nodes in result");
    }

    @Test
    @DisplayName("Verify remote configuration defaults")
    public void testRemoteConfigurationDefaults() {
        SolverOptions options = new SolverOptions();

        // Verify default values
        assertFalse(options.config.remote, "Remote should be disabled by default");
        assertEquals("http://localhost:8080", options.config.remote_url,
            "Default remote URL should be http://localhost:8080");
    }

    @Test
    @DisplayName("Verify remote configuration can be modified")
    public void testRemoteConfigurationModification() {
        SolverOptions options = new SolverOptions();

        // Modify remote settings
        options.config.remote = true;
        options.config.remote_url = "http://example.com:9999";

        // Verify changes
        assertTrue(options.config.remote, "Remote should be enabled");
        assertEquals("http://example.com:9999", options.config.remote_url,
            "Remote URL should be updated");
    }
}
