package jline.solvers.lqns;

import jline.VerboseLevel;
import jline.examples.java.basic.LayeredModel;
import jline.lang.layered.LayeredNetwork;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.util.Maths;
import org.junit.jupiter.api.*;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.concurrent.TimeUnit;

import static jline.TestTools.withSuppressedOutput;
import static org.junit.jupiter.api.Assertions.*;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

/**
 * Test remote LQNS execution via the lqns-rest Docker container.
 *
 * <p>When Docker is available, the container is automatically started before tests
 * and stopped afterwards. If the Docker image is not present locally, the system
 * property {@code -Dlqns.docker.pull=true} must be set to authorize pulling it.</p>
 *
 * <p>To manually start the container instead:
 * <pre>docker run -p 8082:8080 imperialqore/lqns-rest:latest</pre></p>
 */
@Tag("remote")
public class SolverLQNSRemoteTest {

    private static final String REMOTE_URL = "http://localhost:8082";
    private static final String DOCKER_IMAGE = "imperialqore/lqns-rest:latest";
    private static final String CONTAINER_NAME = "lqns-rest-test";
    private static final int HOST_PORT = 8082;
    private static final int CONTAINER_PORT = 8080;
    private static final int HEALTH_CHECK_TIMEOUT_SECONDS = 60;
    private static final String PULL_PROPERTY = "lqns.docker.pull";

    /** True if this test class started the Docker container (and should stop it). */
    private static boolean containerStartedByUs = false;

    @BeforeAll
    public static void setUp() {
        Maths.setRandomNumbersMatlab(true);

        // If the service is already reachable, nothing to do
        if (isServiceReachable()) {
            System.out.println("[LQNS Docker] Service already reachable at " + REMOTE_URL);
            return;
        }

        // Check if Docker daemon is available
        if (!isDockerInstalled()) {
            System.err.println("[LQNS Docker] Docker is not installed or not running. " +
                "Remote LQNS tests will be skipped.");
            return;
        }

        // Clean up any stale container from a previous crashed run
        cleanupStaleContainer();

        // Check if the image is available locally
        if (!isImageAvailable()) {
            String pullProp = System.getProperty(PULL_PROPERTY, "false");
            if (!"true".equalsIgnoreCase(pullProp)) {
                System.err.println("[LQNS Docker] Image '" + DOCKER_IMAGE + "' not found locally.\n" +
                    "  To allow automatic download, re-run with: -D" + PULL_PROPERTY + "=true\n" +
                    "  Or pull manually: docker pull " + DOCKER_IMAGE + "\n" +
                    "  Remote LQNS tests will be skipped.");
                return;
            }
            System.out.println("[LQNS Docker] Pulling image " + DOCKER_IMAGE + " ...");
            if (!pullImage()) {
                System.err.println("[LQNS Docker] Failed to pull image. Remote LQNS tests will be skipped.");
                return;
            }
            System.out.println("[LQNS Docker] Image pulled successfully.");
        }

        // Start the container
        System.out.println("[LQNS Docker] Starting container " + CONTAINER_NAME +
            " on port " + HOST_PORT + " ...");
        if (!startContainer()) {
            System.err.println("[LQNS Docker] Failed to start container. Remote LQNS tests will be skipped.");
            return;
        }

        // Wait for the service to become healthy
        System.out.println("[LQNS Docker] Waiting for health check (up to " +
            HEALTH_CHECK_TIMEOUT_SECONDS + "s) ...");
        if (!waitForHealthCheck()) {
            System.err.println("[LQNS Docker] Container started but health check failed. Stopping container.");
            stopContainer();
            return;
        }

        containerStartedByUs = true;
        System.out.println("[LQNS Docker] Container ready at " + REMOTE_URL);
    }

    @AfterAll
    public static void tearDown() {
        if (containerStartedByUs) {
            System.out.println("[LQNS Docker] Stopping container " + CONTAINER_NAME + " ...");
            stopContainer();
            containerStartedByUs = false;
            System.out.println("[LQNS Docker] Container stopped and removed.");
        }
    }

    // ===== Docker lifecycle helpers =====

    /**
     * Check if Docker CLI is installed and the daemon is running.
     */
    private static boolean isDockerInstalled() {
        try {
            Process p = new ProcessBuilder("docker", "info")
                .redirectErrorStream(true).start();
            consumeStream(p);
            boolean finished = p.waitFor(10, TimeUnit.SECONDS);
            return finished && p.exitValue() == 0;
        } catch (Exception e) {
            return false;
        }
    }

    /**
     * Check if the Docker image is available locally.
     */
    private static boolean isImageAvailable() {
        try {
            Process p = new ProcessBuilder("docker", "images", "-q", DOCKER_IMAGE).start();
            String output = readProcessOutput(p);
            boolean finished = p.waitFor(10, TimeUnit.SECONDS);
            return finished && p.exitValue() == 0 && !output.isEmpty();
        } catch (Exception e) {
            return false;
        }
    }

    /**
     * Pull the Docker image from the registry.
     */
    private static boolean pullImage() {
        try {
            Process p = new ProcessBuilder("docker", "pull", DOCKER_IMAGE)
                .inheritIO().start();
            boolean finished = p.waitFor(300, TimeUnit.SECONDS);
            return finished && p.exitValue() == 0;
        } catch (Exception e) {
            return false;
        }
    }

    /**
     * Remove any stale container with our name from a previous run.
     */
    private static void cleanupStaleContainer() {
        try {
            // Check if a container with our name exists (running or stopped)
            Process check = new ProcessBuilder("docker", "container", "inspect", CONTAINER_NAME)
                .redirectErrorStream(true).start();
            consumeStream(check);
            boolean finished = check.waitFor(10, TimeUnit.SECONDS);
            if (finished && check.exitValue() == 0) {
                System.out.println("[LQNS Docker] Removing stale container " + CONTAINER_NAME + " ...");
                Process stop = new ProcessBuilder("docker", "rm", "-f", CONTAINER_NAME)
                    .redirectErrorStream(true).start();
                consumeStream(stop);
                stop.waitFor(15, TimeUnit.SECONDS);
            }
        } catch (Exception e) {
            // Ignore - no stale container
        }
    }

    /**
     * Start the Docker container in detached mode.
     */
    private static boolean startContainer() {
        try {
            Process p = new ProcessBuilder(
                "docker", "run", "-d",
                "--name", CONTAINER_NAME,
                "-p", HOST_PORT + ":" + CONTAINER_PORT,
                DOCKER_IMAGE
            ).redirectErrorStream(true).start();
            consumeStream(p);
            boolean finished = p.waitFor(30, TimeUnit.SECONDS);
            return finished && p.exitValue() == 0;
        } catch (Exception e) {
            return false;
        }
    }

    /**
     * Stop and remove the Docker container.
     */
    private static void stopContainer() {
        try {
            Process p = new ProcessBuilder("docker", "rm", "-f", CONTAINER_NAME)
                .redirectErrorStream(true).start();
            consumeStream(p);
            p.waitFor(15, TimeUnit.SECONDS);
        } catch (Exception e) {
            System.err.println("[LQNS Docker] Warning: failed to stop container: " + e.getMessage());
        }
    }

    /**
     * Wait for the service health check to pass.
     */
    private static boolean waitForHealthCheck() {
        long deadline = System.currentTimeMillis() + HEALTH_CHECK_TIMEOUT_SECONDS * 1000L;
        while (System.currentTimeMillis() < deadline) {
            if (isServiceReachable()) {
                return true;
            }
            try {
                Thread.sleep(1000);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                return false;
            }
        }
        return false;
    }

    /**
     * Check if the REST service is reachable via HTTP.
     */
    private static boolean isServiceReachable() {
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
     * Read all output from a process's stdout.
     */
    private static String readProcessOutput(Process process) throws Exception {
        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
        StringBuilder sb = new StringBuilder();
        String line;
        while ((line = reader.readLine()) != null) {
            sb.append(line).append("\n");
        }
        return sb.toString().trim();
    }

    /**
     * Consume (discard) a process's stdout to prevent blocking.
     */
    private static void consumeStream(final Process process) {
        Thread t = new Thread(new Runnable() {
            @Override
            public void run() {
                try {
                    BufferedReader reader = new BufferedReader(
                        new InputStreamReader(process.getInputStream()));
                    while (reader.readLine() != null) {
                        // discard
                    }
                } catch (Exception e) {
                    // ignore
                }
            }
        });
        t.setDaemon(true);
        t.start();
    }

    // ===== Test guard =====

    private static final String DOCKER_UNAVAILABLE_WARNING =
        "WARNING: LQNS Docker container not available at " + REMOTE_URL + ". " +
        "Ensure Docker is installed and either start the container manually " +
        "(docker run -p 8082:8080 imperialqore/lqns-rest:latest) " +
        "or run with -D" + PULL_PROPERTY + "=true to allow automatic image download.";

    /**
     * Check Docker availability and skip test with warning if not available.
     */
    private void assumeDockerAvailable() {
        boolean available = isServiceReachable();
        if (!available) {
            System.err.println(DOCKER_UNAVAILABLE_WARNING);
        }
        assumeTrue(available, DOCKER_UNAVAILABLE_WARNING);
    }

    // ===== Tests =====

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
