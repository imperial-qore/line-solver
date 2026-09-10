package jline.api.sym;

import com.google.gson.JsonObject;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.ServerSocket;
import java.util.ArrayList;
import java.util.List;

/**
 * Resolves the symbolic backend to use, and owns the container that serves it.
 *
 * <p>Resolution order, the same in MATLAB (SAGE.m) and Python
 * (line_solver.api.sym):</p>
 * <ol>
 *   <li>an explicit URL, from solver options or the {@code requested} argument;</li>
 *   <li>the {@code LINE_SAGE_URL} environment variable;</li>
 *   <li>a line-sage-rest service already listening on a conventional port;</li>
 *   <li>a container started here from a locally present image;</li>
 *   <li>nothing, in which case the caller falls back to whatever native
 *       algebra it has, or reports that no backend is configured.</li>
 * </ol>
 *
 * <p>Step 3 verifies identity through {@code /api/v1/info} rather than trusting
 * the port: every imperialqore line-*-rest service listens on 8080 by
 * convention, so a health probe alone would happily accept the LQNS service.</p>
 *
 * <p>The container started in step 4 is reused for the life of the JVM and
 * stopped by a shutdown hook. It is bound to an ephemeral host port, so
 * several JVMs, or a JVM alongside a hand-started service, do not collide.</p>
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public class SymEngines {

    /** Image serving the symbolic REST API. */
    public static final String DOCKER_IMAGE = "imperialqore/line-sage-rest:latest";
    /** Fallback tags, tried in order after {@link #DOCKER_IMAGE}. */
    public static final String[] DOCKER_IMAGE_CANDIDATES = {
            "imperialqore/line-sage-rest:latest", "imperialqore/line-sage-rest"};
    /** Environment variable naming a service to use. */
    public static final String URL_ENV = "LINE_SAGE_URL";
    /** Ports probed for an already running service, in order. */
    public static final int[] PROBE_PORTS = {8085, 8080};
    /** Seconds to wait for a container to report healthy. */
    public static final int STARTUP_TIMEOUT_SECONDS = 120;

    private static SageRestEngine started;
    private static String startedContainer;

    private SymEngines() {
    }

    /**
     * Resolves an engine with the default request, i.e. "auto".
     *
     * @return an engine, or null if no backend could be resolved
     */
    public static SymEngine resolve() {
        return resolve("auto");
    }

    /**
     * Resolves an engine.
     *
     * @param requested "" or "auto" to search, a URL to use a specific service,
     *                  "none" to disable the backend, or an image name to start
     * @return an engine, or null if no backend could be resolved
     */
    public static synchronized SymEngine resolve(String requested) {
        String req = requested == null ? "auto" : requested.trim();
        if ("none".equalsIgnoreCase(req) || "off".equalsIgnoreCase(req)) {
            return null;
        }
        if (req.startsWith("http://") || req.startsWith("https://")) {
            SageRestEngine engine = new SageRestEngine(req);
            return engine.isAvailable() && engine.isUsable() ? engine : null;
        }

        String env = System.getenv(URL_ENV);
        if (env != null && env.trim().length() > 0) {
            SageRestEngine engine = new SageRestEngine(env.trim());
            if (engine.isAvailable() && engine.isUsable()) {
                return engine;
            }
        }

        if (started != null && started.isAvailable() && started.isUsable()) {
            return started;
        }

        for (int i = 0; i < PROBE_PORTS.length; i++) {
            SageRestEngine engine = new SageRestEngine("http://localhost:" + PROBE_PORTS[i]);
            if (isSageService(engine) && engine.isUsable()) {
                return engine;
            }
        }

        // A container this JVM did not start, on the EPHEMERAL port startContainer
        // gives it. PROBE_PORTS only names the two conventional ones, so a healthy
        // line-sage-rest published on a free port was invisible and a second
        // container was started beside it -- and since a JVM that dies without
        // running its shutdown hook never stops the first, they accumulate: two
        // were found on picard05 on 2026-09-09, up SEVEN DAYS and still serving.
        // Asking docker what is already listening reuses it instead. NOT stopped
        // by stopContainer: startedContainer stays null, because a container this
        // process did not start is not this process's to remove.
        SageRestEngine running = runningContainer();
        if (running != null) {
            return running;
        }

        // see _kb/03-api-layer.md for rationale (sym/ section)
        boolean search = req.length() == 0 || "auto".equalsIgnoreCase(req)
                || "true".equalsIgnoreCase(req) || "sage".equalsIgnoreCase(req);
        String image = search ? findImage() : req;
        // Auto-pull only on an explicit opt-in: the "sage" keyword or a named
        // image. Bare "auto"/"true"/"" keep the native backend unless the image
        // is already local, so leaving symbolic on auto never triggers a pull.
        if (image == null && "sage".equalsIgnoreCase(req)) {
            image = pullImage(DOCKER_IMAGE);
        } else if (!search && image != null && !jline.io.DockerImage.hasLocalImage(image)) {
            image = pullImage(image);
        }
        if (image == null) {
            return null;
        }
        try {
            return startContainer(image);
        } catch (IOException e) {
            return null;
        }
    }

    /**
     * A line-sage-rest container already running on this host, or null.
     *
     * <p>Reads the published host port out of {@code docker ps} and verifies the
     * service the same way the port probe does, so a container that is up but
     * unhealthy (or is some other line-*-rest) is not returned.</p>
     *
     * @return a usable engine for an already-running container, or null
     */
    private static SageRestEngine runningContainer() {
        String out = run(15, "docker", "ps", "--filter", "name=line-sage-rest-",
                "--format", "{{.Ports}}");
        if (out == null) {
            return null;
        }
        String[] lines = out.split("\n");
        for (int i = 0; i < lines.length; i++) {
            // e.g. "8888/tcp, 0.0.0.0:38572->8080/tcp, [::]:38572->8080/tcp"
            String[] parts = lines[i].split(",");
            for (int j = 0; j < parts.length; j++) {
                String part = parts[j].trim();
                int arrow = part.indexOf("->8080/tcp");
                if (arrow < 0) {
                    continue;
                }
                int colon = part.lastIndexOf(':', arrow);
                if (colon < 0) {
                    continue;
                }
                String port = part.substring(colon + 1, arrow).trim();
                try {
                    Integer.parseInt(port);
                } catch (NumberFormatException e) {
                    continue;
                }
                SageRestEngine engine = new SageRestEngine("http://localhost:" + port);
                if (isSageService(engine) && engine.isUsable()) {
                    return engine;
                }
            }
        }
        return null;
    }

    /**
     * Checks that a service is line-sage-rest and not another line-*-rest
     * service on the same port.
     *
     * @param engine the candidate
     * @return true if the service identifies as the symbolic API
     */
    private static boolean isSageService(SageRestEngine engine) {
        try {
            JsonObject info = engine.info();
            return info.has("sage_version");
        } catch (IOException e) {
            return false;
        }
    }

    /**
     * Pull {@code target} if the Docker storage location has room; return the
     * tag on success, else null. Storage-guarded via {@link jline.io.DockerImage}
     * (same guard as the LQNS/QNS/JMT wrappers), so an opt-in symbolic request
     * never silently fills the Docker disk; on refusal the caller keeps its
     * native algebra.
     */
    private static String pullImage(String target) {
        if (!jline.io.DockerImage.hasStorageFor(target)) {
            System.err.println("[LINE] Skipping docker pull of " + target
                    + ": insufficient free space at the Docker storage location; "
                    + "keeping the native symbolic backend.");
            return null;
        }
        System.out.println("[LINE] Pulling Docker image " + target + " (this may take a while)...");
        if (jline.io.DockerImage.pull(target) && jline.io.DockerImage.hasLocalImage(target)) {
            return target;
        }
        return null;
    }

    /**
     * @return the first locally present image tag, or null if none is
     */
    public static String findImage() {
        for (int i = 0; i < DOCKER_IMAGE_CANDIDATES.length; i++) {
            String out = run(30, "docker", "images", "-q", DOCKER_IMAGE_CANDIDATES[i]);
            if (out != null && out.trim().length() > 0) {
                return DOCKER_IMAGE_CANDIDATES[i];
            }
        }
        return null;
    }

    /**
     * Starts the service in a container and waits for it to report healthy.
     *
     * @param image the image to run
     * @return the engine talking to it
     * @throws IOException if the container does not become healthy in time
     */
    private static SageRestEngine startContainer(String image) throws IOException {
        int port = freePort();
        String name = "line-sage-rest-" + port;
        String id = run(120, "docker", "run", "-d", "--rm", "--name", name,
                "-p", port + ":8080", image);
        if (id == null || id.trim().length() == 0) {
            throw new IOException("could not start " + image);
        }
        startedContainer = name;
        Runtime.getRuntime().addShutdownHook(new Thread(new Runnable() {
            @Override
            public void run() {
                stopContainer();
            }
        }));

        SageRestEngine engine = new SageRestEngine("http://localhost:" + port);
        long deadline = System.currentTimeMillis() + STARTUP_TIMEOUT_SECONDS * 1000L;
        while (System.currentTimeMillis() < deadline) {
            if (engine.isAvailable()) {
                if (engine.isUsable()) {
                    started = engine;
                    return engine;
                }
                // Booted, but its arithmetic dies on this CPU. Keeping it
                // running would only cost memory, and returning it would hand
                // the caller a backend that kills every request.
                stopContainer();
                throw new IOException("container " + name + " answers but cannot evaluate on "
                        + "this CPU");
            }
            try {
                Thread.sleep(500);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                break;
            }
        }
        stopContainer();
        throw new IOException("container " + name + " did not become healthy within "
                + STARTUP_TIMEOUT_SECONDS + " s");
    }

    /** Stops the container started by this JVM, if any. */
    public static synchronized void stopContainer() {
        if (startedContainer != null) {
            run(30, "docker", "stop", "-t", "1", startedContainer);
            startedContainer = null;
            started = null;
        }
    }

    private static int freePort() throws IOException {
        ServerSocket socket = new ServerSocket(0);
        try {
            return socket.getLocalPort();
        } finally {
            socket.close();
        }
    }

    /**
     * Runs a command and returns its stdout, or null if it fails.
     *
     * @param timeoutSeconds how long to wait before giving up
     * @param command        the command and its arguments
     * @return the trimmed stdout, or null on failure
     */
    private static String run(int timeoutSeconds, String... command) {
        List<String> cmd = new ArrayList<String>();
        for (int i = 0; i < command.length; i++) {
            cmd.add(command[i]);
        }
        try {
            ProcessBuilder pb = new ProcessBuilder(cmd);
            pb.redirectErrorStream(false);
            final Process process = pb.start();
            StringBuilder sb = new StringBuilder();
            BufferedReader reader = new BufferedReader(
                    new InputStreamReader(process.getInputStream()));
            try {
                String line;
                while ((line = reader.readLine()) != null) {
                    sb.append(line).append("\n");
                }
            } finally {
                reader.close();
            }
            // see _kb/03-api-layer.md for rationale (sym/ section)
            int code = process.waitFor();
            return code == 0 ? sb.toString().trim() : null;
        } catch (IOException e) {
            return null;
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            return null;
        }
    }
}
