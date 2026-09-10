package jline.solvers.wrappers.jmt;

import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import jline.io.SysUtils;
import jline.solvers.SolverOptions;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;

/**
 * Backend dispatch for the JMT command line.
 *
 * One analysis of jmt.commandline.Jmt is run on a model file and the result is
 * left where the JMT CLI itself would leave it, that is at
 * &lt;model&gt;-result.jsim for mode "sim" and &lt;model&gt;-result.jmva for mode
 * "mva". Every backend satisfies that contract, so the result parsers do not
 * know or care which one ran.
 *
 * Backend selection, in order:
 * <ol>
 *   <li>options.restUrl non-empty: POST to a JMT REST server (the
 *       imperialqore/jmt-rest container). Nothing is executed locally.</li>
 *   <li>a local JVM plus common/JMT.jar: the default, unchanged.</li>
 *   <li>no local JVM, but Docker is usable: ask once per session whether to
 *       pull and use the JMT image, and dispatch through it.</li>
 * </ol>
 *
 * see _kb/06-solver-catalog.md (Wrappers: three ways to reach an external binary)
 */
public class JmtBackend {

    /** Docker images tried, in order, when none is named by the options. */
    public static final String[] DEFAULT_IMAGES = {
        "imperialqore/jmt-rest:latest", "imperialqore/jmt-rest"
    };

    // Session state: the JVM probe and the Docker consent are both asked once.
    private static Boolean hasJava = null;
    private static String dockerDecision = null;

    private JmtBackend() {
    }

    /**
     * Returns the path JMT writes its result to for this model and mode.
     *
     * @param modelPath path of the model file handed to JMT
     * @param mode      "sim" or "mva"
     * @return the result path
     */
    public static String resultPath(String modelPath, String mode) {
        return modelPath + resultSuffix(mode);
    }

    private static String resultSuffix(String mode) {
        if ("sim".equals(mode)) {
            return "-result.jsim";
        } else if ("mva".equals(mode)) {
            return "-result.jmva";
        }
        throw new IllegalArgumentException("Unknown JMT analysis mode: " + mode);
    }

    /**
     * Reports whether a JVM is on the path.
     *
     * Probed once per session: the check costs a process launch, and a JVM
     * does not appear or vanish mid-session.
     *
     * @return true if "java -version" runs and exits zero
     */
    public static synchronized boolean hasJava() {
        if (hasJava == null) {
            hasJava = Boolean.valueOf(runQuiet(new String[]{SysUtils.javaLauncher(), "-version"}) == 0);
        }
        return hasJava.booleanValue();
    }

    /**
     * Builds the shell command that runs one JMT analysis locally.
     *
     * @param jmtPath   path of JMT.jar
     * @param mode      "sim" or "mva"
     * @param modelPath path of the model file
     * @param seed      simulation seed, ignored by the JMVA engine
     * @return the command line, as passed to SysUtils.system
     */
    public static String localCommand(String jmtPath, String mode, String modelPath, int seed) {
        resultSuffix(mode);
        // Quoted: the launcher, the jar and the model all sit under paths with
        // spaces on a default Windows install, and SysUtils.system keeps a
        // quoted span whole.
        return String.format("\"%s\" -cp \"%s\" jmt.commandline.Jmt %s \"%s\" -seed %s --illegal-access=permit",
                SysUtils.javaLauncher(), jmtPath, mode, modelPath, seed);
    }

    /**
     * Runs one JMT analysis through the backend the options select, when that
     * backend is not the local JVM.
     *
     * @param mode      "sim" or "mva"
     * @param modelPath path of the model file
     * @param seed      simulation seed, ignored by the JMVA engine
     * @param options   solver options, read for restUrl and container
     * @return true if a non-local backend ran the analysis and wrote its
     *         result file; false if the caller should run the local command
     * @throws RuntimeException if no backend can run JMT, or a remote one failed
     */
    public static boolean runRemote(String mode, String modelPath, int seed, SolverOptions options) {
        String restUrl = (options != null) ? options.restUrl : null;
        if (restUrl != null && !restUrl.isEmpty()) {
            solveRest(restUrl, mode, modelPath, seed, options);
            return true;
        }

        if (hasJava()) {
            return false;
        }

        String image = resolveDockerImage(options);
        if (image == null || image.isEmpty()) {
            throw new RuntimeException("SolverJMT requires a Java runtime and JMT.jar. No JVM was "
                    + "found on the path, and Docker is not usable either. Install Java, or start a "
                    + "JMT REST server and set options.restUrl.");
        }
        runDocker(image, mode, modelPath, seed);
        return true;
    }

    /**
     * Resolves the JMT Docker image to dispatch through, with user consent.
     *
     * The decision is remembered for the session: a sweep over many models must
     * not ask once per model.
     *
     * @param options solver options, read for the container override
     * @return the image name, or null when Docker is unusable, no image can be
     *         obtained, or the user declines
     */
    public static synchronized String resolveDockerImage(SolverOptions options) {
        // Bind-mount dispatch is supported on unix hosts only, as for LQNS.
        if (System.getProperty("os.name").startsWith("Windows")) {
            return null;
        }
        if (runQuiet(new String[]{"docker", "info"}) != 0) {
            return null;
        }

        String requested = (options != null) ? options.container : null;
        if (requested == null || requested.isEmpty()) {
            requested = System.getenv("LINE_JMT_IMAGE");
        }
        String[] candidates = (requested != null && !requested.isEmpty())
                ? new String[]{requested} : DEFAULT_IMAGES;

        // An image already on the host needs no pull and no question.
        for (int i = 0; i < candidates.length; i++) {
            if (imagePresent(candidates[i])) {
                return candidates[i];
            }
        }

        if ("no".equals(dockerDecision)) {
            return null;
        }

        String target = candidates[0];
        if (!"yes".equals(dockerDecision) && !askDockerConsent(target)) {
            dockerDecision = "no";
            return null;
        }
        dockerDecision = "yes";

        // Storage check before pulling (shared guard used by the LQNS/QNS wrappers).
        if (!jline.io.DockerImage.hasStorageFor(target)) {
            dockerDecision = "no";
            System.err.println("[LINE] Skipping docker pull of " + target
                    + ": insufficient free space at the Docker storage location. Install Java, "
                    + "or set options.restUrl to a running JMT REST server.");
            return null;
        }

        System.out.println("Pulling " + target + " (this happens once)...");
        if (runQuiet(new String[]{"docker", "pull", target}) != 0) {
            dockerDecision = "no";
            throw new RuntimeException("Could not pull " + target + ". Install Java, or set "
                    + "options.restUrl to a running JMT REST server.");
        }
        return target;
    }

    private static boolean imagePresent(String image) {
        try {
            Process p = new ProcessBuilder("docker", "images", "-q", image)
                    .redirectErrorStream(true).start();
            String out = readAll(p.getInputStream());
            p.waitFor();
            return p.exitValue() == 0 && out.trim().length() > 0;
        } catch (IOException e) {
            return false;
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            return false;
        }
    }

    /**
     * Asks whether the JMT Docker image may be pulled and used.
     *
     * LINE_JMT_DOCKER answers for the user in unattended runs: "1" consents,
     * "0" refuses. Without it, a session with no console refuses rather than
     * blocking: a solver must never hang on a prompt nobody can answer, and
     * must never download without being asked.
     *
     * @param image the image that would be pulled
     * @return true if the pull may proceed
     */
    private static boolean askDockerConsent(String image) {
        String env = System.getenv("LINE_JMT_DOCKER");
        if (env != null && env.trim().length() > 0) {
            String v = env.trim().toLowerCase();
            return v.equals("1") || v.equals("true") || v.equals("yes") || v.equals("y");
        }

        if (System.console() == null) {
            return false;
        }
        System.out.println();
        System.out.println("SolverJMT needs Java, which was not found on this host.");
        System.out.println("Docker is available and can run JMT from the image " + image + " instead.");
        String answer = System.console().readLine("Pull and use that image? [y/N]: ");
        if (answer == null) {
            return false;
        }
        String v = answer.trim().toLowerCase();
        return v.equals("y") || v.equals("yes");
    }

    /**
     * Runs the analysis in the JMT container and copies the result back next to
     * the model, so the caller sees the layout the local JVM would have left.
     *
     * The model is staged under HOME rather than run in place: snap-confined
     * Docker cannot bind-mount the system temp dir where the JMT model normally
     * lives. In mva mode JMT rewrites the model file itself, which is why the
     * staged copy, not the original, is what the container touches.
     *
     * @param image     the Docker image to run
     * @param mode      "sim" or "mva"
     * @param modelPath path of the model file
     * @param seed      simulation seed
     */
    private static void runDocker(String image, String mode, String modelPath, int seed) {
        Path stageRoot = Paths.get(System.getProperty("user.home"), ".line", "line_workspace", "jmt-docker");
        Path workdir = null;
        try {
            Files.createDirectories(stageRoot);
            workdir = Files.createTempDirectory(stageRoot, "tmp_");
            String base = new File(modelPath).getName();
            Path staged = workdir.resolve(base);
            Files.copy(Paths.get(modelPath), staged, StandardCopyOption.REPLACE_EXISTING);

            String uid = readAll(new ProcessBuilder("id", "-u").start().getInputStream()).trim();
            String gid = readAll(new ProcessBuilder("id", "-g").start().getInputStream()).trim();

            ProcessBuilder pb = new ProcessBuilder("docker", "run", "--rm",
                    "--user", uid + ":" + gid,
                    "-v", workdir + ":" + workdir, "-w", workdir.toString(),
                    image, mode, base, "-seed", String.valueOf(seed));
            pb.redirectErrorStream(true);
            Process p = pb.start();
            String out = readAll(p.getInputStream());
            p.waitFor();

            Path stagedResult = Paths.get(resultPath(staged.toString(), mode));
            if (!Files.exists(stagedResult)) {
                throw new RuntimeException("JMT produced no result in the container: " + out.trim());
            }
            Files.copy(stagedResult, Paths.get(resultPath(modelPath, mode)),
                    StandardCopyOption.REPLACE_EXISTING);
        } catch (IOException e) {
            throw new RuntimeException("JMT Docker dispatch failed: " + e.getMessage(), e);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("JMT Docker dispatch was interrupted", e);
        } finally {
            deleteRecursively(workdir);
        }
    }

    /**
     * Solves through a JMT REST server and writes its result next to the model.
     *
     * The request carries the same JSIM or JMVA document the CLI reads, and the
     * response carries the same result document the CLI writes, so a fixed seed
     * gives the same numbers as the local backend.
     *
     * @param restUrl   base URL of the server, e.g. "http://localhost:8080"
     * @param mode      "sim" or "mva"
     * @param modelPath path of the model file
     * @param seed      simulation seed, sent only for the simulation route
     * @param options   solver options, read for the timeout
     */
    private static void solveRest(String restUrl, String mode, String modelPath, int seed,
                                  SolverOptions options) {
        String url = restUrl.replaceAll("/+$", "");
        if (!url.matches(".*/api/v\\d+/solve/(sim|mva)$")) {
            url = url + "/api/v1/solve/" + mode;
        }

        try {
            String modelText = new String(Files.readAllBytes(Paths.get(modelPath)),
                    StandardCharsets.UTF_8);

            JsonObject model = new JsonObject();
            model.addProperty("content", modelText);
            model.addProperty("base64", false);
            JsonObject request = new JsonObject();
            request.add("model", model);
            // JMVA takes its algorithm and tolerance from the model document, so
            // the seed is only meaningful for the simulation route; sending it to
            // /solve/mva would be rejected by the server's option allow-list.
            if ("sim".equals(mode)) {
                JsonObject opts = new JsonObject();
                opts.addProperty("seed", seed);
                request.add("options", opts);
            }

            JsonObject response = postJson(url, request.toString(), timeoutMillis(options));
            String status = response.has("status") ? response.get("status").getAsString() : "";
            if (!"completed".equals(status)) {
                String message = (response.has("error") && !response.get("error").isJsonNull())
                        ? response.get("error").getAsString() : "unspecified error";
                throw new RuntimeException("JMT REST solve failed: " + message.trim());
            }
            if (!response.has("raw_output") || response.get("raw_output").isJsonNull()) {
                throw new RuntimeException("JMT REST response carries no raw output. The server was "
                        + "asked to include it; check that include_raw_output is not disabled.");
            }
            JsonObject raw = response.getAsJsonObject("raw_output");
            if (!raw.has("result_xml") || raw.get("result_xml").isJsonNull()) {
                throw new RuntimeException("JMT REST response carries no result document.");
            }
            Files.write(Paths.get(resultPath(modelPath, mode)),
                    raw.get("result_xml").getAsString().getBytes(StandardCharsets.UTF_8));
        } catch (IOException e) {
            throw new RuntimeException("JMT REST request to " + url + " failed: " + e.getMessage(), e);
        }
    }

    private static int timeoutMillis(SolverOptions options) {
        double t = (options != null) ? options.timeout : Double.POSITIVE_INFINITY;
        if (Double.isFinite(t) && t > 0) {
            return (int) Math.min(Integer.MAX_VALUE, (t + 30.0) * 1000.0);
        }
        return 3600 * 1000;
    }

    private static JsonObject postJson(String url, String body, int timeoutMs) throws IOException {
        HttpURLConnection conn = (HttpURLConnection) new URL(url).openConnection();
        try {
            conn.setRequestMethod("POST");
            conn.setDoOutput(true);
            conn.setConnectTimeout(30000);
            conn.setReadTimeout(timeoutMs);
            conn.setRequestProperty("Content-Type", "application/json");
            byte[] payload = body.getBytes(StandardCharsets.UTF_8);
            OutputStream os = conn.getOutputStream();
            try {
                os.write(payload);
            } finally {
                os.close();
            }

            int code = conn.getResponseCode();
            InputStream stream = (code >= 400) ? conn.getErrorStream() : conn.getInputStream();
            String text = (stream != null) ? readAll(stream) : "";
            // The service answers 200 for both "completed" and "failed", and 500
            // only when it could not run JMT at all; both bodies are JSON.
            if (text.isEmpty()) {
                throw new IOException("empty response, HTTP " + code);
            }
            return JsonParser.parseString(text).getAsJsonObject();
        } finally {
            conn.disconnect();
        }
    }

    private static int runQuiet(String[] command) {
        try {
            ProcessBuilder pb = new ProcessBuilder(command);
            pb.redirectErrorStream(true);
            Process p = pb.start();
            readAll(p.getInputStream());
            p.waitFor();
            return p.exitValue();
        } catch (IOException e) {
            return -1;
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            return -1;
        }
    }

    private static String readAll(InputStream stream) throws IOException {
        BufferedReader reader = new BufferedReader(
                new InputStreamReader(stream, StandardCharsets.UTF_8));
        StringBuilder sb = new StringBuilder();
        String line;
        try {
            while ((line = reader.readLine()) != null) {
                sb.append(line).append('\n');
            }
        } finally {
            reader.close();
        }
        return sb.toString();
    }

    private static void deleteRecursively(Path root) {
        if (root == null) {
            return;
        }
        try {
            File[] files = root.toFile().listFiles();
            if (files != null) {
                for (int i = 0; i < files.length; i++) {
                    files[i].delete();
                }
            }
            root.toFile().delete();
        } catch (SecurityException e) {
            // A leftover staging directory is harmless; failing the solve is not.
        }
    }
}
