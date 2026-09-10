/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Docker primitives shared by the backends that legitimately ship an image: the
 * JMT backend ({@code imperialqore/jmt-rest}) and the Sage symbolic engine
 * ({@code imperialqore/line-sage-rest}). Each owns its image name; this class
 * only answers whether the daemon is up, whether an image is present, whether
 * there is room to pull it, and performs the pull.
 *
 * <p>LQNS, lqsim and qnsolver are deliberately absent: their licence is an
 * evaluation agreement that forbids redistribution, so LINE runs them only from
 * a binary the user installed. {@code run-tests.sh --lqns-docker} puts a shim on
 * the PATH when a containerised build is what should be exercised.
 *
 * <p>Every method is a no-op / returns empty on Windows.
 */
public final class DockerImage {
    private DockerImage() {}

    /** Conservative free-space floor required before a pull (2 GiB). */
    private static final long DEFAULT_MIN_FREE_BYTES = 2L * 1024 * 1024 * 1024;

    private static boolean isWindows() {
        return System.getProperty("os.name", "").toLowerCase().contains("win");
    }

    /** True if the Docker daemon is reachable (unix only). */
    public static boolean daemonAvailable() {
        if (isWindows()) return false;
        return runQuiet(new String[]{"docker", "info"}, 15) == 0;
    }

    /** True if the named image is already present in the local Docker store. */
    public static boolean hasLocalImage(String image) {
        if (isWindows() || image == null || image.isEmpty()) return false;
        String out = capture(new String[]{"docker", "images", "-q", image}, 15);
        return out != null && !out.trim().isEmpty();
    }

    /** Pull an image, streaming Docker's progress to stdout/stderr. No timeout. */
    public static boolean pull(String image) {
        if (isWindows() || image == null || image.isEmpty()) return false;
        try {
            Process p = new ProcessBuilder("docker", "pull", image).inheritIO().start();
            return p.waitFor() == 0;
        } catch (Exception e) {
            return false;
        }
    }

    /** True if the Docker storage location has room for {@code image}. */
    public static boolean hasStorageFor(String image) {
        long free = freeBytesAtDockerRoot();
        if (free < 0) return true; // could not determine; do not block the pull
        return free >= requiredBytes(image);
    }

    private static long requiredBytes(String image) {
        long override = overrideMinFreeBytes();
        if (override > 0) return override;
        long est = estimateImageBytes(image);
        return Math.max(DEFAULT_MIN_FREE_BYTES, est);
    }

    private static long overrideMinFreeBytes() {
        String v = System.getProperty("line.docker.minFreeBytes");
        if (v == null || v.isEmpty()) v = System.getenv("LINE_DOCKER_MIN_FREE_BYTES");
        if (v != null && !v.isEmpty()) {
            try { return Long.parseLong(v.trim()); } catch (NumberFormatException e) { /* ignore */ }
        }
        return -1;
    }

    /** Usable bytes on the filesystem backing the Docker root dir; -1 if unknown. */
    static long freeBytesAtDockerRoot() {
        String root = capture(new String[]{"docker", "info", "--format", "{{.DockerRootDir}}"}, 15);
        if (root != null) root = root.trim();
        if (root == null || root.isEmpty()) root = "/var/lib/docker";
        File f = new File(root);
        // The root dir may not exist for / be readable by this user: walk up to an
        // existing ancestor so getUsableSpace() reports a real filesystem.
        while (f != null && !f.exists()) f = f.getParentFile();
        if (f == null) f = new File("/");
        try {
            long usable = f.getUsableSpace();
            return usable > 0 ? usable : -1;
        } catch (Exception e) {
            return -1;
        }
    }

    /** On-disk estimate (compressed layer sizes x3), or 0 if it cannot be determined. */
    static long estimateImageBytes(String image) {
        String json = capture(new String[]{"docker", "manifest", "inspect", image}, 30);
        if (json == null || json.isEmpty()) return 0L;
        long sum = 0L;
        Matcher m = Pattern.compile("\"size\"\\s*:\\s*(\\d+)").matcher(json);
        while (m.find()) {
            try { sum += Long.parseLong(m.group(1)); } catch (NumberFormatException e) { /* skip */ }
        }
        return sum > 0 ? sum * 3L : 0L;
    }

    // --- process helpers -----------------------------------------------------

    /** Run discarding all output; return the exit code (124 on timeout, 127 on error). */
    private static int runQuiet(String[] cmd, int timeoutSeconds) {
        try {
            ProcessBuilder pb = new ProcessBuilder(cmd);
            pb.redirectOutput(new File("/dev/null"));
            pb.redirectError(new File("/dev/null"));
            Process p = pb.start();
            if (timeoutSeconds > 0) {
                if (!p.waitFor(timeoutSeconds, TimeUnit.SECONDS)) { p.destroyForcibly(); return 124; }
                return p.exitValue();
            }
            return p.waitFor();
        } catch (Exception e) {
            return 127;
        }
    }

    /** Run capturing stdout; return it on exit 0, else null. stderr is discarded. */
    private static String capture(String[] cmd, int timeoutSeconds) {
        try {
            ProcessBuilder pb = new ProcessBuilder(cmd);
            pb.redirectError(new File("/dev/null"));
            Process p = pb.start();
            StringBuilder sb = new StringBuilder();
            BufferedReader r = new BufferedReader(new InputStreamReader(p.getInputStream()));
            try {
                String line;
                while ((line = r.readLine()) != null) sb.append(line).append('\n');
            } finally {
                r.close();
            }
            if (timeoutSeconds > 0) {
                if (!p.waitFor(timeoutSeconds, TimeUnit.SECONDS)) { p.destroyForcibly(); return null; }
            } else {
                p.waitFor();
            }
            return p.exitValue() == 0 ? sb.toString() : null;
        } catch (Exception e) {
            return null;
        }
    }
}
