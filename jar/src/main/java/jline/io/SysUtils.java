/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.io;

import jline.GlobalConstants;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URI;
import java.nio.file.FileVisitResult;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.SimpleFileVisitor;
import java.nio.file.StandardCopyOption;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.concurrent.TimeUnit;

public final class SysUtils {
    private SysUtils() {}

    public static String getRootFolder() {
        try {
            URI uri = GlobalConstants.class.getProtectionDomain().getCodeSource().getLocation().toURI();
            return Paths.get(uri).getParent().toString();
        } catch (Exception e) {
            try {
                ClassLoader classLoader = Thread.currentThread().getContextClassLoader();
                java.net.URL resource = classLoader.getResource("jline/Scratch.class");
                if (resource != null) {
                    String path = resource.getPath();
                    if (path.contains(".jar!")) {
                        String jarPath = path.substring(0, path.indexOf(".jar!") + 4);
                        if (jarPath.startsWith("file:")) {
                            return Paths.get(jarPath.substring(5)).getParent().toString();
                        } else {
                            return Paths.get(jarPath).getParent().toString();
                        }
                    } else {
                        return Paths.get(path).getParent().getParent().getParent().toString();
                    }
                }
            } catch (Exception e2) {
                // fall through
            }
            try {
                String userDir = System.getProperty("user.dir");
                if (userDir != null) return userDir;
            } catch (Exception e3) {
                // fall through
            }
            try {
                ClassLoader loader = Thread.currentThread().getContextClassLoader();
                if (loader == null) loader = ClassLoader.getSystemClassLoader();
                java.net.URL resource = loader.getResource("");
                if (resource != null) return Paths.get(resource.toURI()).toString();
            } catch (Exception e4) {
                // fall through
            }
            String userDir = System.getProperty("user.dir");
            return (userDir != null) ? userDir : System.getProperty("java.io.tmpdir");
        }
    }

    public static String jlqnGetPath() {
        String rf = getRootFolder();
        String path;
        if (rf != null) {
            Path commonPath = Paths.get(rf, "common", "JLQN.jar");
            if (Files.exists(commonPath)) {
                path = commonPath.toString();
            } else {
                path = Paths.get(rf, "JLQN.jar").toString();
            }
        } else {
            String currentDir = System.getProperty("user.dir");
            if (currentDir == null) currentDir = System.getProperty("java.io.tmpdir");
            if (currentDir == null) currentDir = ".";
            Path commonPath = Paths.get(currentDir, "common", "JLQN.jar");
            if (Files.exists(commonPath)) {
                path = commonPath.toString();
            } else {
                path = Paths.get(currentDir, "JLQN.jar").toString();
            }
        }
        return jlqnGetPath(path);
    }

    public static String jlqnGetPath(String jlqnPath) {
        File jlqnFile = new File(jlqnPath);
        if (!jlqnFile.exists()) {
            System.out.println("\nJLQN GUI cannot be found. LINE will try to download the latest JLQN version (download approx. 50MB).");
            String m = "Y";
            if (m.equalsIgnoreCase("Y")) {
                System.out.println("Download started, please wait - this may take several minutes.");
                try {
                    java.net.URL url = new URI("https://github.com/imperial-qore/JLQN/raw/main/target/jlqn-singlejar.jar").toURL();
                    Files.copy(url.openStream(), jlqnFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
                    System.out.println("Download completed. JLQN jar now located at: " + jlqnPath);
                } catch (Exception e) {
                    jlqnFile.delete();
                    e.printStackTrace();
                }
            } else {
                System.out.println("JLQN was not found. Please download it manually and place it in the root folder.");
            }
        }
        return jlqnPath;
    }

    /**
     * Directory holding the running jline.jar, or null when it cannot be
     * resolved (e.g. an exploded classpath with no code source).
     *
     * <p>This is the JAR's analogue of the anchor the other two codebases use:
     * MATLAB {@code jmtGetPath.m} walks up from {@code mfilename('fullpath')}
     * and native Python from {@code __file__}. Neither consults the working
     * directory, because the caller's cwd says nothing about where the
     * installation is.</p>
     *
     * @return the directory containing the running jar, or null
     */
    private static Path runningJarDir() {
        try {
            URI jarLocation = GlobalConstants.class.getProtectionDomain().getCodeSource().getLocation().toURI();
            Path jarPath = Paths.get(jarLocation);
            return Files.isDirectory(jarPath) ? jarPath : jarPath.getParent();
        } catch (Exception e) {
            return null;
        }
    }

    /** The JMT.jar under {@code dir} or {@code dir/common}, if one is really there. */
    private static Path jmtJarUnder(Path dir) {
        if (dir == null) {
            return null;
        }
        Path here = dir.resolve("JMT.jar");
        if (Files.isRegularFile(here)) {
            return here;
        }
        Path nested = dir.resolve("common").resolve("JMT.jar");
        return Files.isRegularFile(nested) ? nested : null;
    }

    /**
     * Resolves the Java launcher used to spawn a JVM tool (JMT, LDES).
     *
     * <p>Resolution order: {@code $LINE_JAVA}; {@code $JAVA_HOME/bin/java};
     * the JVM running this code ({@code java.home}); then the bare name
     * {@code java}, left to PATH. The {@code java.home} step matters where
     * PATH carries no java at all, which is the common case on Windows: a
     * process already inside a JVM always has one launcher it can name.</p>
     *
     * @return an absolute launcher path, or "java" when only PATH can answer
     */
    public static String javaLauncher() {
        String exeName = System.getProperty("os.name", "").startsWith("Windows") ? "java.exe" : "java";
        String override = System.getenv("LINE_JAVA");
        if (override != null && !override.trim().isEmpty() && new File(override.trim()).isFile()) {
            return override.trim();
        }
        String javaHome = System.getenv("JAVA_HOME");
        if (javaHome != null && !javaHome.trim().isEmpty()) {
            File cand = new File(new File(javaHome.trim(), "bin"), exeName);
            if (cand.isFile()) {
                return cand.getAbsolutePath();
            }
        }
        String runningHome = System.getProperty("java.home");
        if (runningHome != null && !runningHome.isEmpty()) {
            File cand = new File(new File(runningHome, "bin"), exeName);
            if (cand.isFile()) {
                return cand.getAbsolutePath();
            }
        }
        return "java";
    }

    /**
     * Locates {@code JMT.jar}.
     *
     * <p>Resolution order: {@code $LINE_JMT_JAR}; the directory of the running
     * jar and its ancestors; then the working directory and its ancestors. A
     * candidate counts only when {@code JMT.jar} is actually in it. The walk
     * used to accept the first ancestor merely NAMED {@code common}, so an
     * empty {@code python/common} shadowed the real {@code common/} for every
     * invocation whose cwd was under {@code python/} and the CLI went off to
     * download 50MB instead. See {@code _kb/12-interfaces-and-docs.md}.</p>
     *
     * @return the path of JMT.jar, existing where one was found
     */
    public static String jmtGetPath() {
        String override = System.getenv("LINE_JMT_JAR");
        if (override != null && !override.trim().isEmpty()) {
            return jmtGetPath(override.trim());
        }

        for (Path dir = runningJarDir(); dir != null; dir = dir.getParent()) {
            Path found = jmtJarUnder(dir);
            if (found != null) {
                return found.toString();
            }
        }

        String currentDir = System.getProperty("user.dir");
        if (currentDir == null) currentDir = System.getProperty("java.io.tmpdir");
        if (currentDir == null) currentDir = ".";
        for (Path dir = Paths.get(currentDir).toAbsolutePath(); dir != null; dir = dir.getParent()) {
            Path found = jmtJarUnder(dir);
            if (found != null) {
                return found.toString();
            }
        }

        // Nothing found: name the canonical location without creating it. Appending
        // "common" to a directory already called that is what made common/common.
        Path jarDir = runningJarDir();
        Path base = (jarDir != null) ? jarDir : Paths.get(currentDir).toAbsolutePath();
        Path target = "common".equals(String.valueOf(base.getFileName()))
                ? base.resolve("JMT.jar")
                : base.resolve("common").resolve("JMT.jar");
        return jmtGetPath(target.toString());
    }

    public static String jmtGetPath(String jmtPath) {
        File jmtFile = new File(jmtPath);
        if (!jmtFile.exists()) {
            // LINE_JMT_DOWNLOAD=0 refuses the fetch, for unattended runs
            String consent = System.getenv("LINE_JMT_DOWNLOAD");
            if (consent != null && (consent.trim().equals("0") || consent.trim().equalsIgnoreCase("false"))) {
                System.out.println("\nJava Modelling Tools cannot be found at " + jmtPath
                        + " and LINE_JMT_DOWNLOAD forbids downloading it. Place JMT.jar there,"
                        + " or point LINE_JMT_JAR at an existing copy.");
                return jmtPath;
            }
            System.out.println("\nJava Modelling Tools cannot be found. LINE will try to download the latest JMT version (download approx. 50MB).");
            String m = "Y";
            if (m.equalsIgnoreCase("Y")) {
                System.out.println("Download started, please wait - this may take several minutes.");
                try {
                    java.net.URL url = new URI("https://line-solver.sourceforge.net/latest/JMT.jar").toURL();
                    Files.copy(url.openStream(), jmtFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
                    System.out.println("Download completed. JMT.jar now located at: " + jmtPath);
                } catch (Exception e) {
                    System.out.println("Download failed. JMT.jar could not be saved at: " + jmtPath);
                    jmtFile.delete();
                    e.printStackTrace();
                }
            } else {
                System.out.println("JMT was not found. Please download it manually and place it in the root folder.");
            }
        }
        return jmtPath;
    }

    public static String lineViewerGetPath() {
        String currentDir = System.getProperty("user.dir");
        if (currentDir == null) currentDir = System.getProperty("java.io.tmpdir");
        if (currentDir == null) currentDir = ".";
        Path currentPath = Paths.get(currentDir);
        while (currentPath != null && currentPath.getNameCount() > 0) {
            Path commonPath = currentPath.resolve("common");
            if (Files.exists(commonPath) && Files.isDirectory(commonPath)) {
                Path jarPath = commonPath.resolve("line-viewer.jar");
                if (Files.exists(jarPath)) return jarPath.toString();
            }
            currentPath = currentPath.getParent();
        }

        String rf = getRootFolder();
        try {
            if (rf != null) {
                Path commonPath = Paths.get(rf, "common");
                if (!Files.exists(commonPath)) Files.createDirectories(commonPath);
                return commonPath.resolve("line-viewer.jar").toString();
            } else {
                URI jarLocation = GlobalConstants.class.getProtectionDomain().getCodeSource().getLocation().toURI();
                Path jarPath = Paths.get(jarLocation);
                Path jarDir = Files.isDirectory(jarPath) ? jarPath : jarPath.getParent();
                Path commonPath = jarDir.resolve("common");
                if (!Files.exists(commonPath)) Files.createDirectories(commonPath);
                return commonPath.resolve("line-viewer.jar").toString();
            }
        } catch (Exception e) {
            try {
                Path commonPath = Paths.get(currentDir, "common");
                if (!Files.exists(commonPath)) Files.createDirectories(commonPath);
                return commonPath.resolve("line-viewer.jar").toString();
            } catch (IOException ioe) {
                return Paths.get(currentDir, "common", "line-viewer.jar").toString();
            }
        }
    }

    public static String system(String cmd) {
        return system(cmd, 0L);
    }

    /**
     * Splits a command line into an argument vector, keeping a double-quoted
     * span together.
     *
     * <p>Splitting on whitespace alone tore apart every path with a space in
     * it, which on Windows is where both the JVM ({@code C:\Program
     * Files\...\java.exe}) and the model normally live.</p>
     *
     * @param cmd the command line
     * @return its arguments, with the quotes removed
     */
    public static java.util.List<String> splitCommand(String cmd) {
        java.util.List<String> argv = new java.util.ArrayList<String>();
        StringBuilder current = new StringBuilder();
        boolean quoted = false;
        boolean started = false;
        for (int i = 0; i < cmd.length(); i++) {
            char c = cmd.charAt(i);
            if (c == '"') {
                quoted = !quoted;
                started = true;
            } else if (!quoted && Character.isWhitespace(c)) {
                if (started) {
                    argv.add(current.toString());
                    current.setLength(0);
                    started = false;
                }
            } else {
                current.append(c);
                started = true;
            }
        }
        if (started) {
            argv.add(current.toString());
        }
        return argv;
    }

    public static String system(String cmd, long timeoutSeconds) {
        final StringBuilder output = new StringBuilder();
        final StringBuilder errorOutput = new StringBuilder();
        final Process process;
        try {
            process = new ProcessBuilder(splitCommand(cmd)).start();
        } catch (IOException e) {
            e.printStackTrace();
            return e.getMessage() == null ? "Failed to execute command" : e.getMessage();
        }

        Thread stdoutThread = new Thread(new Runnable() {
            public void run() {
                try {
                    BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                    try {
                        String line;
                        while ((line = reader.readLine()) != null) {
                            output.append(line).append('\n');
                        }
                    } finally {
                        reader.close();
                    }
                } catch (Exception e) {
                    // ignore
                }
            }
        });
        Thread stderrThread = new Thread(new Runnable() {
            public void run() {
                try {
                    BufferedReader reader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
                    try {
                        String line;
                        while ((line = reader.readLine()) != null) {
                            errorOutput.append(line).append('\n');
                        }
                    } finally {
                        reader.close();
                    }
                } catch (Exception e) {
                    // ignore
                }
            }
        });

        stdoutThread.start();
        stderrThread.start();

        boolean completed;
        try {
            if (timeoutSeconds > 0) {
                completed = process.waitFor(timeoutSeconds, TimeUnit.SECONDS);
            } else {
                process.waitFor();
                completed = true;
            }
        } catch (InterruptedException ie) {
            Thread.currentThread().interrupt();
            return "Interrupted";
        }

        if (!completed) {
            process.destroyForcibly();
            try { process.getInputStream().close(); } catch (Exception e) { }
            try { process.getErrorStream().close(); } catch (Exception e) { }
            try { process.waitFor(5, TimeUnit.SECONDS); } catch (Exception e) { }
            try { stdoutThread.join(1000); } catch (Exception e) { }
            try { stderrThread.join(1000); } catch (Exception e) { }
            return "TIMEOUT: Command exceeded " + timeoutSeconds + "s limit";
        }

        try {
            stdoutThread.join();
            stderrThread.join();
        } catch (InterruptedException ie) {
            Thread.currentThread().interrupt();
        }

        int exitCode = process.exitValue();
        return (exitCode == 0) ? output.toString() : errorOutput.toString();
    }

    public static String lineTempName() throws IOException {
        return lineTempName("");
    }

    public static String lineTempName(String solverName) throws IOException {
        return lineTempName(solverName, false);
    }

    /**
     * As {@link #lineTempName(String)}, but when {@code mountable} is true the
     * workspace is rooted under {@code $HOME/.line} instead of the system temp
     * dir. Snap-confined Docker cannot bind-mount the system temp dir (typically
     * {@code /tmp}), so solvers that dispatch through Docker request a HOME-based,
     * mount-accessible location (mirrors the MATLAB {@code lineTempName}).
     */
    public static String lineTempName(String solverName, boolean mountable) throws IOException {
        Path baseRoot;
        // LINE_WORKSPACE_ROOT relocates every staged model. run-tests.sh sets it when
        // wrapping a solver in a container, so the staging dir is one the container
        // can bind-mount; no solver needs to know a container is involved.
        String envRoot = System.getenv("LINE_WORKSPACE_ROOT");
        if (envRoot != null && !envRoot.trim().isEmpty()) {
            baseRoot = Paths.get(envRoot.trim());
        } else if (mountable) {
            String home = System.getProperty("user.home");
            baseRoot = (home != null && !home.isEmpty())
                    ? Paths.get(home, ".line")
                    : Files.createTempFile("", ".tmp").toAbsolutePath().getParent();
        } else {
            baseRoot = Files.createTempFile("", ".tmp").toAbsolutePath().getParent();
        }
        Path workspacePath = baseRoot.resolve("workspace").resolve(solverName == null ? "" : solverName);
        if (!Files.exists(workspacePath)) Files.createDirectories(workspacePath);
        return Files.createTempDirectory(workspacePath, "").toString();
    }

    public static void removeDirectory(final Path dirPath) throws IOException {
        Files.walkFileTree(dirPath, new SimpleFileVisitor<Path>() {
            @Override
            public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
                Files.delete(file);
                return FileVisitResult.CONTINUE;
            }

            @Override
            public FileVisitResult postVisitDirectory(Path dir, IOException exc) throws IOException {
                Files.delete(dir);
                return FileVisitResult.CONTINUE;
            }
        });
    }
}
