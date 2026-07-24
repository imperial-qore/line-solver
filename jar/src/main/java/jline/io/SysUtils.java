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

    public static String jmtGetPath() {
        String currentDir = System.getProperty("user.dir");
        if (currentDir == null) currentDir = System.getProperty("java.io.tmpdir");
        if (currentDir == null) currentDir = ".";
        Path currentPath = Paths.get(currentDir);
        while (currentPath != null && currentPath.getNameCount() > 0) {
            Path commonPath = currentPath.resolve("common");
            if (Files.exists(commonPath) && Files.isDirectory(commonPath)) {
                return jmtGetPath(commonPath.resolve("JMT.jar").toString());
            }
            currentPath = currentPath.getParent();
        }

        String rf = getRootFolder();
        String path;
        try {
            if (rf != null) {
                Path commonPath = Paths.get(rf, "common");
                if (!Files.exists(commonPath)) Files.createDirectories(commonPath);
                path = commonPath.resolve("JMT.jar").toString();
            } else {
                URI jarLocation = GlobalConstants.class.getProtectionDomain().getCodeSource().getLocation().toURI();
                Path jarPath = Paths.get(jarLocation);
                Path jarDir = Files.isDirectory(jarPath) ? jarPath : jarPath.getParent();
                Path commonPath = jarDir.resolve("common");
                if (!Files.exists(commonPath)) Files.createDirectories(commonPath);
                path = commonPath.resolve("JMT.jar").toString();
            }
        } catch (Exception e) {
            try {
                Path commonPath = Paths.get(currentDir, "common");
                if (!Files.exists(commonPath)) Files.createDirectories(commonPath);
                path = commonPath.resolve("JMT.jar").toString();
            } catch (IOException ioe) {
                path = Paths.get(currentDir, "common", "JMT.jar").toString();
            }
        }
        return jmtGetPath(path);
    }

    public static String jmtGetPath(String jmtPath) {
        File jmtFile = new File(jmtPath);
        if (!jmtFile.exists()) {
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

    public static String system(String cmd, long timeoutSeconds) {
        final StringBuilder output = new StringBuilder();
        final StringBuilder errorOutput = new StringBuilder();
        final Process process;
        try {
            process = new ProcessBuilder(cmd.split("\\s+")).start();
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
        Path tempDir = Files.createTempFile("", ".tmp").toAbsolutePath().getParent();
        Path workspacePath = tempDir.resolve("workspace").resolve(solverName == null ? "" : solverName);
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
