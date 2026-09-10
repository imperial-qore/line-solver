/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.cli;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import jline.api.sym.SymEngines;

/**
 * Environment check for the JAR, mirroring MATLAB's {@code lineInstall}.
 *
 * <p>It verifies that the optional external dependencies used by some solvers
 * are reachable and warns, without failing, when one is missing. In particular
 * it reports whether the SageMath symbolic backend (the
 * {@code imperialqore/line-sage-rest} Docker image, the JAR's only computer
 * algebra system, required by the symbolic methods of SolverCTMC/SolverFluid)
 * can be started on first use.</p>
 *
 * <p>Run it with {@code java -cp jline.jar jline.cli.LineInstall}.</p>
 */
public final class LineInstall {

    private LineInstall() {
    }

    /**
     * Runs the checks and prints warnings to stderr.
     *
     * @return true when everything is in place, false when a warning was issued
     */
    public static boolean check() {
        boolean hasWarnings = false;

        System.out.println("Checking Java...");
        String javaVersion = System.getProperty("java.version");
        System.out.println("  Java " + (javaVersion == null ? "unknown" : javaVersion));

        System.out.println("Checking symbolic backend (line-sage-rest)...");
        if (!commandOk(30, "docker", "info")) {
            System.err.println("WARNING: Docker is not available, so the SageMath symbolic "
                    + "backend cannot start. It is the JAR's only computer algebra system and "
                    + "is required by the symbolic methods of SolverCTMC/SolverFluid. Install "
                    + "Docker, then run: docker run -d -p 8080:8080 " + SymEngines.DOCKER_IMAGE);
            hasWarnings = true;
        } else if (SymEngines.findImage() == null) {
            System.err.println("WARNING: the line-sage-rest image is not present locally, this "
                    + "may be required by some LINE methods. Pull it with: "
                    + "docker run -d -p 8080:8080 " + SymEngines.DOCKER_IMAGE);
            hasWarnings = true;
        }

        if (hasWarnings) {
            System.out.println("Completed. LINE has warnings.");
        } else {
            System.out.println("Success. LINE is ready to use.");
        }
        return !hasWarnings;
    }

    /**
     * Runs a command and reports whether it exits zero.
     *
     * @param timeoutSeconds how long to wait before giving up
     * @param command        the command and its arguments
     * @return true if the command ran and exited zero, false on any failure
     */
    private static boolean commandOk(int timeoutSeconds, String... command) {
        List<String> cmd = new ArrayList<String>();
        for (int i = 0; i < command.length; i++) {
            cmd.add(command[i]);
        }
        try {
            ProcessBuilder pb = new ProcessBuilder(cmd);
            pb.redirectErrorStream(true);
            Process process = pb.start();
            BufferedReader reader = new BufferedReader(
                    new InputStreamReader(process.getInputStream()));
            try {
                while (reader.readLine() != null) {
                    // drain output so the process does not block on a full pipe
                }
            } finally {
                reader.close();
            }
            if (!process.waitFor(timeoutSeconds, java.util.concurrent.TimeUnit.SECONDS)) {
                process.destroy();
                return false;
            }
            return process.exitValue() == 0;
        } catch (Exception e) {
            return false;
        }
    }

    public static void main(String[] args) {
        check();
    }
}
