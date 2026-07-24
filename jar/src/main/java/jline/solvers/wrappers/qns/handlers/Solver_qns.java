package jline.solvers.wrappers.qns.handlers;

import java.io.File;
import java.util.HashMap;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.api.sn.SnDeaggregateChainResults;
import jline.api.sn.SnGetDemandsChain;
import jline.io.SysUtils;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.wrappers.jmt.SolverJMT;
import jline.solvers.wrappers.qns.QNSResult;
import jline.util.matrix.Matrix;

/**
 * Core handler for the QNS solver.
 */
public class Solver_qns {
    private final NetworkStruct sn;
    private final SolverOptions options;
    private String actualMethod;

    public Solver_qns(NetworkStruct sn, SolverOptions options) {
        this.sn = sn;
        this.options = options;
        this.actualMethod = (options.method != null) ? options.method : "default";
    }

    public QNSResult solve() {
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix QN = Matrix.zeros(M, K);
        Matrix UN = Matrix.zeros(M, K);
        Matrix RN = Matrix.zeros(M, K);
        Matrix TN = Matrix.zeros(M, K);
        Matrix AN = Matrix.zeros(M, K);
        Matrix WN = Matrix.zeros(M, K);
        Matrix CN = Matrix.zeros(1, K);
        Matrix XN = Matrix.zeros(1, K);

        actualMethod = (options.method != null) ? options.method : "default";

        if ("conway".equals(options.method)) options.config.multiserver = "conway";
        else if ("reiser".equals(options.method)) options.config.multiserver = "reiser";
        else if ("rolia".equals(options.method)) options.config.multiserver = "rolia";
        else if ("zhou".equals(options.method)) options.config.multiserver = "zhou";
        else if ("suri".equals(options.method)) options.config.multiserver = "suri";
        else if ("schmidt".equals(options.method)) options.config.multiserver = "schmidt";

        String tempDirPath;
        try {
            tempDirPath = SysUtils.lineTempName("qns");
        } catch (java.io.IOException ex) {
            throw new RuntimeException("QNS: Cannot allocate temp directory", ex);
        }
        File tempDir = new File(tempDirPath);

        try {
            File modelFile = new File(tempDir, "model.jmva");
            writeJMVAFile(modelFile);
            File resultFile = new File(tempDir, "result.jmva");
            File logFile = new File(tempDir, "console.out");
            String cmd = buildCommand(modelFile, resultFile, logFile);
            if (GlobalConstants.Verbose == VerboseLevel.DEBUG) {
                System.out.println("SolverQNS command: " + cmd);
            }
            File devNull = new File(isWindows() ? "NUL" : "/dev/null");
            int exitCode;
            if (isWindows()) {
                Process process = new ProcessBuilder("cmd", "/c", cmd)
                        .directory(tempDir)
                        .redirectOutput(devNull)
                        .redirectError(devNull)
                        .start();
                exitCode = process.waitFor();
            } else {
                Process process = new ProcessBuilder("sh", "-c", cmd)
                        .directory(tempDir)
                        .redirectOutput(devNull)
                        .redirectError(devNull)
                        .start();
                exitCode = process.waitFor();
            }
            if (exitCode != 0) {
                String logContent = logFile.exists() ? readFile(logFile) : "No log file";
                throw new RuntimeException("QNS solver failed with exit code: " + exitCode + "\nLog: " + logContent);
            }

            ParsedResults parsed = parseResults(resultFile, sn.nchains);
            Ret.snGetDemands demandResults = SnGetDemandsChain.snGetDemandsChain(sn);
            Matrix Lchain = demandResults.Dchain;
            Matrix STchain = demandResults.STchain;
            Matrix Vchain = demandResults.Vchain;
            Matrix alpha = demandResults.alpha;
            Matrix Xchain = Matrix.zeros(1, sn.nchains);
            for (int c = 0; c < sn.nchains; c++) {
                int refstat = (int) sn.refstat.get(c);
                double tChainRefstat = parsed.Tchain.get(refstat, c);
                if (tChainRefstat > 0) Xchain.set(0, c, tChainRefstat);
                else {
                    for (int i = 0; i < sn.nstations; i++) {
                        if (Vchain.get(i, c) > 0 && parsed.Tchain.get(i, c) > 0) {
                            Xchain.set(0, c, parsed.Tchain.get(i, c) / Vchain.get(i, c));
                            break;
                        }
                    }
                }
            }
            Matrix Rchain = parsed.Wchain.copy();
            Rchain.removeNaN();
            for (int i = 0; i < sn.nstations; i++) {
                double servers = sn.nservers.get(i);
                if (!Double.isInfinite(servers) && servers > 1) {
                    for (int c = 0; c < sn.nchains; c++) {
                        parsed.Uchain.set(i, c, parsed.Uchain.get(i, c) / servers);
                    }
                }
            }
            Ret.snDeaggregateChainResults results = SnDeaggregateChainResults.snDeaggregateChainResults(
                    sn, Lchain, null, STchain, Vchain, alpha, parsed.Qchain, parsed.Uchain, Rchain, parsed.Tchain, null, Xchain);
            QN.setTo(results.Q);
            UN.setTo(results.U);
            RN.setTo(results.R);
            TN.setTo(results.T);
            CN.setTo(results.C);
            XN.setTo(results.X);
            WN.setTo(results.R);
        } catch (Exception ex) {
            if (ex instanceof RuntimeException) throw (RuntimeException) ex;
            throw new RuntimeException(ex);
        } finally {
            deleteDirectory(tempDir);
        }
        return new QNSResult(QN, UN, RN, TN, AN, WN, CN, XN, 0.0, actualMethod, 0);
    }

    private static String readFile(File f) {
        try {
            return new String(java.nio.file.Files.readAllBytes(f.toPath()));
        } catch (Exception e) {
            return "";
        }
    }

    private void writeJMVAFile(File modelFile) {
        SolverJMT.writeJMVA(sn, modelFile.getAbsolutePath(), options);
    }

    private String buildCommand(File modelFile, File resultFile, File logFile) {
        StringBuilder cmd = new StringBuilder();
        cmd.append("qnsolver");
        cmd.append(" -l ").append(modelFile.getAbsolutePath());
        if (hasMultiServer() && options.config.multiserver != null) {
            String ms = options.config.multiserver;
            if ("default".equals(ms) || "conway".equals(ms)) {
                cmd.append(" -mconway");
                actualMethod = "conway";
            } else if ("reiser".equals(ms)) {
                cmd.append(" -mreiser");
                actualMethod = "reiser";
            } else if ("rolia".equals(ms)) {
                cmd.append(" -mrolia");
                actualMethod = "rolia";
            } else if ("zhou".equals(ms)) {
                cmd.append(" -mzhou");
                actualMethod = "zhou";
            }
        }
        cmd.append(" -o ").append(resultFile.getAbsolutePath());
        cmd.append(" > ").append(logFile.getAbsolutePath()).append(" 2>&1");
        return cmd.toString();
    }

    private boolean hasMultiServer() {
        for (int i = 0; i < sn.nstations; i++) {
            int servers = (int) sn.nservers.get(i);
            if (servers > 1 && servers != Integer.MAX_VALUE) return true;
        }
        return false;
    }

    private boolean isWindows() {
        return System.getProperty("os.name").toLowerCase().contains("win");
    }

    private ParsedResults parseResults(File resultFile, int nchains) {
        Matrix Uchain = new Matrix(sn.nstations, nchains);
        Matrix Qchain = new Matrix(sn.nstations, nchains);
        Matrix Wchain = new Matrix(sn.nstations, nchains);
        Matrix Tchain = new Matrix(sn.nstations, nchains);
        Uchain.fill(0.0); Qchain.fill(0.0); Wchain.fill(0.0); Tchain.fill(0.0);

        if (!resultFile.exists()) {
            throw new RuntimeException("QNS result file not found: " + resultFile.getAbsolutePath());
        }
        try (java.io.BufferedReader reader = new java.io.BufferedReader(new java.io.FileReader(resultFile))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.contains(",") && !line.contains("$")) {
                    ParsedLine parsed = (sn.nclasses == 1)
                            ? parseDollarOutputSingleClass(line, nchains)
                            : parseDollarOutput(line, nchains);
                    if (parsed != null) {
                        int stationIdx = -1;
                        for (int i = 0; i < sn.nstations; i++) {
                            int nodeIdx = (int) sn.stationToNode.get(i);
                            String stationNodeName = sn.nodenames.get(nodeIdx);
                            if (stationNodeName.equals(parsed.statName)) { stationIdx = i; break; }
                        }
                        if (stationIdx != -1) {
                            for (int c = 0; c < nchains; c++) {
                                Qchain.set(stationIdx, c, parsed.Q[c]);
                                Wchain.set(stationIdx, c, parsed.W[c]);
                                Uchain.set(stationIdx, c, parsed.U[c]);
                                Tchain.set(stationIdx, c, parsed.T[c]);
                            }
                        }
                    }
                }
            }
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
        return new ParsedResults(Uchain, Qchain, Wchain, Tchain);
    }

    private ParsedLine parseDollarOutput(String line, int nchains) {
        String[] parts = line.replace(" ", "").split(",");
        if (parts.length < 1 + 4 * (nchains + 1)) return null;
        String statName = parts[0];
        double[] Q = new double[nchains];
        double[] W = new double[nchains];
        double[] U = new double[nchains];
        double[] T = new double[nchains];
        int ptr = 1;
        for (int r = 0; r < nchains; r++) Q[r] = parseDoubleOrZero(parts[ptr + r]);
        ptr += nchains + 1;
        for (int r = 0; r < nchains; r++) W[r] = parseDoubleOrZero(parts[ptr + r]);
        ptr += nchains + 1;
        for (int r = 0; r < nchains; r++) U[r] = parseDoubleOrZero(parts[ptr + r]);
        ptr += nchains + 1;
        for (int r = 0; r < nchains; r++) T[r] = parseDoubleOrZero(parts[ptr + r]);
        return new ParsedLine(statName, Q, W, U, T);
    }

    private ParsedLine parseDollarOutputSingleClass(String line, int nchains) {
        String[] parts = line.replace(" ", "").split(",");
        if (parts.length < 5) return null;
        String statName = parts[0];
        double[] Q = new double[nchains];
        double[] W = new double[nchains];
        double[] U = new double[nchains];
        double[] T = new double[nchains];
        int ptr = 1;
        for (int r = 0; r < nchains; r++) Q[r] = parseDoubleOrZero(parts[ptr + r]);
        ptr += 1;
        for (int r = 0; r < nchains; r++) W[r] = parseDoubleOrZero(parts[ptr + r]);
        ptr += 1;
        for (int r = 0; r < nchains; r++) U[r] = parseDoubleOrZero(parts[ptr + r]);
        ptr += 1;
        for (int r = 0; r < nchains; r++) T[r] = parseDoubleOrZero(parts[ptr + r]);
        return new ParsedLine(statName, Q, W, U, T);
    }

    private static double parseDoubleOrZero(String s) {
        try { return Double.parseDouble(s); } catch (Exception e) { return 0.0; }
    }

    private boolean deleteDirectory(File dir) {
        if (dir.isDirectory()) {
            File[] children = dir.listFiles();
            if (children != null) for (File child : children) deleteDirectory(child);
        }
        return dir.delete();
    }

    private static final class ParsedLine {
        final String statName;
        final double[] Q;
        final double[] W;
        final double[] U;
        final double[] T;
        ParsedLine(String statName, double[] Q, double[] W, double[] U, double[] T) {
            this.statName = statName; this.Q = Q; this.W = W; this.U = U; this.T = T;
        }
    }

    private static final class ParsedResults {
        final Matrix Uchain;
        final Matrix Qchain;
        final Matrix Wchain;
        final Matrix Tchain;
        ParsedResults(Matrix Uchain, Matrix Qchain, Matrix Wchain, Matrix Tchain) {
            this.Uchain = Uchain; this.Qchain = Qchain; this.Wchain = Wchain; this.Tchain = Tchain;
        }
    }
}
