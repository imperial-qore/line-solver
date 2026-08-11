package jline.solvers.wrappers.qns.analyzers;

import java.io.File;

import jline.api.sn.SnGetArvRFromTput;
import jline.lang.NetworkStruct;
import jline.solvers.AvgHandle;
import jline.solvers.SolverOptions;
import jline.solvers.wrappers.qns.QNSResult;
import jline.solvers.wrappers.qns.SolverQNS;
import jline.solvers.wrappers.qns.handlers.Solver_qns;
import jline.util.matrix.Matrix;

/**
 * Analyzer for the QNS solver.
 */
public class Solver_qns_analyzer {
    private final SolverQNS solver;

    public Solver_qns_analyzer(SolverQNS solver) {
        this.solver = solver;
    }

    /**
     * Check if the external qnsolver tool is available.
     */
    public static boolean isQNSolverAvailable() {
        try {
            File devNull = new File(System.getProperty("os.name").toLowerCase().contains("win") ? "NUL" : "/dev/null");
            ProcessBuilder pb = new ProcessBuilder("qnsolver", "--help");
            pb.redirectOutput(devNull);
            pb.redirectError(devNull);
            Process process = pb.start();
            process.waitFor();
            return true;
        } catch (Exception e) {
            return false;
        }
    }

    /**
     * Run the QNS analyzer.
     */
    public QNSResult runAnalyzer() {
        long startTime = System.nanoTime();

        NetworkStruct sn = solver.sn;
        SolverOptions options = solver.options;

        if (options.method == null) {
            options.method = "default";
        }

        if ("conway".equals(options.method)) {
            options.config.multiserver = "conway";
        } else if ("rolia".equals(options.method)) {
            options.config.multiserver = "rolia";
        } else if ("zhou".equals(options.method)) {
            options.config.multiserver = "zhou";
        } else if ("suri".equals(options.method)) {
            options.config.multiserver = "suri";
        } else if ("reiser".equals(options.method)) {
            options.config.multiserver = "reiser";
        } else if ("schmidt".equals(options.method)) {
            options.config.multiserver = "schmidt";
        } else if ("default".equals(options.method)) {
            options.config.multiserver = "default";
        }

        Solver_qns handler = new Solver_qns(sn, options);
        QNSResult result = handler.solve();

        AvgHandle T = solver.getAvgTputHandles();
        Matrix AN = SnGetArvRFromTput.snGetArvRFromTput(sn, result.TN, T);

        long endTime = System.nanoTime();
        double runtime = (endTime - startTime) / 1e9;

        String actualMethod;
        if ("default".equals(options.method)) {
            actualMethod = "default/" + result.method;
        } else {
            actualMethod = result.method;
        }

        return new QNSResult(
                result.QN,
                result.UN,
                result.RN,
                result.TN,
                AN,
                result.WN,
                result.CN,
                result.XN,
                runtime,
                actualMethod,
                0);
    }
}
