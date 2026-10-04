package jline.solvers.wrappers.lqns.qnsolver;

import java.io.File;

import jline.api.sn.SnGetArvRFromTput;
import jline.api.sn.SnGetResidTFromRespT;
import jline.lang.NetworkStruct;
import jline.solvers.AvgHandle;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.wrappers.lqns.SolverLQNS;
import jline.util.matrix.Matrix;

/**
 * Analyzer for the qnsolver path of SolverLQNS (its qns methods on a flat Network).
 */
public class Solver_qns_analyzer {
    private final SolverLQNS solver;

    public Solver_qns_analyzer(SolverLQNS solver) {
        this.solver = solver;
    }

    /**
     * Check if a native {@code qnsolver} binary is available on the PATH.
     */
    public static boolean hasNativeQNSolver() {
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
     * Check if qnsolver can be run: a native binary on the PATH. LINE never runs
     * it from a container image, since qnsolver ships with LQNS and that licence
     * forbids redistribution.
     */
    public static boolean isQNSolverAvailable() {
        return hasNativeQNSolver();
    }

    /**
     * Run the qnsolver analyzer.
     */
    public SolverResult runAnalyzer() {
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
        SolverResult result = handler.solve();

        AvgHandle T = solver.getAvgTputHandles();
        Matrix AN = SnGetArvRFromTput.snGetArvRFromTput(sn, result.TN, T);
        // The residence time is the response time TIMES THE VISITS. The handler
        // returns RN in both slots, and copying it through reported ResidT ==
        // RespT at every station visited more than once per system passage.
        // MATLAB leaves WN empty here and lets getAvg derive it the same way.
        Matrix WN = SnGetResidTFromRespT.snGetResidTFromRespT(sn, result.RN,
                solver.getAvgResidTHandles());

        long endTime = System.nanoTime();
        double runtime = (endTime - startTime) / 1e9;

        String actualMethod;
        if ("default".equals(options.method)) {
            actualMethod = "default/" + result.method;
        } else {
            actualMethod = result.method;
        }

        return Solver_qns.result(
                result.QN,
                result.UN,
                result.RN,
                result.TN,
                AN,
                WN,
                result.CN,
                result.XN,
                runtime,
                actualMethod,
                0);
    }
}
