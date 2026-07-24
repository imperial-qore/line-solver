package jline.solvers.ln;

import jline.TestTools;
import jline.solvers.LayeredNetworkAvgTable;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.logging.Level;
import java.util.logging.Logger;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Abstract base class for SolverLN tests, providing shared utilities.
 */
abstract class SolverLNTestBase {

    static {
        // Suppress all logging during tests
        Logger.getLogger("").setLevel(Level.OFF);
        Logger.getLogger("jline").setLevel(Level.OFF);
        Logger.getLogger("org.apache.commons.io.FileUtils").setLevel(Level.OFF);
    }

    protected static final boolean COMPARE_WITH_LQSIM = false;

    // Store original streams for restoration
    private final PrintStream originalOut = System.out;
    private final PrintStream originalErr = System.err;

    /**
     * Suppress System.out and System.err during execution of a Runnable
     */
    protected void suppressOutput(Runnable action) {
        ByteArrayOutputStream devNull = new ByteArrayOutputStream();
        PrintStream nullStream = new PrintStream(devNull);

        // Suppress Java logger output
        Logger rootLogger = Logger.getLogger("");
        Logger jlineLogger = Logger.getLogger("jline");
        Logger fileUtilsLogger = Logger.getLogger("org.apache.commons.io.FileUtils");
        Level originalRootLevel = rootLogger.getLevel();
        Level originalJlineLevel = jlineLogger.getLevel();
        Level originalFileUtilsLevel = fileUtilsLogger.getLevel();

        try {
            System.setOut(nullStream);
            System.setErr(nullStream);
            rootLogger.setLevel(Level.OFF);
            jlineLogger.setLevel(Level.OFF);
            fileUtilsLogger.setLevel(Level.OFF);
            action.run();
        } finally {
            System.setOut(originalOut);
            System.setErr(originalErr);
            rootLogger.setLevel(originalRootLevel);
            jlineLogger.setLevel(originalJlineLevel);
            fileUtilsLogger.setLevel(originalFileUtilsLevel);
            nullStream.close();
        }
    }

    /**
     * Execute a function with suppressed output and return its result
     */
    protected <T> T suppressOutput(java.util.function.Supplier<T> supplier) {
        ByteArrayOutputStream devNull = new ByteArrayOutputStream();
        PrintStream nullStream = new PrintStream(devNull);

        // Suppress Java logger output
        Logger rootLogger = Logger.getLogger("");
        Logger jlineLogger = Logger.getLogger("jline");
        Logger fileUtilsLogger = Logger.getLogger("org.apache.commons.io.FileUtils");
        Level originalRootLevel = rootLogger.getLevel();
        Level originalJlineLevel = jlineLogger.getLevel();
        Level originalFileUtilsLevel = fileUtilsLogger.getLevel();

        try {
            System.setOut(nullStream);
            System.setErr(nullStream);
            rootLogger.setLevel(Level.OFF);
            jlineLogger.setLevel(Level.OFF);
            fileUtilsLogger.setLevel(Level.OFF);
            return supplier.get();
        } finally {
            System.setOut(originalOut);
            System.setErr(originalErr);
            rootLogger.setLevel(originalRootLevel);
            jlineLogger.setLevel(originalJlineLevel);
            fileUtilsLogger.setLevel(originalFileUtilsLevel);
            nullStream.close();
        }
    }

    protected void callAssertTableMetrics(LayeredNetworkAvgTable avgTable, double[][] expected) {
        if (expected == null) {
            assertNotNull(avgTable);
            return;
        }

        int nRows = expected.length;
        double[] Q = new double[nRows];
        double[] U = new double[nRows];
        double[] R = new double[nRows];
        double[] Res = new double[nRows];
        double[] Arv = new double[nRows];
        double[] T = new double[nRows];

        for (int i = 0; i < nRows; i++) {
            Q[i] = expected[i][0];
            U[i] = expected[i][1];
            R[i] = expected[i][2];
            Res[i] = expected[i][3];
            Arv[i] = expected[i][4];
            T[i] = expected[i][5];
        }

        TestTools.assertTableMetrics(avgTable, Q, U, R, Res, Arv, T);
    }
}
