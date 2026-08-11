package jline.solvers.mam;

import jline.io.InputOutput;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.lang.processes.ME;
import jline.util.matrix.Matrix;
import jline.GlobalConstants;
import jline.VerboseLevel;
import org.apache.commons.io.FileUtils;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import java.util.logging.Logger;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * A multiclass station whose service is matrix-exponential falls outside the
 * exact RAP/RAP/1 QBD and is answered by the phase-type approximation
 * MMAPPH1FCFS, which is not exact for a non-phase-type service. The user is
 * told so, and must be told so for every model concerned: this warning reports
 * a correctness limitation, so suppressing it as a repeat would let the second
 * and later models in a session read as clean.
 */
public class MamMEWarningTest {

    private VerboseLevel savedVerbose;

    // Surefire runs the suite with -Djline.verbose=SILENT, under which
    // line_warning_always correctly emits nothing. Restore the default STD
    // verbosity for these tests, which are about what a user sees.
    @BeforeEach
    public void raiseVerbosity() {
        savedVerbose = GlobalConstants.Verbose;
        GlobalConstants.Verbose = VerboseLevel.STD;
    }

    @AfterEach
    public void restoreVerbosity() {
        GlobalConstants.Verbose = savedVerbose;
    }

    /** Collects the records emitted through InputOutput's logger. */
    private static class CapturingHandler extends Handler {
        final List<String> messages = new ArrayList<String>();

        @Override
        public void publish(LogRecord record) {
            messages.add(record.getMessage());
        }

        @Override
        public void flush() {
        }

        @Override
        public void close() {
        }
    }

    /**
     * A matrix-exponential with a negative entry in alpha, whose density has an
     * interior zero, so it admits no phase-type representation of any order.
     * Same representation as MEDistributionTest.nonPhaseTypeME.
     */
    private static ME nonPhaseTypeME() {
        Matrix alpha = new Matrix(new double[]{
                0.61058991931158258, -0.15547146730086722, 0.54488154798928464});
        Matrix A = new Matrix(new double[][]{
                {-1.0, 0.0, 0.0},
                {0.0, -2.0, 2.0},
                {0.0, -2.0, -2.0}
        });
        return new ME(alpha, A);
    }

    private static Network twoClassMeQueue() {
        Network model = new Network("M/ME/1 two classes");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass class1 = new OpenClass(model, "Class1", 0);
        OpenClass class2 = new OpenClass(model, "Class2", 0);
        source.setArrival(class1, Exp.fitRate(0.2));
        source.setArrival(class2, Exp.fitRate(0.2));
        queue.setService(class1, nonPhaseTypeME());
        queue.setService(class2, nonPhaseTypeME());
        RoutingMatrix rm = model.initRoutingMatrix();
        rm.addConnection(class1, class1, source, queue, 1.0);
        rm.addConnection(class1, class1, queue, sink, 1.0);
        rm.addConnection(class2, class2, source, queue, 1.0);
        rm.addConnection(class2, class2, queue, sink, 1.0);
        model.link(rm);
        return model;
    }

    private static List<String> solveCapturingWarnings(int numSolves) {
        Logger logger = Logger.getLogger(FileUtils.class.getName());
        CapturingHandler handler = new CapturingHandler();
        handler.setLevel(Level.ALL);
        logger.addHandler(handler);
        try {
            for (int i = 0; i < numSolves; i++) {
                new SolverMAM(twoClassMeQueue()).getAvgTable();
            }
        } finally {
            logger.removeHandler(handler);
        }
        return handler.messages;
    }

    private static int countMEWarnings(List<String> messages) {
        int n = 0;
        for (String m : messages) {
            if (m != null && m.contains("matrix-exponential or rational service process")) {
                n++;
            }
        }
        return n;
    }

    @Test
    public void testMEMulticlassWarningIsRaised() {
        List<String> messages = solveCapturingWarnings(1);
        assertEquals(1, countMEWarnings(messages));

        String warning = null;
        for (String m : messages) {
            if (m != null && m.contains("matrix-exponential or rational service process")) {
                warning = m;
            }
        }
        // The caller resolves: mfilename(new Object() {}) reads the enclosing
        // method, which exists here because the anonymous class is declared
        // inside solver_mam_basic and not inside a constructor.
        assertTrue(warning.startsWith("[solver_mam_basic] "), warning);
        // The station is named, not indexed: station indices are 1-based in
        // MATLAB and 0-based in the JAR and in Python, so an index would make
        // the same message read differently in each codebase.
        assertTrue(warning.contains("Station Queue has a matrix-exponential "
                + "or rational service process"), warning);
        assertTrue(warning.contains("here 2 classes, 1 servers"), warning);
        assertTrue(warning.contains("not exact for this service process"), warning);
    }

    /**
     * The fixed-point iteration visits the station several times per solve; the
     * user needs the limitation stated once, not once per sweep.
     */
    @Test
    public void testMEWarningRaisedOncePerSolveNotPerIteration() {
        assertEquals(1, countMEWarnings(solveCapturingWarnings(1)));
    }

    /**
     * Three models solved back to back produce three warnings. line_warning
     * would have collapsed these into one plus a suppression notice, leaving the
     * second and third models looking clean.
     */
    @Test
    public void testMEWarningRepeatsAcrossModels() {
        assertEquals(3, countMEWarnings(solveCapturingWarnings(3)));
    }

    /**
     * line_warning_always is exempt from the repeat suppression that
     * line_warning applies, but not from the verbosity gate.
     */
    @Test
    public void testLineWarningAlwaysDoesNotSuppressRepeats() {
        Logger logger = Logger.getLogger(FileUtils.class.getName());
        CapturingHandler handler = new CapturingHandler();
        handler.setLevel(Level.ALL);
        logger.addHandler(handler);
        try {
            for (int i = 0; i < 3; i++) {
                InputOutput.line_warning_always("MamMEWarningTest", "identical message %d", 7);
            }
        } finally {
            logger.removeHandler(handler);
        }
        int n = 0;
        for (String m : handler.messages) {
            if (m != null && m.equals("[MamMEWarningTest] identical message 7")) {
                n++;
            }
        }
        assertEquals(3, n);
    }
}
