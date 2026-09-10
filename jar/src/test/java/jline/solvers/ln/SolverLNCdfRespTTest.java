package jline.solvers.ln;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.Activity;
import jline.lang.layered.Entry;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Processor;
import jline.lang.layered.Task;
import jline.lang.processes.Exp;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

/**
 * {@link SolverLN#getCdfRespT()}: the per-entry response time distribution.
 *
 * <p>The law exists only under the {@code moment3} update pass, which fits an APH
 * to the first three moments of every activity and call term and convolves them.
 * The getter re-runs the ensemble under {@code moment3} when the solver was built
 * for another method, exactly as MATLAB {@code @SolverLN/getCdfRespT.m} does, and
 * returns one (n x 2) [F(t), t] table per entry in the entry-local index space.</p>
 *
 * <p>Entry 0 is the case that used to be lost: the port carried MATLAB's 1-based
 * entry numbering into a 0-based index space and wrote the first entry's law to
 * slot -1, and it kept F(0) -- zero by construction -- rather than the table.</p>
 */
class SolverLNCdfRespTTest {

    @BeforeAll
    public static void setUp() {
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    /** T1 -> T2 -> T3, so the middle task is at once a server and a caller. */
    private static LayeredNetwork serialChain() {
        LayeredNetwork m = new LayeredNetwork("serialchain");
        Processor p1 = new Processor(m, "P1", 1, SchedStrategy.PS);
        Processor p2 = new Processor(m, "P2", 1, SchedStrategy.PS);
        Processor p3 = new Processor(m, "P3", 1, SchedStrategy.PS);
        Task t1 = new Task(m, "T1", 4, SchedStrategy.REF).on(p1).setThinkTime(new Exp(1));
        Task t2 = new Task(m, "T2", 2, SchedStrategy.FCFS).on(p2);
        Task t3 = new Task(m, "T3", 1, SchedStrategy.FCFS).on(p3);
        Entry e1 = new Entry(m, "E1").on(t1);
        Entry e2 = new Entry(m, "E2").on(t2);
        Entry e3 = new Entry(m, "E3").on(t3);
        new Activity(m, "A1", new Exp(5)).on(t1).boundTo(e1).synchCall(e2, 1);
        new Activity(m, "A2", new Exp(5)).on(t2).boundTo(e2).synchCall(e3, 1).repliesTo(e2);
        new Activity(m, "A3", new Exp(5)).on(t3).boundTo(e3).repliesTo(e3);
        return m;
    }

    private static SolverLN solver(String method) {
        SolverOptions o = new SolverOptions();
        if (method != null) {
            o.method = method;
        }
        o.verbose = VerboseLevel.SILENT;
        return new SolverLN(serialChain(), o);
    }

    /** [F(t), t], NOT [t, F(t)]: column 0 is the probability, column 1 the time. */
    private static void assertIsCdfTable(Matrix cdf) {
        assertNotNull(cdf);
        assertEquals(2, cdf.getNumCols());
        assertTrue(cdf.getNumRows() > 1);
        assertEquals(0.0, cdf.get(0, 1), 1e-12);
        assertEquals(0.0, cdf.get(0, 0), 1e-9);
        for (int i = 1; i < cdf.getNumRows(); i++) {
            assertTrue(cdf.get(i, 1) >= cdf.get(i - 1, 1), "time column must be nondecreasing");
            assertTrue(cdf.get(i, 0) >= cdf.get(i - 1, 0) - 1e-9, "CDF column must be nondecreasing");
            assertTrue(cdf.get(i, 0) <= 1 + 1e-12);
        }
        assertTrue(cdf.get(cdf.getNumRows() - 1, 0) > 0.99, "the grid must reach into the tail");
    }

    @Test
    void moment3ReturnsOneLawPerEntry() {
        SolverLN s = solver("moment3");
        s.getAvgTable();
        List<Matrix> RD = s.getCdfRespT();
        assertEquals(3, RD.size());
        for (int e = 0; e < RD.size(); e++) {
            assertIsCdfTable(RD.get(e));
        }
    }

    /**
     * The tabulated law and the entry's own mean are the SAME quantity. They part
     * company when the interlock rescale of updatePopulations is let at the entry
     * servt after the moment pass has formed it -- BUG-97.
     */
    @Test
    void lawAgreesWithTheEntryMean() {
        SolverLN s = solver("moment3");
        s.getAvgTable();
        List<Matrix> RD = s.getCdfRespT();
        for (int e = 0; e < RD.size(); e++) {
            Matrix cdf = RD.get(e);
            // E[X] = int (1-F) dt over the tabulated grid, trapezoidal.
            double mean = 0;
            for (int i = 1; i < cdf.getNumRows(); i++) {
                double dt = cdf.get(i, 1) - cdf.get(i - 1, 1);
                mean += 0.5 * dt * ((1 - cdf.get(i, 0)) + (1 - cdf.get(i - 1, 0)));
            }
            double aphMean = s.entryproc.get(e).getMean();
            assertEquals(aphMean, mean, 2e-2 * aphMean);
        }
        // E3 holds Exp(5) alone and waits on nothing below it.
        assertTrue(s.entryproc.get(2).getMean() >= 0.2);
    }

    /** The mean-based update forms no law, so the getter flips the method itself. */
    @Test
    void getterRerunsUnderMoment3AndRestoresTheMethod() {
        SolverLN s = solver("srvn.cs");
        s.getAvgTable();
        assertEquals("srvn.cs", s.lnmethod);
        List<Matrix> RD = s.getCdfRespT();
        assertIsCdfTable(RD.get(0));
        assertEquals("srvn.cs", s.lnmethod, "the caller's method must be put back");
    }

    /** srvn.ph carries no activity-graph routing, so the getter refuses by name. */
    @Test
    void phEncodingIsRefusedByName() {
        SolverLN s = solver("srvn.ph");
        s.getAvgTable();
        if (!s.isPHEncoding()) {
            return; // the alias fell back to a routing encoding: nothing to refuse
        }
        RuntimeException ex = assertThrows(RuntimeException.class, () -> s.getCdfRespT());
        assertTrue(ex.getMessage().contains("routing encoding"), ex.getMessage());
    }
}
