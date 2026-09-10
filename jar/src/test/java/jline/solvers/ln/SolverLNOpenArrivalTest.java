package jline.solvers.ln;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.examples.java.basic.LayeredModel;
import jline.lang.layered.LayeredNetwork;
import jline.solvers.LayeredNetworkAvgTable;

import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

/**
 * Entry-level open arrivals in SolverLN.
 *
 * <p>An entry may take an exogenous Poisson stream that is not a call
 * (Entry.setArrival). Where the entry ALSO has a caller, the stream rides an
 * OpenClass through the entry's host layer in one pass, walking no activity
 * graph. Where the arrival is the ONLY way into the task, as here, there is no
 * open class at all: that task has no task layer, so nothing would ever set the
 * surrogate delay of its caller chain, and a stream on top of that unthrottled
 * chain loads the host twice. SolverLN closes the chain on the known arrival
 * rate instead, the construction a forwarding target gets.</p>
 *
 * <p>The goldens below are lqns 6.2.28's, exactly, and lqsim (0.1996 / 0.319 /
 * 1.580) and LDES (0.19986 / 0.31957 / 1.599) confirm them; MATLAB, the C++ port
 * and native Python agree to 1e-6. They were rebased on 2026-08-11 from
 * 0.425 / 0.68 / 2.3529, which all four codebases produced and which this file
 * used to defend as "the layer decomposition's answer" -- the four agreed only
 * because three were ported from the fourth's construction. Twin of C++
 * cpp/tests/test_ln_openarrival.cpp and matlab
 * examples/basic/layeredModel/lqn_open_arrival.m.</p>
 */
class SolverLNOpenArrivalTest {

    private static final double TOL = 1e-6;

    @BeforeAll
    public static void setUp() {
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    private static double rowOf(LayeredNetworkAvgTable t, String nodeName, List<Double> col) {
        List<String> names = t.getNodeNames();
        int i = names.indexOf(nodeName);
        assertTrue(i >= 0, "row " + nodeName + " not in " + names);
        return col.get(i);
    }

    @Test
    void entryOpenArrivalMatchesReference() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_open_arrival();
        LayeredNetworkAvgTable t = (LayeredNetworkAvgTable) new SolverLN(model).getAvgTable();

        assertEquals(0.32, rowOf(t, "P1", t.getUtil()), TOL);
        assertEquals(0.32, rowOf(t, "T1", t.getUtil()), TOL);

        assertEquals(0.2, rowOf(t, "T1", t.getTput()), TOL);
        assertEquals(0.2, rowOf(t, "E1", t.getTput()), TOL);
        assertEquals(0.2, rowOf(t, "A1", t.getTput()), TOL);

        assertEquals(0.32, rowOf(t, "T1", t.getQLen()), TOL);
        assertEquals(0.32, rowOf(t, "E1", t.getQLen()), TOL);

        double residt = rowOf(t, "T1", t.getResidT());
        assertEquals(1.6, residt, 1e-6);
        assertEquals(1.6, rowOf(t, "E1", t.getRespT()), 1e-6);
        assertEquals(1.6, rowOf(t, "A1", t.getRespT()), 1e-6);

        // Little's law on the task chain, which is what the 0-metric defects broke
        assertEquals(rowOf(t, "T1", t.getQLen()), rowOf(t, "T1", t.getTput()) * residt, 1e-6);
    }

    @Test
    void forkPlusOpenStreamIsRefused() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_fork_open_arrival();
        // the refusal is the ROUTING encoding's: fj_mmt mints its own Source, which
        // detaches the one this layer already routes the stream through. 'srvn.ph'
        // composes the fork into the entry law and mints none, so it serves the
        // model and the 'srvn' alias no longer reaches this path -- hence the
        // encoding is named rather than taken. MATLAB refuses and serves it the
        // same way round; twin of cpp/tests/test_ln_fork_openarrival.cpp.
        LNOptions options = new LNOptions();
        options.method = "srvn.cs";
        RuntimeException e = assertThrows(RuntimeException.class,
                () -> new SolverLN(model, options).getAvgTable());
        assertTrue(e.getMessage().contains("carries both an AND fork and an open stream"),
                "unexpected message: " + e.getMessage());
    }

    @Test
    void forkPlusOpenStreamIsServedByThePhaseTypeEncoding() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_fork_open_arrival();
        SolverLN solver = new SolverLN(model);
        LayeredNetworkAvgTable t = (LayeredNetworkAvgTable) solver.getAvgTable();
        assertEquals("srvn.ph", solver.lnmethod);

        // MATLAB SolverLN with method='srvn.ph' on the same model, which the C++
        // port reproduces: the reference task, the server task and the open entry
        assertEquals(0.4, rowOf(t, "Client", t.getTput()), 4e-4);
        assertEquals(0.554545, rowOf(t, "Server", t.getTput()), 5.6e-4);
        assertEquals(0.110909, rowOf(t, "OE", t.getTput()), 1.2e-4);

        // the open entry takes its own arrival rate, the server carries both streams
        assertEquals(0.554545, rowOf(t, "P2", t.getUtil()), 5.6e-4);
        assertEquals(0.443636, rowOf(t, "SE", t.getTput()), 4.5e-4);
    }
}
