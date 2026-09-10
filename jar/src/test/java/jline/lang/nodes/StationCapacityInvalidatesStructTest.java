package jline.lang.nodes;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Exp;
import jline.solvers.ag.SolverAG;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * A capacity set AFTER the struct was materialized must still be honoured.
 * <p>
 * {@code Network} caches the compiled {@code sn} behind {@code hasStruct} and
 * rebuilds it only when that flag is down. {@code sn.cap}/{@code sn.classcap} are
 * DERIVED from the station's declared capacity, so a {@code setCapacity} that does
 * not invalidate is not applied late -- it is DROPPED, and every sn-reading solver
 * then answers the UNBOUNDED model with no error and a plausible number.
 * </p><p>
 * The model is the smallest one that makes the difference visible: two FCFS queues
 * in a closed cycle, mu = 1 and 0.8, N = 2, buffer 1 at the second queue. Capped,
 * the chain has two states and is exact by hand --
 * (2,0) --1.0--&gt; (1,1) --0.8--&gt; (2,0), so pi = [0.8, 1.0]/1.8, X = 4/9 and
 * QLen = [13/9, 5/9]. Uncapped it is the product form over D = [1, 1.25]:
 * X = 0.590164 and QLen2 = 1.1475, MORE JOBS THAN THE BUFFER HOLDS, which is what
 * the stale struct used to return.
 * </p>
 */
public class StationCapacityInvalidatesStructTest {

    private static final double TOL = 1e-9;
    private static final double Q1_EXACT = 13.0 / 9.0;
    private static final double Q2_EXACT = 5.0 / 9.0;
    private static final double X_EXACT = 4.0 / 9.0;
    private static final double Q2_UNBOUNDED = 4.375 / 3.8125;

    /** cap &lt; 0 leaves the model uncapped; capLate defers the setter past getStruct(). */
    private static Network cycle(int cap, boolean capLate) {
        Network model = new Network("cycle");
        Queue q1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        ClosedClass c1 = new ClosedClass(model, "Class1", 2, q1, 0);
        q1.setService(c1, new Exp(1.0));
        q2.setService(c1, new Exp(0.8));
        if (cap >= 0 && !capLate) {
            q2.setCapacity(cap);
        }
        model.link(model.serialRouting(q1, q2));
        if (cap >= 0 && capLate) {
            model.getStruct();     // materialize the cache FIRST
            q2.setCapacity(cap);   // the setter under test
        }
        return model;
    }

    @Test
    public void lateCapacityReachesTheStruct() {
        NetworkStruct sn = cycle(1, true).getStruct();
        assertEquals(1.0, sn.cap.get(1, 0), TOL);
        assertEquals(1.0, sn.classcap.get(1, 0), TOL);
    }

    @Test
    public void ctmcHonoursALateCapacity() {
        for (boolean late : new boolean[]{false, true}) {
            SolverCTMC solver = new SolverCTMC(cycle(1, late));
            Matrix QN = solver.getAvgQLen();
            Matrix TN = solver.getAvgTput();
            assertEquals(Q1_EXACT, QN.get(0, 0), 1e-8, "QLen1, capLate=" + late);
            assertEquals(Q2_EXACT, QN.get(1, 0), 1e-8, "QLen2, capLate=" + late);
            assertEquals(X_EXACT, TN.get(0, 0), 1e-8, "Tput, capLate=" + late);
            // and specifically NOT the unbounded answer a stale struct gave
            assertTrue(Math.abs(QN.get(1, 0) - Q2_UNBOUNDED) > 1e-3,
                    "returned the UNBOUNDED answer, capLate=" + late);
            // a station never holds more than its buffer
            assertTrue(QN.get(1, 0) <= 1.0 + TOL);
        }
    }

    @Test
    public void numberOfServersReachesTheStruct() {
        Network model = cycle(-1, false);
        model.getStruct();
        ((Station) model.getNodeByName("Queue2")).setNumberOfServers(3);
        assertEquals(3.0, model.getStruct().nservers.get(1, 0), TOL);
    }

    /**
     * RCAT has no representation of a finite buffer: before the gate SolverAG
     * returned the same figures with and without the cap.
     */
    @Test
    public void agRefusesABindingCapacity() {
        String reason = new SolverAG(cycle(1, false)).supportsModelMethod("inap");
        assertFalse(reason.isEmpty(), "SolverAG accepted a binding finite capacity");
        assertTrue(reason.contains("apacity"), reason);
    }

    /** Only a capacity that can BIND is refused. */
    @Test
    public void agStillSolvesTheUncappedModel() {
        String reason = new SolverAG(cycle(-1, false)).supportsModelMethod("inap");
        assertTrue(reason.isEmpty(), reason);
    }
}
