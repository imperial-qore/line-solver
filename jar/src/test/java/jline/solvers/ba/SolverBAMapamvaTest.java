package jline.solvers.ba;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.Test;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.MAP;
import jline.solvers.NetworkAvgTable;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;

/**
 * MAP-AMVA LP bounds (Casale-Smirni, IEEE/IFIP DSN 2009) in SolverBA.
 *
 * <p>The family is the only one in SolverBA derived FOR a correlated service
 * process, so these tests pin the two things that distinguish it from its
 * neighbours: that its bracket actually contains the exact solution of a MAP
 * model, and that the feature gate lets a MAP reach it and no other family.
 */
public class SolverBAMapamvaTest {

    /** D0 + D1 is a proper generator; D1 has an off-diagonal entry, so successive services correlate. */
    private static Matrix mapD0() {
        Matrix d0 = new Matrix(2, 2);
        d0.set(0, 0, -3.0);
        d0.set(0, 1, 0.5);
        d0.set(1, 0, 0.2);
        d0.set(1, 1, -0.4);
        return d0;
    }

    private static Matrix mapD1() {
        Matrix d1 = new Matrix(2, 2);
        d1.set(0, 0, 2.5);
        d1.set(0, 1, 0.0);
        d1.set(1, 0, 0.2);
        d1.set(1, 1, 0.0);
        return d1;
    }

    private static Network tandem(boolean mapFirst, int njobs) {
        Network model = new Network("mapamvaTandem");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(model, "C", njobs, q1);
        if (mapFirst) {
            q1.setService(c, new MAP(mapD0(), mapD1()));
            q2.setService(c, new Exp(2.0));
        } else {
            q1.setService(c, new Exp(2.0));
            q2.setService(c, new MAP(mapD0(), mapD1()));
        }
        model.link(model.serialRouting(q1, q2));
        return model;
    }

    private static List<List<Double>> metrics(NetworkAvgTable t) {
        return Arrays.asList(t.getQLen(), t.getUtil(), t.getRespT(), t.getTput());
    }

    private void assertBrackets(Network model) {
        List<List<Double>> lo = metrics(new SolverBA(model, "mapamva.lower").getAvgTable());
        List<List<Double>> up = metrics(new SolverBA(model, "mapamva.upper").getAvgTable());
        List<List<Double>> ex = metrics(new SolverCTMC(model).getAvgTable());
        String[] names = {"QLen", "Util", "RespT", "Tput"};
        for (int c = 0; c < names.length; c++) {
            for (int i = 0; i < ex.get(c).size(); i++) {
                double l = lo.get(c).get(i), u = up.get(c).get(i), e = ex.get(c).get(i);
                assertTrue(l <= e + 1e-9,
                        names[c] + " station " + i + ": lower " + l + " above exact " + e);
                assertTrue(e <= u + 1e-9,
                        names[c] + " station " + i + ": exact " + e + " above upper " + u);
            }
        }
    }

    @Test
    public void aMapAtTheLastStationIsBracketed() {
        assertBrackets(tandem(false, 8));
    }

    @Test
    public void aMapAtTheFirstStationIsBracketed() {
        // The LP requires the phase-carrying queue to be index M; the analyzer
        // PERMUTES it there and inverts the permutation before reporting, so the
        // bracket must hold with the stations the other way round too.
        assertBrackets(tandem(true, 8));
    }

    @Test
    public void anAllExponentialModelDegeneratesToOneLevel() {
        // No phase-carrying station: K collapses to 1 and the balances become the
        // product-form ones. Still a bracket, and the exact solution is in it.
        Network model = new Network("mapamvaExp");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(model, "C", 5, q1);
        q1.setService(c, new Exp(2.0));
        q2.setService(c, new Exp(1.0));
        model.link(model.serialRouting(q1, q2));
        assertBrackets(model);
    }

    @Test
    public void aMapModelReachesMapamvaAndNothingElse() {
        // "MAP" sits in the BASE feature envelope so that this family can accept
        // it, and getMethodFeatureSet takes it back from every other: each of the
        // rest reads the service mean alone and would bracket a DIFFERENT system.
        SolverBA solver = new SolverBA(tandem(false, 4));
        List<String> supported = new ArrayList<String>();
        for (String m : solver.listValidMethods()) {
            if (solver.supportsModelMethod(m).isEmpty()) {
                supported.add(m);
            }
        }
        Collections.sort(supported);
        assertEquals(Arrays.asList("mapamva.lower", "mapamva.upper"), supported);
    }

    @Test
    public void aRenewalModelDoesNotLoseItsOwnFamilies() {
        // The mirror check: taking MAP away from the others must not take
        // anything else with it.
        Network model = new Network("mapamvaRenewal");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(model, "C", 5, q1);
        q1.setService(c, new Exp(2.0));
        q2.setService(c, Erlang.fitMeanAndOrder(1.0, 2));
        model.link(model.serialRouting(q1, q2));
        SolverBA solver = new SolverBA(model);
        List<String> supported = new ArrayList<String>();
        for (String m : solver.listValidMethods()) {
            if (solver.supportsModelMethod(m).isEmpty()) {
                supported.add(m);
            }
        }
        for (String name : new String[]{"gb.upper", "gb.lower", "aba.upper", "mapamva.upper"}) {
            assertTrue(supported.contains(name), name + " is no longer offered");
        }
    }

    @Test
    public void aDelayStationIsRefused() {
        // The LP is a network of queues, and Casale-Smirni name the delay
        // extension as open work.
        Network model = new Network("mapamvaDelay");
        Delay d = new Delay(model, "Think");
        Queue q = new Queue(model, "Q1", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(model, "C", 4, d);
        d.setService(c, new Exp(1.0));
        q.setService(c, new MAP(mapD0(), mapD1()));
        model.link(model.serialRouting(d, q));
        assertTrue(!Arrays.asList(new SolverBA(model).listValidMethods()).contains("mapamva.upper"));
        assertThrows(Exception.class, () -> new SolverBA(model, "mapamva.upper").getAvgTable());
    }

    @Test
    public void twoPhaseCarryingStationsAreRefused() {
        // The LP gives queue M the (D0,D1) pair and every other queue a SCALAR
        // rate, so it has nowhere to put a second phase-type station.
        Network model = new Network("mapamvaTwoPh");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(model, "C", 4, q1);
        q1.setService(c, Erlang.fitMeanAndOrder(0.5, 2));
        q2.setService(c, new MAP(mapD0(), mapD1()));
        model.link(model.serialRouting(q1, q2));
        assertThrows(Exception.class, () -> new SolverBA(model, "mapamva.upper").getAvgTable());
    }
}
