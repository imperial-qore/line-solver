package jline.solvers.ssa;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.Mode;
import jline.lang.nodes.Place;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.nodes.Transition;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates OPEN stochastic Petri nets (Source -> Place -> Transition -> Sink)
 * on the SERIAL (Gillespie) SSA engine. Before the fix the serial afterEvent
 * engine mis-measured Place markings: a Place is INF-scheduled, so
 * State.toMarginal read only the server slot while the FIRE-PRE relocated the
 * surviving tokens to the (invisible) buffer slot, and the consumption scaled
 * by the enabling degree instead of the arc weight. The result was garbage
 * marking means (M/M/inf collapsed, single-server read as infinite-server).
 *
 * Ground truth here is the exact closed form, which SolverJMT reproduces (the
 * cross-check against JMT is done in the MATLAB test): an M/M/1-as-SPN has mean
 * tokens rho/(1-rho); an M/M/inf-as-SPN has mean tokens lambda/mu; a tandem of
 * two single-server transitions is a Jackson network whose places are each
 * M/M/1 at the shared arrival rate.
 *
 * Each test asserts (a) the run used the serial engine (result.method ==
 * "serial", never a silent NRM fallback) and (b) the marking mean and
 * transition throughput match the closed form. The single-seed tolerance is 5%
 * relative, well above the noise floor at this sample count but tight enough to
 * fail on a systematic firing-accounting defect.
 */
public class SolverSSASpnOpenSerialTest {

    private static final int SAMPLES = 200000;
    private static final int SEED = 23000;
    private static final double RTOL = 0.05;

    /** Open net Source(Exp lambda) -> P1 -> T1(nservers, Exp mu) -> Sink. */
    private static Network singleQueue(double lambda, double mu, int nservers) {
        Network model = new Network("spn_open");
        Source source = new Source(model, "Source");
        Sink sink = new Sink(model, "Sink");
        Place P1 = new Place(model, "P1");
        Transition T1 = new Transition(model, "T1");
        OpenClass jobclass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobclass, new Exp(lambda));
        Mode m1 = T1.addMode("Mode1");
        T1.setNumberOfServers(m1, nservers);
        T1.setDistribution(m1, new Exp(mu));
        T1.setEnablingConditions(m1, jobclass, P1, 1);
        T1.setFiringOutcome(m1, jobclass, sink, 1);
        model.link(Network.serialRouting(source, P1, T1, sink));
        return model;
    }

    /** Tandem Source -> P1 -> T1(1 srv) -> P2 -> T2(1 srv) -> Sink (all Exp). */
    private static Network tandem(double lambda, double mu1, double mu2) {
        Network model = new Network("spn_tandem");
        Source source = new Source(model, "Source");
        Sink sink = new Sink(model, "Sink");
        Place P1 = new Place(model, "P1");
        Place P2 = new Place(model, "P2");
        Transition T1 = new Transition(model, "T1");
        Transition T2 = new Transition(model, "T2");
        OpenClass jobclass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobclass, new Exp(lambda));
        Mode m1 = T1.addMode("Mode1");
        T1.setNumberOfServers(m1, 1);
        T1.setDistribution(m1, new Exp(mu1));
        T1.setEnablingConditions(m1, jobclass, P1, 1);
        T1.setFiringOutcome(m1, jobclass, P2, 1);
        Mode m2 = T2.addMode("Mode2");
        T2.setNumberOfServers(m2, 1);
        T2.setDistribution(m2, new Exp(mu2));
        T2.setEnablingConditions(m2, jobclass, P2, 1);
        T2.setFiringOutcome(m2, jobclass, sink, 1);
        model.link(Network.serialRouting(source, P1, T1, P2, T2, sink));
        return model;
    }

    private static SolverSSA serial(Network model) {
        return new SolverSSA(model, "method", "serial", "samples", SAMPLES,
                "seed", SEED, "cutoff", 60);
    }

    private static void assertRanSerial(SolverSSA solver) {
        assertTrue(solver.result != null && solver.result.method != null
                        && solver.result.method.contains("serial"),
                "open SPN must run the serial method, got: "
                        + (solver.result == null ? "<no result>" : solver.result.method));
    }

    private static void assertClose(double expected, double got, String what) {
        double denom = Math.max(Math.abs(expected), 1e-9);
        double err = Math.abs(got - expected) / denom;
        assertTrue(err < RTOL, what + ": serial SSA " + got + " deviates from exact "
                + expected + " by " + (100.0 * err) + "%");
    }

    /** P1 is station index 1 (Source is 0). */
    @Test
    public void testMM1SingleServer() {
        double lambda = 0.5, mu = 1.0;
        double rho = lambda / mu;
        double qExact = rho / (1.0 - rho);   // = 1.0
        SolverSSA solver = serial(singleQueue(lambda, mu, 1));
        NetworkAvgTable t = solver.getAvgTable();
        assertRanSerial(solver);
        List<Double> q = t.getQLen();
        List<Double> tp = t.getTput();
        assertClose(qExact, q.get(1), "M/M/1 P1 mean tokens");
        assertClose(lambda, tp.get(1), "M/M/1 P1 throughput");
    }

    @Test
    public void testMM1HeavierLoad() {
        double lambda = 0.8, mu = 1.0;
        double rho = lambda / mu;
        double qExact = rho / (1.0 - rho);   // = 4.0
        SolverSSA solver = serial(singleQueue(lambda, mu, 1));
        NetworkAvgTable t = solver.getAvgTable();
        assertRanSerial(solver);
        List<Double> q = t.getQLen();
        List<Double> tp = t.getTput();
        assertClose(qExact, q.get(1), "M/M/1 rho=0.8 P1 mean tokens");
        assertClose(lambda, tp.get(1), "M/M/1 rho=0.8 P1 throughput");
    }

    @Test
    public void testMMInfiniteServer() {
        double lambda = 2.0, mu = 1.0;
        double qExact = lambda / mu;         // = 2.0
        SolverSSA solver = serial(singleQueue(lambda, mu, Integer.MAX_VALUE));
        NetworkAvgTable t = solver.getAvgTable();
        assertRanSerial(solver);
        List<Double> q = t.getQLen();
        List<Double> tp = t.getTput();
        assertClose(qExact, q.get(1), "M/M/inf P1 mean tokens");
        assertClose(lambda, tp.get(1), "M/M/inf P1 throughput");
    }

    @Test
    public void testTandemJackson() {
        double lambda = 0.5, mu1 = 1.0, mu2 = 0.8;
        double q1 = (lambda / mu1) / (1.0 - lambda / mu1);   // = 1.0
        double q2 = (lambda / mu2) / (1.0 - lambda / mu2);   // = 1.5
        SolverSSA solver = serial(tandem(lambda, mu1, mu2));
        NetworkAvgTable t = solver.getAvgTable();
        assertRanSerial(solver);
        List<Double> q = t.getQLen();   // Source(0), P1(1), P2(2)
        List<Double> tp = t.getTput();
        assertClose(q1, q.get(1), "tandem P1 mean tokens");
        assertClose(q2, q.get(2), "tandem P2 mean tokens");
        assertClose(lambda, tp.get(1), "tandem P1 throughput");
        assertClose(lambda, tp.get(2), "tandem P2 throughput");
    }
}
