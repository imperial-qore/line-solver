package jline.solvers.ctmc;

import java.util.List;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Regression tests for the decision-diagram aggregation method of SolverCTMC,
 * SolverCTMC(model, "mdd"), which stores the reachable set in an MDD and solves
 * K coupled level-CTMCs instead of the |S|-state generator.
 *
 * The reference is SolverCTMC with its default enumerating method on the same
 * model. The method is EXACT on product-form networks (Miner-Ciardo-Donatelli
 * Sec. 5), so those cases assert equality to numerical tolerance rather than an
 * approximation bound; the non-product-form cases assert only the conservation
 * laws that hold regardless. The model-class gate must refuse every model
 * outside closed single-class rather than approximate it.
 *
 * The MATLAB twin is line-test.git and the python twin is the mdd cases in
 * python/tests; the three assert the same bounds.
 */
public class SolverCTMCMddTest {

    /** Product-form cases must agree with the exact solve to numerical noise. */
    private static final double EXACT_TOL = 1e-9;

    private static Network cyclic(int K, int N, double[] rates, SchedStrategy[] sched,
                                  Object[] procs, int[] servers) {
        Network model = new Network("mdd");
        jline.lang.nodes.Station[] nodes = new jline.lang.nodes.Station[K];
        for (int i = 0; i < K; i++) {
            if (sched != null && sched[i] == SchedStrategy.INF) {
                nodes[i] = new Delay(model, "S" + i);
            } else {
                Queue q = new Queue(model, "S" + i,
                        sched == null ? SchedStrategy.FCFS : sched[i]);
                if (servers != null) {
                    q.setNumberOfServers(servers[i]);
                }
                nodes[i] = q;
            }
        }
        ClosedClass job = new ClosedClass(model, "Jobs", N, nodes[0]);
        for (int i = 0; i < K; i++) {
            if (procs != null && procs[i] != null) {
                nodes[i].setService(job, (jline.lang.processes.Distribution) procs[i]);
            } else {
                nodes[i].setService(job, new Exp(rates[i]));
            }
        }
        model.link(model.serialRouting(nodes));
        return model;
    }

    private static double maxAbsDiff(Matrix a, Matrix b) {
        double d = 0;
        for (int i = 0; i < a.getNumRows(); i++) {
            for (int j = 0; j < a.getNumCols(); j++) {
                d = Math.max(d, Math.abs(a.get(i, j) - b.get(i, j)));
            }
        }
        return d;
    }

    /**
     * Select a method explicitly, and assert it was actually selected.
     *
     * <p>A method-selection test that never checks the method was selected
     * cannot fail for the right reason: before Solver.parseOptions learned to
     * read a positional method name, new SolverCTMC(model, "mdd") silently ran
     * the DEFAULT method, so comparing "mdd" against "exact" compared the
     * default with itself and every numerical assertion here passed vacuously.
     * The bare form works now; the explicit form plus the assertion is what
     * keeps that failure mode from returning.</p>
     */
    private static SolverCTMC solver(Network model, String method) {
        SolverOptions options = new SolverOptions(SolverType.CTMC);
        options.method = method;
        SolverCTMC s = new SolverCTMC(model, options);
        assertEquals(method, s.getOptions().method, "the requested method must be selected");
        return s;
    }

    private static void assertAgreesWithExact(Network model, String label) {
        SolverCTMC exact = solver(model, "exact");
        SolverCTMC mdd = solver(model, "mdd");
        assertEquals(0.0, maxAbsDiff(exact.getAvgQLen(), mdd.getAvgQLen()), EXACT_TOL,
                label + ": queue lengths must match the exact solve");
        assertEquals(0.0, maxAbsDiff(exact.getAvgUtil(), mdd.getAvgUtil()), EXACT_TOL,
                label + ": utilizations must match the exact solve");
        assertEquals(0.0, maxAbsDiff(exact.getAvgTput(), mdd.getAvgTput()), EXACT_TOL,
                label + ": throughputs must match the exact solve");
    }

    @Test
    public void mddIsListedAndExactIsAnAlias() {
        Network model = cyclic(3, 4, new double[]{1.0, 1.5, 2.0}, null, null, null);
        List<String> methods = new SolverCTMC(model).listValidMethods();
        assertTrue(methods.contains("mdd"), "mdd must be a valid SolverCTMC method");
        assertTrue(methods.contains("exact"), "exact must be a valid SolverCTMC method");
        // a positional method name must reach options.method, not be dropped
        assertEquals("mdd", new SolverCTMC(model, "mdd").getOptions().method,
                "a bare method string must select the method");
        // the alias must select the same code as the default, not merely be accepted
        SolverCTMC def = new SolverCTMC(model);
        SolverCTMC alias = solver(model, "exact");
        assertEquals(0.0, maxAbsDiff(def.getAvgQLen(), alias.getAvgQLen()), 0.0,
                "'exact' must be bit-identical to the default method");
    }

    @Test
    public void exponentialClosedNetworksAreExact() {
        // exponential service is product-form under any work-conserving
        // discipline, so the single approximation of the method is not active
        assertAgreesWithExact(cyclic(3, 4, new double[]{1.0, 1.5, 2.0}, null, null, null),
                "K=3 N=4 FCFS");
        assertAgreesWithExact(cyclic(4, 5, new double[]{1.0, 1.5, 2.0, 2.5}, null, null, null),
                "K=4 N=5 FCFS");
        assertAgreesWithExact(cyclic(5, 5, new double[]{1.0, 1.5, 2.0, 2.5, 3.0}, null, null, null),
                "K=5 N=5 FCFS");
    }

    @Test
    public void delayAndMultiserverStationsAreExact() {
        SchedStrategy[] sched = {SchedStrategy.INF, SchedStrategy.FCFS,
                SchedStrategy.FCFS, SchedStrategy.FCFS};
        assertAgreesWithExact(cyclic(4, 5, new double[]{1.0, 1.5, 2.0, 2.5}, sched, null, null),
                "K=4 N=5 with a delay station");
        assertAgreesWithExact(cyclic(4, 5, new double[]{1.0, 1.5, 2.0, 2.5}, null, null,
                new int[]{1, 2, 1, 1}), "K=4 N=5 with a 2-server station");
    }

    @Test
    public void processorSharingWithPhaseTypeIsInsensitive() {
        // BCMP type 2: the queue lengths depend on the service law only through
        // its mean, so Erlang and HyperExp must reproduce the exponential answer
        double[] rates = {1.0, 1.5, 2.0};
        SchedStrategy[] ps = {SchedStrategy.PS, SchedStrategy.PS, SchedStrategy.PS};
        Object[] erl = new Object[3];
        Object[] hyp = new Object[3];
        double[] scv = {4.0, 9.0, 2.0};   // a UNIFORM scv would only rescale every
        // rate by one factor, leaving the closed-network queue lengths invariant,
        // and so could not detect a wrong entry law
        for (int i = 0; i < 3; i++) {
            erl[i] = Erlang.fitMeanAndOrder(1.0 / rates[i], 2);
            hyp[i] = HyperExp.fitMeanAndSCV(1.0 / rates[i], scv[i]);
        }
        assertAgreesWithExact(cyclic(3, 4, rates, ps, erl, null), "PS Erlang-2");
        assertAgreesWithExact(cyclic(3, 4, rates, ps, hyp, null), "PS HyperExp");

        Matrix expQ = new SolverCTMC(cyclic(3, 4, rates, ps, null, null)).getAvgQLen();
        Matrix hypQ = solver(cyclic(3, 4, rates, ps, hyp, null), "mdd").getAvgQLen();
        assertEquals(0.0, maxAbsDiff(expQ, hypQ), 1e-9,
                "PS is insensitive: HyperExp must give the exponential queue lengths");
    }

    @Test
    public void nonProductFormStillConservesPopulation() {
        // FCFS with phase-type service is NOT product form, so only the
        // conservation laws are asserted, not agreement with the exact solve
        double[] rates = {1.0, 1.5, 2.0};
        Object[] erl = new Object[3];
        for (int i = 0; i < 3; i++) {
            erl[i] = Erlang.fitMeanAndOrder(1.0 / rates[i], 2);
        }
        Matrix Q = solver(cyclic(3, 4, rates, null, erl, null), "mdd").getAvgQLen();
        double total = 0;
        for (int i = 0; i < Q.getNumRows(); i++) {
            total += Q.get(i, 0);
        }
        assertEquals(4.0, total, 1e-6, "the closed population must be conserved");
    }

    @Test
    public void openModelsAreRejected() {
        Network model = new Network("open");
        Source source = new Source(model, "Source");
        Queue q = new Queue(model, "Q1", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oc = new OpenClass(model, "Class1");
        source.setArrival(oc, new Exp(0.5));
        q.setService(oc, new Exp(1.0));
        model.link(model.serialRouting(source, q, sink));
        assertThrows(RuntimeException.class,
                () -> solver(model, "mdd").getAvgQLen(),
                "an open stream makes the marking unbounded and must be refused");
    }

    @Test
    public void multiclassModelsAreRejected() {
        Network model = new Network("twoclass");
        Delay think = new Delay(model, "Think");
        Queue q = new Queue(model, "Q1", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "C1", 2, think);
        ClosedClass c2 = new ClosedClass(model, "C2", 2, think);
        think.setService(c1, new Exp(1.0));
        think.setService(c2, new Exp(1.0));
        q.setService(c1, new Exp(2.0));
        q.setService(c2, new Exp(3.0));
        model.link(model.serialRouting(think, q));
        assertThrows(RuntimeException.class,
                () -> solver(model, "mdd").getAvgQLen(),
                "the Kronecker descriptor would need one level per (station,class)");
    }

    @Test
    public void phaseTypeAtAPreemptiveStationIsRejected() {
        // the composite level names ONE in-service phase and restarts it at pie,
        // which is non-preemptive; under preemptive resume the suspended jobs'
        // phases would have to be stacked
        double[] rates = {1.0, 1.5, 2.0};
        SchedStrategy[] sched = {SchedStrategy.LCFSPR, SchedStrategy.FCFS, SchedStrategy.FCFS};
        Object[] erl = new Object[3];
        for (int i = 0; i < 3; i++) {
            erl[i] = Erlang.fitMeanAndOrder(1.0 / rates[i], 2);
        }
        assertThrows(RuntimeException.class,
                () -> solver(cyclic(3, 4, rates, sched, erl, null), "mdd").getAvgQLen(),
                "phase-type service at a preemptive-resume station must be refused");
    }
}
