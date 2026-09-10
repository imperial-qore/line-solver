package jline.examples.advanced;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.examples.java.advanced.CyclicPollingModel;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.PollingType;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.ldes.LDESOptions;
import jline.solvers.ldes.SolverLDES;
import jline.solvers.mva.SolverMVA;
import jline.util.Maths;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Validation tests for the decrementing (semiexhaustive) polling discipline.
 * <p>
 * The symmetric decrementing polling system has the exact closed-form mean waiting
 * time of Pittel (1973) and Takagi (1984) (ACM Comput. Surv. 20(1), 1988, eq. 28).
 * These tests check that:
 * <ol>
 *   <li>the MVA analyzer reproduces the closed-form value, and</li>
 *   <li>the LDES simulator agrees with it, and</li>
 *   <li>the mean waits obey the symmetric ordering
 *       E[W]_exhaustive &le; E[W]_decrementing &le; E[W]_1-limited (Takagi eq. 12b).</li>
 * </ol>
 */
public class PollingDecrementingTest {

    private static final int QUEUE_STATION = 1; // stations: 0=Source, 1=Queue

    @BeforeAll
    public static void setUp() {
        Maths.setRandomNumbersMatlab(true);
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    /** Symmetric decrementing polling model with configurable polling type. */
    private static Network symmetricPolling(PollingType type) {
        Network model = new Network("polling-sym");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.POLLING);
        Sink sink = new Sink(model, "mySink");

        OpenClass c1 = new OpenClass(model, "myClass1");
        source.setArrival(c1, new Exp(0.4));
        queue.setService(c1, new Exp(1.0));
        OpenClass c2 = new OpenClass(model, "myClass2");
        source.setArrival(c2, new Exp(0.4));
        queue.setService(c2, new Exp(1.0));

        if (type == PollingType.KLIMITED) {
            queue.setPollingType(PollingType.KLIMITED, 1);
        } else {
            queue.setPollingType(type);
        }
        queue.setSwitchover(c1, new Exp(10.0));
        queue.setSwitchover(c2, new Exp(10.0));

        model.addLink(source, queue);
        model.addLink(queue, sink);
        source.setProbRouting(c1, queue, 1.0);
        queue.setProbRouting(c1, sink, 1.0);
        source.setProbRouting(c2, queue, 1.0);
        queue.setProbRouting(c2, sink, 1.0);
        return model;
    }

    /**
     * Closed-form symmetric decrementing mean waiting time (Takagi eq. 28) for the
     * model above: N=2, lambda=0.4, service Exp(1) (b=1, b2=2), switchover Exp(10)
     * (r=0.1, delta2=0.01).
     */
    private static double analyticDecrementingW() {
        double N = 2, lam = 0.4, b2 = 2.0, r = 0.1, d2 = 0.01;
        double rho = N * lam * 1.0;
        double denom = 2.0 * (1.0 - rho - lam * r * (N - rho));
        return d2 / (2.0 * r)
                + (N * lam * b2 * (1.0 - lam * r) + (r + lam * d2) * (N - rho)) / denom;
    }

    /** Closed-form symmetric exhaustive mean waiting time (Takagi eq. 14). */
    private static double analyticExhaustiveW() {
        double N = 2, lam = 0.4, b2 = 2.0, r = 0.1, d2 = 0.01;
        double rho = N * lam * 1.0;
        return d2 / (2.0 * r) + (N * lam * b2 + r * (N - rho)) / (2.0 * (1.0 - rho));
    }

    private double mvaQueueRespT(Network model) {
        SolverOptions opts = new SolverOptions(SolverType.MVA);
        opts.verbose = VerboseLevel.SILENT;
        NetworkAvgTable t = new SolverMVA(model, opts).getAvgTable();
        // rows: (Source,c1),(Source,c2),(Queue,c1),(Queue,c2); take Queue/c1 RespT
        return t.getRespT().get(2);
    }

    private double ldesQueueRespT(Network model) {
        LDESOptions opts = new LDESOptions();
        opts.verbose = VerboseLevel.SILENT;
        opts.seed = 23000;
        // 500k samples, not 50k: the discipline-ratio assertion below compares a
        // RATIO of two mean waiting times at rho=0.8, whose Monte Carlo error at
        // 50k samples is 2-6% across seeds (measured), exceeding the 3% bound
        // even for a correct simulator. At 500k the error is 0.1-1.7%. The old
        // empty-buffer-skipping walk passed at 50k only because its shorter
        // cycles had lower variance and its bias cancelled in the ratio.
        opts.samples = 500000;
        NetworkAvgTable t = new SolverLDES(model, opts).getAvgTable();
        // Locate the myQueue / myClass1 row.
        for (int i = 0; i < t.getRespT().size(); i++) {
            if ("myQueue".equals(t.getStationNames().get(i))
                    && "myClass1".equals(t.getClassNames().get(i))) {
                return t.getRespT().get(i);
            }
        }
        fail("myQueue/myClass1 row not found in LDES table");
        return Double.NaN;
    }

    @Test
    public void testGatedBaselineMva() {
        // Regression guard: the full SolverMVA pipeline must solve an open polling
        // model end-to-end (previously broke because setPollingType/setSwitchover
        // cached an empty struct before links were added).
        Network model = CyclicPollingModel.polling_gated();
        jline.lang.NetworkStruct sn = model.getStruct();
        assertTrue(sn.nchains > 0, "open polling model must build chains");
        double respT = mvaQueueRespT(model);
        assertTrue(respT > 0);
    }

    @Test
    public void testMvaMatchesClosedForm() {
        Network model = symmetricPolling(PollingType.DECREMENTING);
        double mvaRespT = mvaQueueRespT(model);
        double mvaW = mvaRespT - 1.0; // subtract mean service (1/mu = 1)
        double analyticW = analyticDecrementingW();
        assertEquals(analyticW, mvaW, 1e-9,
                "MVA decrementing waiting time must equal Takagi eq.28 closed form");
    }

    /**
     * Validates the LDES decrementing discipline against the exact analysis.
     * <p>
     * LINE's LDES polling server uses a parking (non-roving) convention: it does not
     * incur a switchover when it would have to cycle back to a lone non-empty queue,
     * whereas Takagi's closed forms assume a continuously roving server. This produces
     * a small, discipline-independent negative offset in the LDES mean waits (the LDES
     * exhaustive value shows the same offset against Takagi eq. 14). The offset cancels
     * in the decrementing-to-exhaustive ratio, which isolates the discipline logic; that
     * ratio must match the analytic ratio, and the raw ordering must hold.
     */
    @Test
    public void testLdesDecrementingDiscipline() {
        double wExh = ldesQueueRespT(symmetricPolling(PollingType.EXHAUSTIVE)) - 1.0;
        double wDec = ldesQueueRespT(symmetricPolling(PollingType.DECREMENTING)) - 1.0;
        double wLim = ldesQueueRespT(symmetricPolling(PollingType.KLIMITED)) - 1.0;

        // (1) Discipline ordering, Takagi eq. 12b: exhaustive <= decrementing <= 1-limited.
        double slack = 0.02 * wDec;
        assertTrue(wExh <= wDec + slack, "E[W]_exhaustive must not exceed E[W]_decrementing");
        assertTrue(wDec <= wLim + slack, "E[W]_decrementing must not exceed E[W]_1-limited");

        // (2) Decrementing-to-exhaustive ratio matches the exact analysis (offset cancels).
        double ldesRatio = wDec / wExh;
        double analyticRatio = analyticDecrementingW() / analyticExhaustiveW();
        double ratioRelErr = Math.abs(ldesRatio - analyticRatio) / analyticRatio;
        assertTrue(ratioRelErr < 0.03,
                "LDES decrementing/exhaustive ratio must match analytic within 3%; relErr=" + ratioRelErr);
    }
}
