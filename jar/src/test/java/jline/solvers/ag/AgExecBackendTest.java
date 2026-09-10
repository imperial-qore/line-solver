package jline.solvers.ag;

import java.util.ArrayList;
import java.util.List;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.mam.SolverMAM;
import jline.util.matrix.Matrix;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The execution backends of SolverAG, and the solver boundary the RCAT move drew.
 *
 * <p>THE IDENTITY ASSERTIONS ARE EXACT ON PURPOSE. Agent k's generator is
 * Q_k(x) = L_k + sum_c x_c Pb_c, so an agent reads the rest of the model only
 * through the scalar reversed rates and writes only its own slot; the sweep is
 * Jacobi, so the agent order is immaterial. 'parallel' therefore has to reproduce
 * 'serial' BIT FOR BIT, and an approximate assertion here would pass on a
 * backend that had quietly started racing.</p>
 */
public class AgExecBackendTest {

    /** An M/M/1 tandem: three agents, and a product form RCAT solves exactly. */
    private static Network tandem() {
        Network model = new Network("Tandem");
        Source source = new Source(model, "Source");
        Queue q1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass cls = new OpenClass(model, "C");
        source.setArrival(cls, new Exp(1.0));
        q1.setService(cls, new Exp(2.0));
        q2.setService(cls, new Exp(3.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(cls, cls, source, q1, 1.0);
        P.set(cls, cls, q1, q2, 1.0);
        P.set(cls, cls, q2, sink, 1.0);
        model.link(P);
        return model;
    }

    private static AGResult solve(String exec, String method, int nworkers) {
        AGOptions opts = SolverAG.defaultOptions();
        opts.method = method;
        opts.exec = exec;
        opts.nworkers = nworkers;
        SolverAG solver = new SolverAG(tandem(), opts);
        solver.runAnalyzer();
        return solver.getAGResult();
    }

    private static void assertBitIdentical(String what, Matrix a, Matrix b) {
        assertEquals(a.getNumRows(), b.getNumRows(), what + ": row count");
        assertEquals(a.getNumCols(), b.getNumCols(), what + ": column count");
        for (int i = 0; i < a.getNumRows(); i++) {
            for (int j = 0; j < a.getNumCols(); j++) {
                assertEquals(Double.doubleToLongBits(a.get(i, j)),
                        Double.doubleToLongBits(b.get(i, j)),
                        what + " differs at (" + i + "," + j + "): "
                                + a.get(i, j) + " vs " + b.get(i, j));
            }
        }
    }

    @Test
    public void testTandemIsTwoIsolatedMM1Queues() {
        // Burke: the departure stream of an M/M/1 in equilibrium is Poisson at
        // the arrival rate, so QLen is rho/(1-rho) at each queue.
        AGResult r = solve(AGOptions.EXEC_SERIAL, "inap", 0);
        assertEquals(1.0, r.QN.get(1, 0), 1e-6);   // rho = 1/2
        assertEquals(0.5, r.QN.get(2, 0), 1e-6);   // rho = 1/3
    }

    @Test
    public void testThreadsBackendIsBitIdenticalToSerial() {
        AGResult ref = solve(AGOptions.EXEC_SERIAL, "inap", 0);
        AGResult got = solve(AGOptions.EXEC_PARALLEL, "inap", 0);
        assertBitIdentical("QN", ref.QN, got.QN);
        assertBitIdentical("UN", ref.UN, got.UN);
        assertBitIdentical("RN", ref.RN, got.RN);
        assertBitIdentical("TN", ref.TN, got.TN);
        assertEquals(ref.iter, got.iter, "sweep count");
    }

    @Test
    public void testThreadsBackendIsBitIdenticalAtEveryPoolSize() {
        // The partition is over agents and every agent is solved wherever it
        // lands, so the pool size cannot enter the answer.
        AGResult ref = solve(AGOptions.EXEC_SERIAL, "inap", 0);
        int[] sizes = new int[]{1, 2, 8};
        for (int i = 0; i < sizes.length; i++) {
            AGResult got = solve(AGOptions.EXEC_PARALLEL, "inap", sizes[i]);
            assertBitIdentical("QN at nworkers=" + sizes[i], ref.QN, got.QN);
        }
    }

    @Test
    public void testThreadsBackendIsBitIdenticalUnderInapplus() {
        AGResult ref = solve(AGOptions.EXEC_SERIAL, "inapplus", 0);
        AGResult got = solve(AGOptions.EXEC_PARALLEL, "inapplus", 0);
        assertBitIdentical("QN", ref.QN, got.QN);
    }

    @Test
    public void testThreadsBackendIsBitIdenticalUnderInapinf() {
        // The matrix-geometric path has its own per-agent solve; it must fan out
        // exactly as the finite one does.
        AGResult ref = solve(AGOptions.EXEC_SERIAL, "inapinf", 0);
        AGResult got = solve(AGOptions.EXEC_PARALLEL, "inapinf", 0);
        assertBitIdentical("QN", ref.QN, got.QN);
    }

    @Test
    public void testUnknownBackendIsRefusedByName() {
        AGOptions opts = SolverAG.defaultOptions();
        opts.exec = "gpu";
        assertThrows(RuntimeException.class, () -> AgExec.create(opts));
    }

    @Test
    public void testClusterWithoutEndpointsIsRefusedByName() {
        // A cluster run with nowhere to send the agents is a configuration
        // error, not a silent local run: the caller asked for distribution and
        // would otherwise be told nothing when they did not get it.
        AGOptions opts = SolverAG.defaultOptions();
        opts.exec = AGOptions.EXEC_CLUSTER;
        assertThrows(RuntimeException.class, () -> AgExec.create(opts));
    }

    @Test
    public void testClusterRefusesInapinfRatherThanSubstitutingAMethod() {
        // The ag-worker carries the FINITE agent solve; 'inapinf' replaces it
        // with the matrix-geometric tail of an open agent, which the worker does
        // not have. Answering with the finite solve would change the method.
        AGOptions opts = SolverAG.defaultOptions();
        opts.method = "inapinf";
        opts.exec = AGOptions.EXEC_CLUSTER;
        opts.endpoints.add("127.0.0.1:9");
        SolverAG solver = new SolverAG(tandem(), opts);
        RuntimeException e = assertThrows(RuntimeException.class, solver::runAnalyzer);
        assertTrue(e.getMessage().contains("inapinf"),
                "the refusal must name the method it cannot carry: " + e.getMessage());
    }

    @Test
    public void testClusterDegradesToLocalWhenNoWorkerAnswers() {
        // A lost worker costs wall clock and nothing else: any agent can be
        // solved anywhere given x, so an unreachable endpoint must still produce
        // the run's answer rather than an error or a wrong number.
        AGResult ref = solve(AGOptions.EXEC_SERIAL, "inap", 0);

        AGOptions opts = SolverAG.defaultOptions();
        opts.method = "inap";
        opts.exec = AGOptions.EXEC_CLUSTER;
        opts.endpoints.add("127.0.0.1:9");   // discard port: refuses every connection
        opts.workerTimeout = 2.0;
        SolverAG solver = new SolverAG(tandem(), opts);
        solver.runAnalyzer();
        assertBitIdentical("QN", ref.QN, solver.getAGResult().QN);
    }

    @Test
    public void testSolverMamRedirectsTheMovedRcatMethods() {
        // A caller carrying an old options.method must be told where the method
        // went, not that it is unknown.
        String[] moved = new String[]{"inap", "inapplus", "inapinf", "exact"};
        for (int i = 0; i < moved.length; i++) {
            SolverMAM mam = new SolverMAM(tandem(), moved[i]);
            String reason = mam.supportsModelMethod(moved[i]);
            assertFalse(reason.isEmpty(), moved[i] + " must be refused by SolverMAM");
            assertTrue(reason.contains("SolverAG"),
                    "the refusal must name SolverAG: " + reason);
        }
    }

    @Test
    public void testTheGNetworkFeaturesBelongToAgAlone() {
        // The RCAT builder is the only code in LINE that reads sn.issignal, so
        // no MAM method may declare the signal names: declaring them is what let
        // a G-network reach a decomposition that ignores the marking.
        String[] names = new String[]{"OpenSignal", "ClosedSignal",
                "SignalType_NEGATIVE", "SignalType_CATASTROPHE", "SignalBatchRemoval"};
        for (int i = 0; i < names.length; i++) {
            assertTrue(SolverAG.getFeatureSet().inspectFeature(names[i]),
                    "SolverAG must declare " + names[i]);
            assertFalse(SolverMAM.getFeatureSet().inspectFeature(names[i]),
                    "SolverMAM must not declare " + names[i]);
        }
    }

    @Test
    public void testTheTwoSolversAdvertiseDisjointMethodSets() {
        List<String> ag = new SolverAG(tandem()).listValidMethods();
        List<String> mam = new SolverMAM(tandem()).listValidMethods();
        List<String> shared = new ArrayList<String>(ag);
        shared.retainAll(mam);
        assertEquals(1, shared.size(), "only 'default' may be shared, got " + shared);
        assertEquals("default", shared.get(0));
        assertTrue(ag.contains("inap") && ag.contains("inapplus") && ag.contains("inapinf"));
        assertFalse(mam.contains("inap"));
    }

    @Test
    public void testExactReportsTheMethodThatActuallyRan() {
        // 'exact' is the vestigial AutoCAT alias and falls back to inap. It is
        // classified globally as an EXACT method, so leaving the name on the
        // result would banner an iterative approximation as exact.
        assertEquals("inap", solve(AGOptions.EXEC_SERIAL, "exact", 0).method);
        assertEquals("inap", solve(AGOptions.EXEC_SERIAL, "default", 0).method);
    }

    /**
     * "para" is an accepted alias of "parallel"; "threads" is refused BY NAME.
     *
     * The alias pair is SolverSSA's, so one spelling convention covers both
     * solvers. "threads" was this backend's name until 2026-08-19 and must not
     * fall through to the generic unknown-backend message, or a caller with an old
     * script is told a backend that still exists does not.
     */
    @Test
    public void testParaAliasAndThreadsRemoved() {
        assertTrue(AgExec.isParallel(AGOptions.EXEC_PARALLEL));
        assertTrue(AgExec.isParallel(AGOptions.EXEC_PARA));
        assertFalse(AgExec.isParallel("threads"));
        assertFalse(AgExec.isParallel(AGOptions.EXEC_SERIAL));
    }
}
