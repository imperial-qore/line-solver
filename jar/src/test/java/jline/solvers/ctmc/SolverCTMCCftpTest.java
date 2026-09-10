package jline.solvers.ctmc;

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
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Regression tests for the perfect-sampling method of SolverCTMC,
 * SolverCTMC(model, "cftp"), which draws iid states from the exact stationary
 * distribution by monotone Coupling From The Past instead of enumerating the
 * state space.
 *
 * The reference is SolverCTMC with its default enumerating method on the same
 * model. Assertions cover the exact invariants that hold at any sample size
 * (population conservation, flow balance, Little's law, C = N/X, U in [0,1]),
 * the accuracy against the exact CTMC bounded by the Monte Carlo error of the
 * configured sample size, and the model-class gate, which must refuse every
 * model outside the closed single-class product form rather than approximate it.
 *
 * The MATLAB twin is line-test.git/test/testsCTMC/test_ctmc_cftp.m and the Python twin is
 * python/tests/test_ctmc_cftp.py; the three assert the same bounds.
 */
public class SolverCTMCCftpTest {

    private static final int SAMPLES = 20000;

    private static Network cqnModel(SchedStrategy schedQ2, int nservers) {
        Network model = new Network("CQN");
        Delay delay = new Delay(model, "Think");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Q2", schedQ2);
        q2.setNumberOfServers(nservers);
        ClosedClass cl = new ClosedClass(model, "C1", 8, delay, 0);
        delay.setService(cl, new Exp(1 / 0.5));
        q1.setService(cl, new Exp(1 / 1.0));
        q2.setService(cl, new Exp(1 / 0.6));
        model.link(model.serialRouting(delay, q1, q2));
        return model;
    }

    private static SolverCTMC cftpSolver(Network model, String method, int seed) {
        SolverOptions options = new SolverOptions(SolverType.CTMC);
        options.method = method;
        options.samples = SAMPLES;
        options.seed = seed;
        return new SolverCTMC(model, options);
    }

    private static double maxAbsDiff(Matrix a, Matrix b) {
        double worst = 0.0;
        for (int i = 0; i < a.getNumRows(); i++) {
            for (int j = 0; j < a.getNumCols(); j++) {
                worst = Math.max(worst, Math.abs(a.get(i, j) - b.get(i, j)));
            }
        }
        return worst;
    }

    @Test
    public void invariantsHoldAtAnySampleSize() throws Exception {
        SolverCTMC solver = cftpSolver(cqnModel(SchedStrategy.FCFS, 2), "cftp", 7);
        Matrix QN = solver.getAvgQLen();
        Matrix UN = solver.getAvgUtil();
        Matrix RN = solver.getAvgRespT();
        Matrix TN = solver.getAvgTput();
        double XN = solver.getAvgSysTput().get(0);
        double CN = solver.getAvgSysRespT().get(0);
        int N = 8;

        // every sampled state holds the whole closed population
        assertEquals(N, QN.elementSum(), 1e-9);
        for (int i = 0; i < TN.getNumRows(); i++) {
            // a single-class series network has V = 1 everywhere, so throughput is common
            assertEquals(TN.get(0, 0), TN.get(i, 0), 1e-12);
            // Little's law per station
            assertEquals(QN.get(i, 0) / TN.get(i, 0), RN.get(i, 0), 1e-12);
            assertTrue(UN.get(i, 0) >= 0.0, "utilization below zero at station " + i);
            if (i > 0) {
                assertTrue(UN.get(i, 0) <= 1.0, "utilization above one at station " + i);
            }
        }
        assertEquals(N / XN, CN, 1e-12);
    }

    @Test
    public void accuracyAgainstExactCtmc() throws Exception {
        SolverCTMC exact = new SolverCTMC(cqnModel(SchedStrategy.FCFS, 2));
        Matrix Qe = exact.getAvgQLen();
        Matrix Ue = exact.getAvgUtil();
        double Xe = exact.getAvgSysTput().get(0);

        SolverCTMC sampled = cftpSolver(cqnModel(SchedStrategy.FCFS, 2), "cftp", 7);
        Matrix Qc = sampled.getAvgQLen();
        Matrix Uc = sampled.getAvgUtil();
        double Xc = sampled.getAvgSysTput().get(0);

        // 20000 iid draws put the standard error of a queue length below 1e-2 of
        // the population; the bounds below were measured at implementation time.
        assertTrue(maxAbsDiff(Qc, Qe) < 0.1, "queue length error " + maxAbsDiff(Qc, Qe));
        assertTrue(maxAbsDiff(Uc, Ue) < 0.05, "utilization error " + maxAbsDiff(Uc, Ue));
        assertTrue(Math.abs(Xc - Xe) / Xe < 0.05, "throughput error " + Math.abs(Xc - Xe) / Xe);
    }

    @Test
    public void approximateSamplerAgreesWithExactSampler() throws Exception {
        Matrix Qe = new SolverCTMC(cqnModel(SchedStrategy.PS, 1)).getAvgQLen();
        Matrix Qa = cftpSolver(cqnModel(SchedStrategy.PS, 1), "cftp.approx", 7).getAvgQLen();
        assertTrue(maxAbsDiff(Qa, Qe) < 0.15, "queue length error " + maxAbsDiff(Qa, Qe));
        assertEquals(8, Qa.elementSum(), 1e-9);
    }

    @Test
    public void multiserverBalanceFunctionIsExact() throws Exception {
        // min(n,c) at the multiserver station and n! at the delay: a wrong server
        // count shows up as a biased queue length
        Matrix Qe = new SolverCTMC(cqnModel(SchedStrategy.FCFS, 3)).getAvgQLen();
        Matrix Qc = cftpSolver(cqnModel(SchedStrategy.FCFS, 3), "cftp", 11).getAvgQLen();
        assertTrue(maxAbsDiff(Qc, Qe) < 0.1, "queue length error " + maxAbsDiff(Qc, Qe));
    }

    @Test
    public void resultCarriesTheDrawnSamples() throws Exception {
        // the sampled states are the only representation of the stationary
        // distribution this method produces, since no state space is enumerated
        SolverCTMC solver = cftpSolver(cqnModel(SchedStrategy.FCFS, 2), "cftp", 7);
        solver.getAvgQLen();
        CTMCResult res = (CTMCResult) solver.result;
        assertEquals(SAMPLES, res.cftpSamples.getNumRows());
        assertEquals(3, res.cftpSamples.getNumCols());
        assertEquals(SAMPLES, res.cftpHorizon.getNumRows());
        for (int s = 0; s < SAMPLES; s++) {
            double total = 0.0;
            for (int i = 0; i < 3; i++) {
                total += res.cftpSamples.get(s, i);
            }
            assertEquals(8, total, 1e-12, "sample " + s + " does not hold the population");
            assertTrue(res.cftpHorizon.get(s, 0) >= 1, "non-positive horizon at sample " + s);
        }
        // the empirical measure over the distinct states is a distribution
        assertEquals(1.0, res.pi.elementSum(), 1e-9);
        assertEquals(res.pi.getNumRows(), res.spaceAggr.getNumRows());
    }

    @Test
    public void rejectsOpenModel() {
        Network model = new Network("OQN");
        Source source = new Source(model, "Source");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass cl = new OpenClass(model, "C1");
        source.setArrival(cl, new Exp(1));
        q1.setService(cl, new Exp(2));
        model.link(model.serialRouting(source, q1, sink));
        RuntimeException err = assertThrows(RuntimeException.class,
                () -> cftpSolver(model, "cftp", 7).getAvgQLen());
        assertTrue(err.getMessage().contains("cftp method supports"), err.getMessage());
    }

    @Test
    public void rejectsMulticlassModel() {
        Network model = new Network("CQN2");
        Delay delay = new Delay(model, "Think");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "C1", 4, delay, 0);
        ClosedClass c2 = new ClosedClass(model, "C2", 3, delay, 0);
        delay.setService(c1, new Exp(1));
        delay.setService(c2, new Exp(1));
        q1.setService(c1, new Exp(2));
        q1.setService(c2, new Exp(3));
        q2.setService(c1, new Exp(2));
        q2.setService(c2, new Exp(3));
        model.link(model.serialRouting(delay, q1, q2));
        RuntimeException err = assertThrows(RuntimeException.class,
                () -> cftpSolver(model, "cftp", 7).getAvgQLen());
        assertTrue(err.getMessage().contains("single-class models only"), err.getMessage());
    }

    @Test
    public void rejectsNonExponentialService() {
        Network model = new Network("CQNph");
        Delay delay = new Delay(model, "Think");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        ClosedClass cl = new ClosedClass(model, "C1", 4, delay, 0);
        delay.setService(cl, new Exp(2));
        q1.setService(cl, Erlang.fitMeanAndOrder(1, 2));
        q2.setService(cl, new Exp(1));
        model.link(model.serialRouting(delay, q1, q2));
        RuntimeException err = assertThrows(RuntimeException.class,
                () -> cftpSolver(model, "cftp", 7).getAvgQLen());
        assertTrue(err.getMessage().contains("exponential service times"), err.getMessage());
    }
}
