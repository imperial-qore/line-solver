package jline.solvers.mam;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Det;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.Gamma;
import jline.lang.processes.MMPP2;
import jline.lang.processes.Uniform;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.util.matrix.Matrix;

/**
 * SolverMAM must answer a station with NON-PHASE-TYPE, class-dependent service
 * through MMAP[K]/G[K]/1 rather than through the phase-type fit.
 *
 * The generic path reads the service law out of sn.proc, which for a Uniform is
 * a 20-phase PH fit: it matches the mean and, above SCV 1, nothing else. He's
 * transform analysis takes the original law. THE ORACLE IS THE MATLAB REFERENCE,
 * which reaches the same numbers through sn.lst, and JMT at 8e6 samples agrees
 * with both to 5e-4 per class.
 */
public class SolverMamGeneralServiceTest {

    @Test
    public void generalClassDependentServiceUsesTheTransformAnalysis() {
        Network model = new Network("gk1");
        Source src = new Source(model, "S");
        Queue q = new Queue(model, "Q", SchedStrategy.FCFS);
        Sink snk = new Sink(model, "K");
        OpenClass c1 = new OpenClass(model, "C1");
        OpenClass c2 = new OpenClass(model, "C2");
        src.setArrival(c1, new Exp(0.35));
        src.setArrival(c2, new MMPP2(0.9, 0.2, 0.25, 0.35));
        q.setService(c1, new Det(0.5));
        q.setService(c2, new Uniform(0.2, 1.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c1, c1, src, q, 1.0);
        P.set(c1, c1, q, snk, 1.0);
        P.set(c2, c2, src, q, 1.0);
        P.set(c2, c2, q, snk, 1.0);
        model.link(P);

        Matrix R = new SolverMAM(model).getAvgRespT();
        assertEquals(0.936113, R.get(1, 0), 1e-5);
        assertEquals(1.121427, R.get(1, 1), 1e-5);
    }

    /**
     * A service law with NO Det beside it, so nothing but the declared tag can
     * select the transform analysis.
     *
     * SnNonmarkovToPh runs before the dispatch and retags the station APH, and
     * DET is the only tag it preserves. A gate reading sn.procid therefore held
     * for the case above, on the strength of its Det alone, while a Gamma on its
     * own fell to the phase-type fit: at rho 0.5 with SCV 3 that read
     * 1.0785538590 against the Pollaczek-Khinchine value 1.5. P-K needs only the
     * first two moments of an M/G/1 service, so it is an oracle independent of
     * every path in the solver.
     */
    @Test
    public void declaredLawSelectsTheTransformAnalysisWithoutADetBesideIt() {
        Network model = new Network("gamma");
        Source src = new Source(model, "S");
        Queue q = new Queue(model, "Q", SchedStrategy.FCFS);
        Sink snk = new Sink(model, "K");
        OpenClass cls = new OpenClass(model, "C");
        src.setArrival(cls, new Exp(0.5));
        q.setService(cls, Gamma.fitMeanAndSCV(1.0, 3.0));
        model.link(model.serialRouting(src, q, snk));
        Matrix Q = new SolverMAM(model).getAvgQLen();
        assertEquals(1.5, Q.get(1, 0), 1e-6);
    }

    @Test
    public void phaseTypeOnlyServiceIsLeftOnTheExistingPath() {
        Network model = new Network("ph");
        Source src = new Source(model, "S");
        Queue q = new Queue(model, "Q", SchedStrategy.FCFS);
        Sink snk = new Sink(model, "K");
        OpenClass cls = new OpenClass(model, "C");
        src.setArrival(cls, new Exp(0.5));
        q.setService(cls, Erlang.fitMeanAndSCV(1.0, 0.5));
        model.link(model.serialRouting(src, q, snk));
        Matrix Q = new SolverMAM(model).getAvgQLen();
        assertEquals(0.875, Q.get(1, 0), 1e-5);
    }
}
