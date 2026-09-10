package jline.solvers;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.io.Ret.DistributionResult;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.fluid.SolverFluid;
import jline.solvers.ldes.LDESOptions;
import jline.solvers.ldes.SolverLDES;
import jline.solvers.mva.SolverMVA;
import jline.solvers.ssa.SolverSSA;
import jline.util.matrix.Matrix;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

/**
 * The get*Cdf* family serves the same support matrix as the MATLAB reference:
 * LDES answers getCdfRespT with the MEASURED ecdf (not the inherited
 * exponential fit), SSA refuses by name, Fluid's handle overload reaches the
 * passage-time law (it used to fall through to the base fallback), MVA
 * inherits the base exponential law, and CTMC's state-set first passage
 * exists, as it does in MATLAB and python.
 */
class CdfRespTFamilyParityTest {

    @BeforeAll
    public static void setUp() {
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    private static Network mm1() {
        Network model = new Network("mm1cdf");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass k = new OpenClass(model, "Class1");
        source.setArrival(k, new Exp(1.0));
        queue.setService(k, new Exp(2.0));
        RoutingMatrix routing = model.initRoutingMatrix();
        routing.set(k, k, source, queue, 1.0);
        routing.set(k, k, queue, sink, 1.0);
        model.link(routing);
        return model;
    }

    private static Network closedCycle() {
        Network model = new Network("closedcdf");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        ClosedClass k = new ClosedClass(model, "Class1", 3, delay);
        delay.setService(k, new Exp(1.0));
        queue.setService(k, new Exp(2.0));
        RoutingMatrix routing = model.initRoutingMatrix();
        routing.set(k, k, delay, queue, 1.0);
        routing.set(k, k, queue, delay, 1.0);
        model.link(routing);
        return model;
    }

    /** Mean of a tabulated [F(t), t] law, integrating t dF. */
    private static double cdfMean(Matrix cdf) {
        double mean = 0.0;
        double prevF = 0.0;
        for (int i = 0; i < cdf.getNumRows(); i++) {
            mean += cdf.get(i, 1) * (cdf.get(i, 0) - prevF);
            prevF = cdf.get(i, 0);
        }
        return mean;
    }

    @Test
    void ldesServesTheMeasuredEcdfUnderGetCdfRespT() {
        LDESOptions o = new LDESOptions();
        o.samples = 20000;
        o.seed = 23000;
        o.verbose = VerboseLevel.SILENT;
        SolverLDES s = new SolverLDES(mm1(), o);
        DistributionResult rd = s.getCdfRespT();
        // Station 1 is the queue; the ecdf must reproduce the reported mean,
        // which is the check that both read the same samples
        Matrix cdf = rd.getCdf(1, 0);
        assertNotNull(cdf);
        assertTrue(cdf.getNumRows() > 100, "a 20000-sample run tabulates a real ecdf, not a 100-point fit");
        assertEquals(1.0, cdf.get(cdf.getNumRows() - 1, 0), 1e-12);
        double rn = s.getAvgTable().getRespT().get(1);
        assertEquals(rn, cdfMean(cdf), 1e-8 * Math.max(1.0, rn));
        // The transient name serves the same curve
        Matrix cdfTran = s.getTranCdfRespT().getCdf(1, 0);
        assertEquals(cdf.getNumRows(), cdfTran.getNumRows());
        assertEquals(cdf.get(0, 1), cdfTran.get(0, 1), 0.0);
    }

    @Test
    void ssaRefusesTheResponseTimeCdfByName() {
        SolverSSA s = new SolverSSA(mm1());
        RuntimeException ex = assertThrows(RuntimeException.class, s::getCdfRespT);
        assertTrue(ex.getMessage().contains("does not record per-job response times"));
    }

    @Test
    void fluidHandleOverloadServesThePassageLawNotTheFallback() {
        SolverFluid s = new SolverFluid(closedCycle());
        DistributionResult noArg = s.getCdfRespT();
        DistributionResult viaHandle = s.getCdfRespT(s.getAvgRespTHandles());
        Matrix a = noArg.getCdf(1, 0);
        Matrix b = viaHandle.getCdf(1, 0);
        assertNotNull(a);
        assertNotNull(b);
        assertEquals(a.getNumRows(), b.getNumRows(),
                "the handle form must reach the same passage-time law as the no-arg form");
        // The base fallback tabulates exactly 100 quantile points; the fluid
        // passage law does not
        assertNotEquals(100, a.getNumRows());
    }

    @Test
    void mvaInheritsTheBaseExponentialFallback() {
        SolverMVA s = new SolverMVA(closedCycle());
        DistributionResult rd = s.getCdfRespT();
        Matrix cdf = rd.getCdf(1, 0);
        assertNotNull(cdf);
        assertEquals(100, cdf.getNumRows());
        // The law is exponential with the reported mean: t = -ln(1-F) * RN
        double rn = s.getAvgTable().getRespT().get(1);
        for (int i = 0; i < cdf.getNumRows(); i += 25) {
            double F = cdf.get(i, 0);
            assertEquals(-Math.log(1.0 - F) * rn, cdf.get(i, 1), 1e-8 * Math.max(1.0, rn));
        }
    }

    @Test
    void ctmcFirstPassageIsAProperCdf() {
        SolverCTMC s = new SolverCTMC(closedCycle());
        Matrix space = s.getStateSpace().stateSpace;
        assertTrue(space.getNumRows() > 1);
        // Into the last state, from the conditional stationary law on the rest
        Matrix B = new Matrix(1, 1);
        B.set(0, 0, space.getNumRows());
        SolverCTMC.FirstPassageResult out = s.getCdfFirstPassT(null, B);
        assertEquals(1000, out.tset.length);
        assertEquals(out.RD.getNumRows(), out.tset.length);
        double prev = 0.0;
        for (int i = 0; i < out.RD.getNumRows(); i++) {
            double F = out.RD.get(i, 0);
            assertTrue(F >= prev - 1e-12, "F(t) must be nondecreasing");
            assertTrue(F >= 0.0 && F <= 1.0);
            prev = F;
        }
        assertTrue(out.RD.get(out.RD.getNumRows() - 1, 0) > 0.99,
                "the passage into a positive-rate state completes within the horizon");
    }
}
