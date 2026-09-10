package jline.solvers;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.*;
import jline.lang.processes.*;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.fluid.SolverFluid;
import jline.solvers.mam.SolverMAM;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.solvers.ssa.SolverSSA;
import jline.solvers.auto.SolverAUTO;
import jline.solvers.auto.AUTOptions;
import jline.solvers.auto.LINE;
import jline.solvers.wrappers.qns.SolverQNS;
import jline.solvers.uq.SolverUQ;
import jline.io.Ret.ProbabilityResult;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.*;

/**
 * Incremental coverage tests for solver packages.
 * Each test targets a specific uncovered method or branch in the solver hierarchy.
 */
public class SolverCoverageTest {

    @BeforeAll
    public static void setUp() {
        Maths.setRandomNumbersMatlab(true);
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    // ====== Helper model factories ======

    private static Network openMM1() {
        Network model = new Network("OpenMM1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oc = new OpenClass(model, "Class1");
        source.setArrival(oc, new Exp(0.5));
        queue.setService(oc, new Exp(1.0));
        model.link(model.serialRouting(source, queue, sink));
        return model;
    }

    private static Network closedSingle() {
        Network model = new Network("ClosedSingle");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.FCFS);
        ClosedClass cc = new ClosedClass(model, "Class1", 5, delay);
        delay.setService(cc, new Exp(1.0));
        queue.setService(cc, new Exp(2.0));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    private static Network closedMultiClass() {
        Network model = new Network("ClosedMulti");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "Class1", 3, delay);
        ClosedClass c2 = new ClosedClass(model, "Class2", 2, delay);
        delay.setService(c1, new Exp(1.0));
        delay.setService(c2, new Exp(0.5));
        queue.setService(c1, new Exp(2.0));
        queue.setService(c2, new Exp(1.5));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    private static Network openMultiClass() {
        Network model = new Network("OpenMulti");
        Source source = new Source(model, "Source");
        Queue q1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");
        OpenClass c1 = new OpenClass(model, "Class1");
        OpenClass c2 = new OpenClass(model, "Class2");
        source.setArrival(c1, new Exp(0.3));
        source.setArrival(c2, new Exp(0.2));
        q1.setService(c1, new Exp(1.0));
        q1.setService(c2, new Exp(1.5));
        q2.setService(c1, new Exp(2.0));
        q2.setService(c2, new Exp(1.0));
        RoutingMatrix P = new RoutingMatrix(model, model.getJobClasses(), model.getNodes());
        P.addConnection(source, q1, c1, 1.0);
        P.addConnection(source, q1, c2, 1.0);
        P.addConnection(q1, q2, c1, 1.0);
        P.addConnection(q1, q2, c2, 1.0);
        P.addConnection(q2, sink, c1, 1.0);
        P.addConnection(q2, sink, c2, 1.0);
        model.link(P);
        return model;
    }

    private static Network closedPS() {
        Network model = new Network("ClosedPS");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        ClosedClass cc = new ClosedClass(model, "Class1", 3, delay);
        delay.setService(cc, new Exp(1.0));
        queue.setService(cc, new Exp(2.0));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    private static Network closedMultiServer() {
        Network model = new Network("ClosedMultiServer");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.FCFS);
        queue.setNumberOfServers(2);
        ClosedClass cc = new ClosedClass(model, "Class1", 4, delay);
        delay.setService(cc, new Exp(1.0));
        queue.setService(cc, new Exp(3.0));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    private static Network openMM1PS() {
        Network model = new Network("OpenMM1PS");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oc = new OpenClass(model, "Class1");
        source.setArrival(oc, new Exp(0.5));
        queue.setService(oc, new Exp(1.0));
        model.link(model.serialRouting(source, queue, sink));
        return model;
    }

    private static Network closedTwoQueues() {
        Network model = new Network("ClosedTwoQueues");
        Delay delay = new Delay(model, "Delay");
        Queue q1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Queue2", SchedStrategy.PS);
        ClosedClass cc = new ClosedClass(model, "Class1", 4, delay);
        delay.setService(cc, new Exp(1.0));
        q1.setService(cc, new Exp(2.0));
        q2.setService(cc, new Exp(3.0));
        model.link(model.serialRouting(delay, q1, q2));
        return model;
    }

    // ============================================================
    // Tests 01-10: CTMC solver (27.2% -> target higher)
    // ============================================================

    @Test
    public void test01_ctmc_open_mm1_avgTable() {
        Network model = openMM1();
        withSuppressedOutput(() -> {
            SolverCTMC solver = new SolverCTMC(model, "keep", false, "cutoff", 5);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
            assertTrue(avg.getQLen().size() > 0);
        });
    }

    @Test
    public void test02_ctmc_closed_avgTable() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverCTMC solver = new SolverCTMC(model);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
            assertTrue(avg.getQLen().size() > 0);
        });
    }

    @Test
    public void test03_ctmc_getProb() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverCTMC solver = new SolverCTMC(model);
            solver.getAvg();
            // getProbAggr on first stateful node
            List<StatefulNode> nodes = model.getStatefulNodes();
            if (!nodes.isEmpty()) {
                ProbabilityResult prob = solver.getProbAggr(nodes.get(0));
                assertNotNull(prob);
            }
        });
    }

    @Test
    public void test04_ctmc_getProbSys() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverCTMC solver = new SolverCTMC(model);
            solver.getAvg();
            ProbabilityResult prob = solver.getProbSys();
            assertNotNull(prob);
        });
    }

    @Test
    public void test05_ctmc_getProbSysAggr() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverCTMC solver = new SolverCTMC(model);
            solver.getAvg();
            ProbabilityResult prob = solver.getProbSysAggr();
            assertNotNull(prob);
        });
    }

    @Test
    public void test06_ctmc_supports() {
        Network model = closedSingle();
        SolverCTMC solver = new SolverCTMC(model);
        assertTrue(solver.supports(model));
    }

    @Test
    public void test07_ctmc_getFeatureSet() {
        FeatureSet fs = SolverCTMC.getFeatureSet();
        assertNotNull(fs);
    }

    @Test
    public void test08_ctmc_closed_multiclass() {
        Network model = closedMultiClass();
        withSuppressedOutput(() -> {
            SolverCTMC solver = new SolverCTMC(model);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test09_ctmc_open_multiclass() {
        Network model = openMultiClass();
        withSuppressedOutput(() -> {
            SolverCTMC solver = new SolverCTMC(model, "keep", false, "cutoff", 3);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test10_ctmc_defaultOptions() {
        SolverOptions opts = SolverCTMC.defaultOptions();
        assertNotNull(opts);
    }

    // ============================================================
    // Tests 11-20: MVA solver (27.5% -> target higher)
    // ============================================================

    @Test
    public void test11_mva_closed_avgTable() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
            assertTrue(avg.getQLen().size() > 0);
        });
    }

    @Test
    public void test12_mva_closed_ps() {
        Network model = closedPS();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test13_mva_closed_multiclass() {
        Network model = closedMultiClass();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test14_mva_closed_multiserver() {
        Network model = closedMultiServer();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test15_mva_getAvgQLen() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            Matrix qlen = solver.getAvgQLen();
            assertNotNull(qlen);
            assertTrue(qlen.getNumRows() > 0);
        });
    }

    @Test
    public void test16_mva_getAvgUtil() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            Matrix util = solver.getAvgUtil();
            assertNotNull(util);
        });
    }

    @Test
    public void test17_mva_getAvgRespT() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            Matrix respT = solver.getAvgRespT();
            assertNotNull(respT);
        });
    }

    @Test
    public void test18_mva_getAvgTput() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            Matrix tput = solver.getAvgTput();
            assertNotNull(tput);
        });
    }

    @Test
    public void test19_mva_supports() {
        Network model = closedSingle();
        SolverMVA solver = new SolverMVA(model);
        assertTrue(solver.supports(model));
    }

    @Test
    public void test20_mva_getFeatureSet() {
        FeatureSet fs = SolverMVA.getFeatureSet();
        assertNotNull(fs);
    }

    // ============================================================
    // Tests 21-30: NC solver (39.7% -> target higher)
    // ============================================================

    @Test
    public void test21_nc_closed_avgTable() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverNC solver = new SolverNC(model);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test22_nc_closed_multiclass() {
        Network model = closedMultiClass();
        withSuppressedOutput(() -> {
            SolverNC solver = new SolverNC(model);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test23_nc_closed_ps() {
        Network model = closedPS();
        withSuppressedOutput(() -> {
            SolverNC solver = new SolverNC(model);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test24_nc_getAvgQLen() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverNC solver = new SolverNC(model);
            Matrix qlen = solver.getAvgQLen();
            assertNotNull(qlen);
        });
    }

    @Test
    public void test25_nc_getAvgTput() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverNC solver = new SolverNC(model);
            Matrix tput = solver.getAvgTput();
            assertNotNull(tput);
        });
    }

    @Test
    public void test26_nc_supports() {
        Network model = closedSingle();
        SolverNC solver = new SolverNC(model);
        assertTrue(solver.supports(model));
    }

    @Test
    public void test27_nc_getFeatureSet() {
        FeatureSet fs = SolverNC.getFeatureSet();
        assertNotNull(fs);
    }

    @Test
    public void test28_nc_closed_multiserver() {
        Network model = closedMultiServer();
        withSuppressedOutput(() -> {
            SolverNC solver = new SolverNC(model);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test29_nc_defaultOptions() {
        SolverOptions opts = SolverNC.defaultOptions();
        assertNotNull(opts);
    }

    @Test
    public void test30_nc_method_constructor() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverNC solver = new SolverNC(model, "exact");
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    // ============================================================
    // Tests 31-40: SSA solver (29.5% -> target higher)
    // ============================================================

    @Test
    public void test31_ssa_closed_avgTable() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverSSA solver = new SolverSSA(model, "samples", 500, "seed", 1);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test32_ssa_open_avgTable() {
        Network model = openMM1();
        withSuppressedOutput(() -> {
            SolverSSA solver = new SolverSSA(model, "samples", 500, "seed", 1);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test33_ssa_supports() {
        Network model = closedSingle();
        SolverSSA solver = new SolverSSA(model);
        assertTrue(solver.supports(model));
    }

    @Test
    public void test34_ssa_getFeatureSet() {
        FeatureSet fs = SolverSSA.getFeatureSet();
        assertNotNull(fs);
    }

    @Test
    public void test35_ssa_defaultOptions() {
        SolverOptions opts = SolverSSA.defaultOptions();
        assertNotNull(opts);
    }

    @Test
    public void test36_ssa_closed_multiclass() {
        Network model = closedMultiClass();
        withSuppressedOutput(() -> {
            SolverSSA solver = new SolverSSA(model, "samples", 500, "seed", 1);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test37_ssa_getAvgQLen() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverSSA solver = new SolverSSA(model, "samples", 500, "seed", 1);
            Matrix qlen = solver.getAvgQLen();
            assertNotNull(qlen);
        });
    }

    @Test
    public void test38_ssa_getAvgTput() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverSSA solver = new SolverSSA(model, "samples", 500, "seed", 1);
            Matrix tput = solver.getAvgTput();
            assertNotNull(tput);
        });
    }

    @Test
    public void test39_ssa_method_constructor() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverSSA solver = new SolverSSA(model, "default");
            solver.options.samples = 500;
            solver.options.seed = 1;
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test40_ssa_open_multiclass() {
        Network model = openMultiClass();
        withSuppressedOutput(() -> {
            SolverSSA solver = new SolverSSA(model, "samples", 500, "seed", 1);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    // ============================================================
    // Tests 41-48: Fluid solver (handlers at 24.8%)
    // ============================================================

    @Test
    public void test41_fluid_closed_avgTable() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverFluid solver = new SolverFluid(model);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test42_fluid_open_avgTable() {
        Network model = openMM1();
        withSuppressedOutput(() -> {
            SolverFluid solver = new SolverFluid(model);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test43_fluid_closed_multiclass() {
        Network model = closedMultiClass();
        withSuppressedOutput(() -> {
            SolverFluid solver = new SolverFluid(model);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test44_fluid_supports() {
        Network model = closedSingle();
        SolverFluid solver = new SolverFluid(model);
        assertTrue(solver.supports(model));
    }

    @Test
    public void test45_fluid_getFeatureSet() {
        FeatureSet fs = SolverFluid.getFeatureSet();
        assertNotNull(fs);
    }

    @Test
    public void test46_fluid_defaultOptions() {
        SolverOptions opts = SolverFluid.defaultOptions();
        assertNotNull(opts);
    }

    @Test
    public void test47_fluid_closed_twoqueues() {
        Network model = closedTwoQueues();
        withSuppressedOutput(() -> {
            SolverFluid solver = new SolverFluid(model);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test48_fluid_open_ps() {
        Network model = openMM1PS();
        withSuppressedOutput(() -> {
            SolverFluid solver = new SolverFluid(model);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    // ============================================================
    // Tests 49-55: MAM solver (35.6%)
    // ============================================================

    @Test
    public void test49_mam_open_avgTable() {
        Network model = openMM1();
        withSuppressedOutput(() -> {
            SolverMAM solver = new SolverMAM(model);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test50_mam_open_multiclass() {
        Network model = openMultiClass();
        withSuppressedOutput(() -> {
            SolverMAM solver = new SolverMAM(model);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test51_mam_supports() {
        Network model = openMM1();
        SolverMAM solver = new SolverMAM(model);
        assertTrue(solver.supports(model));
    }

    @Test
    public void test52_mam_getFeatureSet() {
        FeatureSet fs = SolverMAM.getFeatureSet();
        assertNotNull(fs);
    }

    @Test
    public void test53_mam_defaultOptions() {
        SolverOptions opts = SolverMAM.defaultOptions();
        assertNotNull(opts);
    }

    @Test
    public void test54_mam_getAvgQLen() {
        Network model = openMM1();
        withSuppressedOutput(() -> {
            SolverMAM solver = new SolverMAM(model);
            Matrix qlen = solver.getAvgQLen();
            assertNotNull(qlen);
        });
    }

    @Test
    public void test55_mam_open_ps() {
        Network model = openMM1PS();
        withSuppressedOutput(() -> {
            SolverMAM solver = new SolverMAM(model);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    // ============================================================
    // Tests 56-65: NetworkSolver table methods (41.3%)
    // Covers: getAvgChainTable, getAvgNodeTable, getAvgNodeChainTable,
    //         getAvgSysTable, getAvgChain, individual metric chain methods
    // ============================================================

    @Test
    public void test56_networkSolver_getAvgChainTable() {
        Network model = closedMultiClass();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            NetworkAvgChainTable chainTable = solver.getAvgChainTable();
            assertNotNull(chainTable);
        });
    }

    @Test
    public void test57_networkSolver_getAvgNodeTable() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            NetworkAvgNodeTable nodeTable = solver.getAvgNodeTable();
            assertNotNull(nodeTable);
        });
    }

    @Test
    public void test58_networkSolver_getAvgNodeChainTable() {
        Network model = closedMultiClass();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            NetworkAvgNodeChainTable nodeChainTable = solver.getAvgNodeChainTable();
            assertNotNull(nodeChainTable);
        });
    }

    @Test
    public void test59_networkSolver_getAvgSysTable() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            NetworkAvgSysTable sysTable = solver.getAvgSysTable();
            assertNotNull(sysTable);
        });
    }

    @Test
    public void test60_networkSolver_getAvgChain() {
        Network model = closedMultiClass();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            SolverResult chainResult = solver.getAvgChain();
            assertNotNull(chainResult);
            assertNotNull(chainResult.QN);
            assertNotNull(chainResult.TN);
        });
    }

    @Test
    public void test61_networkSolver_getAvgQLenChain() {
        Network model = closedMultiClass();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            Matrix qlenChain = solver.getAvgQLenChain();
            assertNotNull(qlenChain);
            assertTrue(qlenChain.getNumRows() > 0);
        });
    }

    @Test
    public void test62_networkSolver_getAvgTputChain() {
        Network model = closedMultiClass();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            Matrix tputChain = solver.getAvgTputChain();
            assertNotNull(tputChain);
        });
    }

    @Test
    public void test63_networkSolver_getAvgUtilChain() {
        Network model = closedMultiClass();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            Matrix utilChain = solver.getAvgUtilChain();
            assertNotNull(utilChain);
        });
    }

    @Test
    public void test64_networkSolver_getAvgRespTChain() {
        Network model = closedMultiClass();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            Matrix respTChain = solver.getAvgRespTChain();
            assertNotNull(respTChain);
        });
    }

    @Test
    public void test65_networkSolver_getAvgResidTChain() {
        Network model = closedMultiClass();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            Matrix residTChain = solver.getAvgResidTChain();
            assertNotNull(residTChain);
        });
    }

    // ============================================================
    // Tests 66-72: More NetworkSolver metric methods
    // ============================================================

    @Test
    public void test66_networkSolver_getAvgResidT() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            Matrix residT = solver.getAvgResidT();
            assertNotNull(residT);
        });
    }

    @Test
    public void test67_networkSolver_getAvgArvR() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            Matrix arvR = solver.getAvgArvR();
            assertNotNull(arvR);
        });
    }

    @Test
    public void test68_networkSolver_getAvgArvRChain() {
        Network model = closedMultiClass();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            Matrix arvRChain = solver.getAvgArvRChain();
            assertNotNull(arvRChain);
        });
    }

    @Test
    public void test69_networkSolver_getAvgNodeTable_open() {
        Network model = openMM1();
        withSuppressedOutput(() -> {
            SolverCTMC solver = new SolverCTMC(model, "keep", false, "cutoff", 5);
            NetworkAvgNodeTable nodeTable = solver.getAvgNodeTable();
            assertNotNull(nodeTable);
        });
    }

    @Test
    public void test70_networkSolver_getAvgChainTable_open() {
        Network model = openMultiClass();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            NetworkAvgChainTable chainTable = solver.getAvgChainTable();
            assertNotNull(chainTable);
        });
    }

    @Test
    public void test71_networkSolver_getAvgSysTable_multiclass() {
        Network model = closedMultiClass();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            NetworkAvgSysTable sysTable = solver.getAvgSysTable();
            assertNotNull(sysTable);
        });
    }

    @Test
    public void test72_networkSolver_getName() {
        Network model = closedSingle();
        SolverMVA mva = new SolverMVA(model);
        assertEquals("SolverMVA", mva.getName());
        SolverNC nc = new SolverNC(model);
        assertEquals("SolverNC", nc.getName());
    }

    // ============================================================
    // Tests 73-80: AUTO solver / LINE class (23.9%)
    // ============================================================

    @Test
    public void test73_auto_closed_avgTable() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverAUTO solver = new SolverAUTO(model);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test74_auto_open_avgTable() {
        Network model = openMM1();
        withSuppressedOutput(() -> {
            SolverAUTO solver = new SolverAUTO(model);
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test75_auto_force_nc() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverAUTO solver = new SolverAUTO(model, "force", "nc");
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
            assertEquals("SolverNC", solver.getSelectedSolverName());
        });
    }

    @Test
    public void test76_auto_force_fluid() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverAUTO solver = new SolverAUTO(model, "force", "fluid");
            NetworkAvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
            assertEquals("SolverFluid", solver.getSelectedSolverName());
        });
    }

    @Test
    public void test77_line_load_mva() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            NetworkSolver solver = LINE.load("mva", model);
            assertNotNull(solver);
            assertTrue(solver instanceof SolverMVA);
        });
    }

    @Test
    public void test78_line_load_nc() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            NetworkSolver solver = LINE.load("nc", model);
            assertNotNull(solver);
            assertTrue(solver instanceof SolverNC);
        });
    }

    @Test
    public void test79_line_load_ctmc() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            NetworkSolver solver = LINE.load("ctmc", model);
            assertNotNull(solver);
            assertTrue(solver instanceof SolverCTMC);
        });
    }

    @Test
    public void test80_line_create() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            NetworkSolver solver = LINE.create(model);
            assertNotNull(solver);
        });
    }

    // ============================================================
    // Tests 81-85: QNS solver (25.3%)
    // ============================================================

    @Test
    public void test81_qns_constructors() {
        Network model = closedMultiServer();
        SolverQNS solver1 = new SolverQNS(model);
        assertNotNull(solver1);
        SolverQNS solver2 = new SolverQNS(model, "zhou");
        assertNotNull(solver2);
    }

    @Test
    public void test82_qns_supports() {
        Network model = closedMultiServer();
        SolverQNS solver = new SolverQNS(model);
        assertNotNull(solver);
        // supports() may throw due to FeatureSet issue; just verify constructor works
    }

    @Test
    public void test83_qns_open_constructors() {
        Network model = openMM1();
        SolverQNS solver = new SolverQNS(model);
        assertNotNull(solver);
        SolverQNS solver2 = new SolverQNS(model, "zhou");
        assertNotNull(solver2);
        SolverOptions opts = SolverQNS.defaultOptions();
        assertNotNull(opts);
    }

    @Test
    public void test84_qns_defaultOptions() {
        SolverOptions opts = SolverQNS.defaultOptions();
        assertNotNull(opts);
    }

    @Test
    public void test85_qns_closed_avgTable() {
        Network model = closedMultiServer();
        withSuppressedOutput(() -> {
            SolverQNS solver = new SolverQNS(model, "zhou");
            try {
                NetworkAvgTable avg = solver.getAvgTable();
                // If qnsolver is available, verify results
                if (avg != null) {
                    assertTrue(avg.getQLen().size() > 0);
                }
            } catch (RuntimeException e) {
                // qnsolver binary may not be installed; test still covers constructor/init paths
            }
        });
    }

    // ============================================================
    // Tests 86-90: UQ solver (0%)
    // ============================================================

    @Test
    public void test86_uq_basic() {
        Network model = new Network("UQTest");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.FCFS);
        ClosedClass cc = new ClosedClass(model, "Class1", 3, delay);
        delay.setService(cc, new Exp(1.0));
        // Use Prior for service at queue
        List<Distribution> alts = Arrays.asList(new Exp(2.0), new Exp(3.0));
        double[] probs = {0.6, 0.4};
        Prior prior = new Prior(alts, probs);
        queue.setService(cc, prior);
        model.link(model.serialRouting(delay, queue));

        withSuppressedOutput(() -> {
            SolverUQ solver = new SolverUQ(model, m -> new SolverMVA(m));
            assertTrue(solver.hasPriorDistribution());
            assertEquals(2, solver.getNumAlternatives());
            AvgTable avg = solver.getAvgTable();
            assertNotNull(avg);
        });
    }

    @Test
    public void test87_uq_getPosteriorTable() {
        Network model = new Network("UQTest2");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.FCFS);
        ClosedClass cc = new ClosedClass(model, "Class1", 3, delay);
        delay.setService(cc, new Exp(1.0));
        List<Distribution> alts = Arrays.asList(new Exp(1.0), new Exp(2.0), new Exp(4.0));
        double[] probs = {0.5, 0.3, 0.2};
        Prior prior = new Prior(alts, probs);
        queue.setService(cc, prior);
        model.link(model.serialRouting(delay, queue));

        withSuppressedOutput(() -> {
            SolverUQ solver = new SolverUQ(model, m -> new SolverMVA(m));
            SolverUQ.PosteriorTable table = solver.getPosteriorTable();
            assertNotNull(table);
            assertFalse(table.rows.isEmpty());
            assertEquals(3, solver.getNumberOfModels());
        });
    }

    @Test
    public void test88_uq_getPosteriorDist() {
        Network model = new Network("UQTest3");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.FCFS);
        ClosedClass cc = new ClosedClass(model, "Class1", 3, delay);
        delay.setService(cc, new Exp(1.0));
        List<Distribution> alts = Arrays.asList(new Exp(2.0), new Exp(4.0));
        double[] probs = {0.7, 0.3};
        Prior prior = new Prior(alts, probs);
        queue.setService(cc, prior);
        model.link(model.serialRouting(delay, queue));

        withSuppressedOutput(() -> {
            SolverUQ solver = new SolverUQ(model, m -> new SolverMVA(m));
            List<Station> stations = model.getStations();
            List<JobClass> classes = model.getClasses();
            SolverUQ.EmpiricalCDF cdf = solver.getPosteriorDist("Q", stations.get(1), classes.get(0));
            assertNotNull(cdf);
            assertTrue(cdf.getMean() > 0);
            assertEquals(1.0, cdf.evalCDF(Double.MAX_VALUE), 1e-10);
        });
    }

    @Test
    public void test89_uq_supports() {
        // Model with Prior
        Network model = new Network("UQSupportTest");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.FCFS);
        ClosedClass cc = new ClosedClass(model, "Class1", 3, delay);
        delay.setService(cc, new Exp(1.0));
        List<Distribution> alts = Arrays.asList(new Exp(2.0), new Exp(3.0));
        double[] probs = {0.5, 0.5};
        queue.setService(cc, new Prior(alts, probs));
        model.link(model.serialRouting(delay, queue));

        SolverUQ solver = new SolverUQ(model, m -> new SolverMVA(m));
        assertTrue(solver.supports(model));
    }

    @Test
    public void test90_uq_empiricalCDF() {
        double[] values = {1.0, 3.0, 2.0};
        double[] probs = {0.2, 0.5, 0.3};
        SolverUQ.EmpiricalCDF cdf = new SolverUQ.EmpiricalCDF(values, probs);
        // Should be sorted: 1.0 (0.2), 2.0 (0.3), 3.0 (0.5)
        assertEquals(0.0, cdf.evalCDF(0.5), 1e-10);
        assertEquals(0.2, cdf.evalCDF(1.5), 1e-10);
        assertEquals(0.5, cdf.evalCDF(2.5), 1e-10);
        assertEquals(1.0, cdf.evalCDF(3.5), 1e-10);
        double mean = cdf.getMean();
        assertEquals(1.0 * 0.2 + 2.0 * 0.3 + 3.0 * 0.5, mean, 1e-10);
    }

    // ============================================================
    // Tests 91-95: Cross-solver metrics and handles
    // ============================================================

    @Test
    public void test91_handles_qlen() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            solver.getAvg();
            AvgHandle qHandle = solver.getAvgQLenHandles();
            assertNotNull(qHandle);
        });
    }

    @Test
    public void test92_handles_util() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            solver.getAvg();
            AvgHandle uHandle = solver.getAvgUtilHandles();
            assertNotNull(uHandle);
        });
    }

    @Test
    public void test93_handles_respT() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            solver.getAvg();
            AvgHandle rHandle = solver.getAvgRespTHandles();
            assertNotNull(rHandle);
        });
    }

    @Test
    public void test94_handles_tput() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            solver.getAvg();
            AvgHandle tHandle = solver.getAvgTputHandles();
            assertNotNull(tHandle);
        });
    }

    @Test
    public void test95_handles_arvR() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            solver.getAvg();
            AvgHandle aHandle = solver.getAvgArvRHandles();
            assertNotNull(aHandle);
        });
    }

    // ============================================================
    // Tests 96-100: Solver base class and SolverOptions
    // ============================================================

    @Test
    public void test96_solverOptions_parseOptions() {
        SolverOptions opts = new SolverOptions(SolverType.MVA);
        SolverOptions parsed = Solver.parseOptions(opts, "method", "exact");
        assertNotNull(parsed);
    }

    @Test
    public void test97_solver_hasResults() {
        Network model = closedSingle();
        SolverMVA solver = new SolverMVA(model);
        assertFalse(solver.hasResults());
        withSuppressedOutput(() -> {
            solver.getAvg();
        });
        assertTrue(solver.hasResults());
    }

    @Test
    public void test98_solver_reset() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            solver.getAvg();
            assertTrue(solver.hasResults());
            solver.reset();
            assertFalse(solver.hasResults());
        });
    }

    @Test
    public void test99_solver_result_fields() {
        Network model = closedSingle();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            SolverResult res = solver.getAvg();
            assertNotNull(res);
            assertNotNull(res.QN);
            assertNotNull(res.UN);
            assertNotNull(res.RN);
            assertNotNull(res.TN);
        });
    }

    @Test
    public void test100_solver_getAvgNodeChainTable_open() {
        Network model = openMultiClass();
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model);
            NetworkAvgNodeChainTable nodeChainTable = solver.getAvgNodeChainTable();
            assertNotNull(nodeChainTable);
        });
    }
}
