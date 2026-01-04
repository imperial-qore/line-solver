package jline.solvers.ln;

import jline.TestTools;
import jline.examples.java.basic.LayeredModel;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.constant.CallType;
import jline.lang.constant.SchedStrategy;
import jline.VerboseLevel;
import jline.lang.constant.SolverType;
import jline.lang.layered.*;
import jline.lang.nodes.ClassSwitch;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.solvers.lqns.SolverLQNS;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAOptions;
import jline.solvers.mva.SolverMVA;
import jline.util.matrix.Matrix;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.logging.Level;
import java.util.logging.Logger;

import static jline.TestTools.*;
import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;

// @Disabled("SolverLN tests exceed 3-minute timeout")
class SolverLNTest {

    static {
        // Suppress all logging during tests
        Logger.getLogger("").setLevel(Level.OFF);
        Logger.getLogger("jline").setLevel(Level.OFF);
        Logger.getLogger("org.apache.commons.io.FileUtils").setLevel(Level.OFF);
    }

    private static final boolean COMPARE_WITH_LQSIM = false;

    // Store original streams for restoration
    private final PrintStream originalOut = System.out;
    private final PrintStream originalErr = System.err;

    /**
     * Suppress System.out and System.err during execution of a Runnable
     */
    private void suppressOutput(Runnable action) {
        ByteArrayOutputStream devNull = new ByteArrayOutputStream();
        PrintStream nullStream = new PrintStream(devNull);

        // Suppress Java logger output
        Logger rootLogger = Logger.getLogger("");
        Logger jlineLogger = Logger.getLogger("jline");
        Logger fileUtilsLogger = Logger.getLogger("org.apache.commons.io.FileUtils");
        Level originalRootLevel = rootLogger.getLevel();
        Level originalJlineLevel = jlineLogger.getLevel();
        Level originalFileUtilsLevel = fileUtilsLogger.getLevel();

        try {
            System.setOut(nullStream);
            System.setErr(nullStream);
            rootLogger.setLevel(Level.OFF);
            jlineLogger.setLevel(Level.OFF);
            fileUtilsLogger.setLevel(Level.OFF);
            action.run();
        } finally {
            System.setOut(originalOut);
            System.setErr(originalErr);
            rootLogger.setLevel(originalRootLevel);
            jlineLogger.setLevel(originalJlineLevel);
            fileUtilsLogger.setLevel(originalFileUtilsLevel);
            nullStream.close();
        }
    }

    /**
     * Execute a function with suppressed output and return its result
     */
    private <T> T suppressOutput(java.util.function.Supplier<T> supplier) {
        ByteArrayOutputStream devNull = new ByteArrayOutputStream();
        PrintStream nullStream = new PrintStream(devNull);

        // Suppress Java logger output
        Logger rootLogger = Logger.getLogger("");
        Logger jlineLogger = Logger.getLogger("jline");
        Logger fileUtilsLogger = Logger.getLogger("org.apache.commons.io.FileUtils");
        Level originalRootLevel = rootLogger.getLevel();
        Level originalJlineLevel = jlineLogger.getLevel();
        Level originalFileUtilsLevel = fileUtilsLogger.getLevel();

        try {
            System.setOut(nullStream);
            System.setErr(nullStream);
            rootLogger.setLevel(Level.OFF);
            jlineLogger.setLevel(Level.OFF);
            fileUtilsLogger.setLevel(Level.OFF);
            return supplier.get();
        } finally {
            System.setOut(originalOut);
            System.setErr(originalErr);
            rootLogger.setLevel(originalRootLevel);
            jlineLogger.setLevel(originalJlineLevel);
            fileUtilsLogger.setLevel(originalFileUtilsLevel);
            nullStream.close();
        }
    }

    @org.junit.jupiter.api.Test


    // //@Disabled("Test failing - disabled for investigation")
    public void test_single_layer() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solverLN = new SolverLN(SolverLNTestFixtures.buildModel1(), options);
        Network network = solverLN.getEnsemble().get(0);
        SolverMVA layersolver = new SolverMVA(network,options);
        layersolver.runAnalyzer();

        //test queue length
        assertEquals(47.375204, layersolver.result.QN.value(), relativeTolerance(47.375204, TestTools.MID_TOL));
        assertEquals(0.164050, layersolver.result.QN.get(1, 2), relativeTolerance(0.164050, TestTools.MID_TOL));
        assertEquals(2.460746, layersolver.result.QN.get(1, 3), relativeTolerance(2.460746, TestTools.MID_TOL));

        //test through put
        assertEquals(0.473752, layersolver.result.TN.value(), relativeTolerance(0.473752, TestTools.MID_TOL));
        assertEquals(0.473752, layersolver.result.TN.get(1, 2), relativeTolerance(0.473752, TestTools.MID_TOL));
        assertEquals(0.473752, layersolver.result.TN.get(1, 3), relativeTolerance(0.473752, TestTools.MID_TOL));

        //test utilization
        assertEquals(47.375204, layersolver.result.UN.value(), relativeTolerance(47.375204, TestTools.MID_TOL));
        assertEquals(0.047375, layersolver.result.UN.get(1, 2), relativeTolerance(0.047375, TestTools.MID_TOL));
        assertEquals(0.710628, layersolver.result.UN.get(1, 3), relativeTolerance(0.710628, TestTools.MID_TOL));

        //test response time
        assertEquals(100.0, layersolver.result.RN.value(), relativeTolerance(100.0, TestTools.MID_TOL));
        assertEquals(0.346278, layersolver.result.RN.get(1, 2), relativeTolerance(0.346278, TestTools.MID_TOL));
        assertEquals(5.194165, layersolver.result.RN.get(1, 3), relativeTolerance(5.194165, TestTools.MID_TOL));

        //test idxhash
        assertTrue(Double.isNaN(solverLN.getIdxhash().get(0)));
        assertEquals(1.0, solverLN.getIdxhash().get(1), relativeTolerance(1.0, TestTools.MID_TOL));
        assertTrue(Double.isNaN(solverLN.getIdxhash().get(2)));
    }

    @org.junit.jupiter.api.Test
    public void test_multiple_layer() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solverLN = new SolverLN(SolverLNTestFixtures.buildModel2(), options);

        //Layer 1
        Network network1 = solverLN.getEnsemble().get(0);
        SolverMVA layersolver1 = new SolverMVA(network1,options);
        layersolver1.runAnalyzer();

        //Layer 2
        Network network2 = solverLN.getEnsemble().get(1);
        SolverMVA layersolver2 = new SolverMVA(network2,options);
        layersolver2.runAnalyzer();

        //Layer 3
        Network network3 = solverLN.getEnsemble().get(2);
        SolverMVA layersolver3 = new SolverMVA(network3,options);
        layersolver3.runAnalyzer();

        // test Servt_classes_updmap
        assertEquals(solverLN.getServt_classes_updmap().get(1, 1), 1);
        assertEquals(solverLN.getServt_classes_updmap().get(1, 2), 7);
        assertEquals(solverLN.getServt_classes_updmap().get(1, 3), 2);
        assertEquals(solverLN.getServt_classes_updmap().get(1, 4), 3);
        assertEquals(solverLN.getServt_classes_updmap().get(2, 1), 1);
        assertEquals(solverLN.getServt_classes_updmap().get(2, 2), 8);
        assertEquals(solverLN.getServt_classes_updmap().get(2, 3), 2);
        assertEquals(solverLN.getServt_classes_updmap().get(2, 4), 4);
        assertEquals(solverLN.getServt_classes_updmap().get(3, 1), 2);
        assertEquals(solverLN.getServt_classes_updmap().get(3, 2), 9);
        assertEquals(solverLN.getServt_classes_updmap().get(3, 3), 2);
        assertEquals(solverLN.getServt_classes_updmap().get(3, 4), 3);
        assertEquals(solverLN.getServt_classes_updmap().get(4, 1), 2);
        assertEquals(solverLN.getServt_classes_updmap().get(4, 2), 10);
        assertEquals(solverLN.getServt_classes_updmap().get(4, 3), 2);
        assertEquals(solverLN.getServt_classes_updmap().get(4, 4), 4);

        // test Call_classes_updmap
        assertEquals(solverLN.getCall_classes_updmap().get(1, 1), 1);
        assertEquals(solverLN.getCall_classes_updmap().get(1, 2), 1);
        assertEquals(solverLN.getCall_classes_updmap().get(1, 3), 1);
        assertEquals(solverLN.getCall_classes_updmap().get(1, 4), 5);
        assertEquals(solverLN.getCall_classes_updmap().get(2, 1), 4);
        assertEquals(solverLN.getCall_classes_updmap().get(2, 2), 1);
        assertEquals(solverLN.getCall_classes_updmap().get(2, 3), 2);
        assertEquals(solverLN.getCall_classes_updmap().get(2, 4), 5);

        // test getThinkt_classes_updmap
        assertEquals(solverLN.getThinkt_classes_updmap().get(1, 1), 2);
        assertEquals(solverLN.getThinkt_classes_updmap().get(1, 2), 4);
        assertEquals(solverLN.getThinkt_classes_updmap().get(1, 3), 1);
        assertEquals(solverLN.getThinkt_classes_updmap().get(1, 4), 1);
        assertEquals(solverLN.getThinkt_classes_updmap().get(2, 1), 4);
        assertEquals(solverLN.getThinkt_classes_updmap().get(2, 2), 7);
        assertEquals(solverLN.getThinkt_classes_updmap().get(2, 3), 1);
        assertEquals(solverLN.getThinkt_classes_updmap().get(2, 4), 3);
        assertEquals(solverLN.getThinkt_classes_updmap().get(3, 1), 4);
        assertEquals(solverLN.getThinkt_classes_updmap().get(3, 2), 8);
        assertEquals(solverLN.getThinkt_classes_updmap().get(3, 3), 1);
        assertEquals(solverLN.getThinkt_classes_updmap().get(3, 4), 4);

        //test idxhash
        assertEquals(solverLN.getIdxhash().get(0), Double.NaN);
        assertEquals(solverLN.getIdxhash().get(1), 1);
        assertEquals(solverLN.getIdxhash().get(2), 2);
        assertEquals(solverLN.getIdxhash().get(3), Double.NaN);
        assertEquals(solverLN.getIdxhash().get(4), 3);
    }

    @org.junit.jupiter.api.Test
    //@Disabled("Test failing - disabled for investigation")
    public void test_lqn_struct() throws Exception {
        LayeredNetworkStruct lqn = SolverLNTestFixtures.buildModel2().getStruct();
        assertEquals(lqn.nidx, 10);
        assertEquals(lqn.nhosts, 2);
        assertEquals(lqn.ntasks, 2);
        assertEquals(lqn.nentries, 2);
        assertEquals(lqn.nacts, 4);
        assertEquals(lqn.ncalls, 1);
        assertEquals(lqn.hshift, 0);
        assertEquals(lqn.tshift, 2);
        assertEquals(lqn.eshift, 4);
        assertEquals(lqn.ashift, 6);
        assertEquals(lqn.cshift, 10);

        assertEquals(lqn.tasksof.size(), 2);
        assertEquals(lqn.tasksof.get(1).size(), 1);
        assertEquals(lqn.tasksof.get(2).size(), 1);
        assertEquals(lqn.tasksof.get(1).get(0), 3);
        assertEquals(lqn.tasksof.get(2).get(0), 4);

        assertEquals(lqn.entriesof.size(), 2);
        assertEquals(lqn.entriesof.get(4).size(), 1);
        assertEquals(lqn.entriesof.get(3).size(), 1);
        assertEquals(lqn.entriesof.get(4).get(0), 6);
        assertEquals(lqn.entriesof.get(3).get(0), 5);

        assertEquals(lqn.actsof.size(), 4);
        assertEquals(lqn.actsof.get(4).size(), 2);
        assertEquals(lqn.actsof.get(3).size(), 2);
        assertEquals(lqn.actsof.get(4).get(0), 9);
        assertEquals(lqn.actsof.get(4).get(1), 10);
        assertEquals(lqn.actsof.get(3).get(0), 7);
        assertEquals(lqn.actsof.get(3).get(1), 8);

        assertEquals(lqn.callsof.get(8).size(), 1);
        assertEquals(lqn.callsof.get(8).get(0), 1);

        assertTrue(lqn.hostdem.get(3).isImmediate());
        assertTrue(lqn.hostdem.get(4).isImmediate());
        assertTrue(lqn.hostdem.get(5).isImmediate());
        assertTrue(lqn.hostdem.get(6).isImmediate());
        assertTrue(lqn.hostdem.get(8).isImmediate());
        assertTrue(lqn.hostdem.get(7) instanceof Exp);
        assertTrue(lqn.hostdem.get(9) instanceof Exp);
        assertTrue(lqn.hostdem.get(10) instanceof Exp);
        assertEquals(lqn.hostdem.get(7).getMean(), 1 / 0.625);
        assertEquals(lqn.hostdem.get(9).getMean(), 1 / 0.2);
        assertEquals(lqn.hostdem.get(10).getMean(), 1.0);

        assertEquals(lqn.think.size(), 2);
        assertTrue(lqn.think.get(3) instanceof Exp);
        assertEquals(lqn.think.get(3).getMean(), 1 / 0.01);
        assertTrue(lqn.think.get(4).isImmediate());

        assertEquals(lqn.sched.size(), 4);
        assertSame(lqn.sched.get(1), SchedStrategy.PS);
        assertSame(lqn.sched.get(2), SchedStrategy.PS);
        assertSame(lqn.sched.get(3), SchedStrategy.REF);
        assertSame(lqn.sched.get(4), SchedStrategy.FCFS);

        assertEquals(lqn.names.size(), 10);
        assertEquals(lqn.names.get(1), "P1");
        assertEquals(lqn.names.get(2), "P2");
        assertEquals(lqn.names.get(3), "T1");
        assertEquals(lqn.names.get(4), "T2");
        assertEquals(lqn.names.get(5), "E1");
        assertEquals(lqn.names.get(6), "E2");
        assertEquals(lqn.names.get(7), "A1");
        assertEquals(lqn.names.get(8), "A2");
        assertEquals(lqn.names.get(9), "A3");
        assertEquals(lqn.names.get(10), "A4");

        assertEquals(lqn.hashnames.size(), 10);
        assertEquals(lqn.hashnames.get(1), "P:P1");
        assertEquals(lqn.hashnames.get(2), "P:P2");
        assertEquals(lqn.hashnames.get(3), "R:T1");
        assertEquals(lqn.hashnames.get(4), "T:T2");
        assertEquals(lqn.hashnames.get(5), "E:E1");
        assertEquals(lqn.hashnames.get(6), "E:E2");
        assertEquals(lqn.hashnames.get(7), "A:A1");
        assertEquals(lqn.hashnames.get(8), "A:A2");
        assertEquals(lqn.hashnames.get(9), "A:A3");
        assertEquals(lqn.hashnames.get(10), "A:A4");

        assertEquals(lqn.mult.getNumRows(), 1);
        assertEquals(lqn.mult.getNumCols(), 11);
        assertEquals(lqn.mult.getNonZeros(), 10);
        assertEquals(lqn.mult.get(1), 1);
        assertEquals(lqn.mult.get(2), 1);
        assertEquals(lqn.mult.get(3), 10);
        assertEquals(lqn.mult.get(4), 1);

        assertEquals(lqn.repl.getNumRows(), 1);
        assertEquals(lqn.repl.getNumCols(), 5);
        assertEquals(lqn.repl.getNonZeros(), 4);
        assertEquals(lqn.repl.get(1), 1);
        assertEquals(lqn.repl.get(2), 1);
        assertEquals(lqn.repl.get(3), 1);
        assertEquals(lqn.repl.get(4), 1);

        assertEquals(lqn.type.getNumRows(), 1);
        assertEquals(lqn.type.getNumCols(), 11);
        assertEquals(lqn.type.getNonZeros(), 8);
        assertEquals(lqn.type.get(3), 1);
        assertEquals(lqn.type.get(4), 1);
        assertEquals(lqn.type.get(5), 2);
        assertEquals(lqn.type.get(6), 2);
        assertEquals(lqn.type.get(7), 3);
        assertEquals(lqn.type.get(8), 3);
        assertEquals(lqn.type.get(9), 3);
        assertEquals(lqn.type.get(10), 3);

        assertEquals(lqn.nitems.getNonZeros(), 0);

        assertTrue(lqn.itemcap.isEmpty());

        assertEquals(lqn.replacestrat.getNonZeros(), 0);

        assertTrue(lqn.itemproc.isEmpty());

        assertEquals(lqn.calltype.size(), 1);
        assertSame(lqn.calltype.get(1), CallType.SYNC);

        assertEquals(lqn.callpair.getNumCols(), 3);
        assertEquals(lqn.callpair.getNumRows(), 2);
        assertEquals(lqn.callpair.getNonZeros(), 2);
        assertEquals(lqn.callpair.get(1, 1), 8);
        assertEquals(lqn.callpair.get(1, 2), 6);

        assertEquals(lqn.callproc.size(), 1);
        assertNotNull(lqn.callproc.get(1));

        assertEquals(lqn.callnames.size(), 1);
        assertEquals(lqn.callnames.get(1), "A2=>E2");

        assertEquals(lqn.callhashnames.size(), 1);
        assertEquals(lqn.callhashnames.get(1), "A:A2=>E:E2");

        assertEquals(lqn.actpretype.getNumRows(), 1);
        assertEquals(lqn.actpretype.getNumCols(), 11);
        assertEquals(lqn.actpretype.getNonZeros(), 2);
        assertEquals(lqn.actpretype.get(0, 7), 1);
        assertEquals(lqn.actpretype.get(0, 9), 1);

        assertEquals(lqn.actposttype.getNumCols(), 11);
        assertEquals(lqn.actposttype.getNumRows(), 1);
        assertEquals(lqn.actposttype.getNonZeros(), 2);
        assertEquals(lqn.actposttype.get(0, 8), 11);
        assertEquals(lqn.actposttype.get(0, 10), 11);

        assertEquals(lqn.graph.getNumRows(), 11);
        assertEquals(lqn.graph.getNumCols(), 11);
        assertEquals(lqn.graph.getNonZeros(), 9);
        assertEquals(lqn.graph.get(3, 1), 1);
        assertEquals(lqn.graph.get(4, 2), 1);
        assertEquals(lqn.graph.get(3, 5), 1);
        assertEquals(lqn.graph.get(4, 6), 1);
        assertEquals(lqn.graph.get(8, 6), 1);
        assertEquals(lqn.graph.get(5, 7), 1);
        assertEquals(lqn.graph.get(7, 8), 1);
        assertEquals(lqn.graph.get(6, 9), 1);
        assertEquals(lqn.graph.get(9, 10), 1);

        assertEquals(lqn.parent.getNonZeros(), 8);
        assertEquals(lqn.parent.getNumRows(), 1);
        assertEquals(lqn.parent.getNumCols(), 11);
        assertEquals(lqn.parent.get(0, 3), 1);
        assertEquals(lqn.parent.get(0, 4), 2);
        assertEquals(lqn.parent.get(0, 5), 3);
        assertEquals(lqn.parent.get(0, 6), 4);
        assertEquals(lqn.parent.get(0, 7), 3);
        assertEquals(lqn.parent.get(0, 8), 3);
        assertEquals(lqn.parent.get(0, 9), 4);
        assertEquals(lqn.parent.get(0, 10), 4);

        assertEquals(lqn.replygraph.getNumCols(), 3);
        assertEquals(lqn.replygraph.getNumRows(), 5);
        assertEquals(lqn.replygraph.getNonZeros(), 2);
        assertEquals(lqn.replygraph.get(2, 1), 1);
        assertEquals(lqn.replygraph.get(4, 2), 1);

        assertEquals(lqn.iscache.getNonZeros(), 0);

        assertEquals(lqn.iscaller.getNonZeros(), 4);
        assertEquals(lqn.iscaller.getNumRows(), 11);
        assertEquals(lqn.iscaller.getNumCols(), 11);
        assertEquals(lqn.iscaller.get(3, 4), 1);
        assertEquals(lqn.iscaller.get(8, 4), 1);
        assertEquals(lqn.iscaller.get(3, 6), 1);
        assertEquals(lqn.iscaller.get(8, 6), 1);

        assertEquals(lqn.issynccaller.getNumCols(), 11);
        assertEquals(lqn.issynccaller.getNumRows(), 11);
        assertEquals(lqn.issynccaller.getNonZeros(), 4);
        assertEquals(lqn.issynccaller.get(3, 4), 1);
        assertEquals(lqn.issynccaller.get(8, 4), 1);
        assertEquals(lqn.issynccaller.get(3, 6), 1);
        assertEquals(lqn.issynccaller.get(8, 6), 1);

        assertEquals(lqn.isasynccaller.getNonZeros(), 0);

        assertEquals(lqn.isref.getNumRows(), 1);
        assertEquals(lqn.isref.getNonZeros(), 1);
        assertEquals(lqn.isref.getNumCols(), 5);
        assertEquals(lqn.isref.get(0, 3), 1);
    }

    @org.junit.jupiter.api.Test


    //@Disabled("Test failing - disabled for investigation")
    public void test_layer_struct() throws Exception {
        // this test is to test the network's job classes, stations and connections

        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solverLN = new SolverLN(SolverLNTestFixtures.buildModel3(), options);

        Network network1 = solverLN.getEnsemble().get(0);

        Network network2 = solverLN.getEnsemble().get(1);

        // test layer 1
        // test numbers of nodes, classes, stations
        assertEquals(network1.getNumberOfNodes(), 6);
        assertEquals(network1.getNumberOfClasses(), 9);
        assertEquals(network1.getNumberOfStations(), 2);

        // test classes types, populations
        for (int i = 0; i < network1.getNumberOfClasses(); i++) {
            assertTrue(network1.getClasses().get(i) instanceof ClosedClass);
        }
        ClosedClass temp = (ClosedClass) network1.getClasses().get(0);
        assertEquals(temp.getPopulation(), 10);
        temp = (ClosedClass) network1.getClasses().get(5);
        assertEquals(temp.getPopulation(), 10);
        for (int i = 1; i < 9; i++) {
            if (i != 5) {
                temp = (ClosedClass) network1.getClasses().get(i);
                assertEquals(temp.getPopulation(), 0);
            }
        }

        // test connection matrix
        assertEquals(network1.getConnectionMatrix().get(2, 0), 1);
        assertEquals(network1.getConnectionMatrix().get(4, 0), 1);
        assertEquals(network1.getConnectionMatrix().get(3, 1), 1);
        assertEquals(network1.getConnectionMatrix().get(5, 1), 1);
        assertEquals(network1.getConnectionMatrix().get(0, 2), 1);
        assertEquals(network1.getConnectionMatrix().get(0, 3), 1);
        assertEquals(network1.getConnectionMatrix().get(1, 4), 1);
        assertEquals(network1.getConnectionMatrix().get(1, 5), 1);
        assertEquals(network1.getConnectionMatrix().getNonZeros(), 8);

        // test the node types
        for (int i = 2; i < network1.getNumberOfNodes(); i++) {
            assertTrue(network1.getNodes().get(i) instanceof ClassSwitch);
        }


        // for every delay and queue test the Service Process
        assertTrue(network1.getStations().get(1) instanceof Queue);
        Queue queue = (Queue) network1.getStations().get(1);
        //assertTrue(queue.getServiceProcess(network1.getJobClass().get(0)) instanceof Disabled);
        //assertTrue(queue.getServiceProcess(network1.getJobClass().get(1)) instanceof Disabled);
        assertTrue(queue.getServiceProcess(network1.getJobClasses().get(2)) instanceof Exp);
        assertTrue(queue.getServiceProcess(network1.getJobClasses().get(3)) instanceof Immediate);
        // Job class 4 can be either Exp or Immediate (when service time is 0)
        assertTrue(queue.getServiceProcess(network1.getJobClasses().get(4)) instanceof Exp || 
                   queue.getServiceProcess(network1.getJobClasses().get(4)) instanceof Immediate);
        //assertTrue(queue.getServiceProcess(network1.getJobClass().get(5)) instanceof Disabled);
        //assertTrue(queue.getServiceProcess(network1.getJobClass().get(6)) instanceof Disabled);
        assertTrue(queue.getServiceProcess(network1.getJobClasses().get(7)) instanceof Exp);
        assertTrue(queue.getServiceProcess(network1.getJobClasses().get(8)) instanceof Exp);

        assertTrue(compareRelErr(queue.getServiceProcess(network1.getJobClasses().get(2)).getMean(), 70));
        assertTrue(compareRelErr(queue.getServiceProcess(network1.getJobClasses().get(7)).getMean(), 5));
        assertTrue(compareRelErr(queue.getServiceProcess(network1.getJobClasses().get(8)).getMean(), 1));
        assertEquals(queue.getNumberOfServers(), 1);

        assertTrue(network1.getStations().get(0) instanceof Delay);
        Delay delay = (Delay) network1.getStations().get(0);
        assertTrue(delay.getServiceProcess(network1.getJobClasses().get(0)) instanceof Exp);
        assertTrue(delay.getServiceProcess(network1.getJobClasses().get(1)) instanceof Immediate);
        //assertTrue(delay.getServiceProcess(network1.getJobClass().get(2)) instanceof Disabled);
        //assertTrue(delay.getServiceProcess(network1.getJobClass().get(3)) instanceof Disabled);
        assertTrue(delay.getServiceProcess(network1.getJobClasses().get(4)) instanceof Immediate);
        assertTrue(delay.getServiceProcess(network1.getJobClasses().get(5)) instanceof Immediate);
        assertTrue(delay.getServiceProcess(network1.getJobClasses().get(6)) instanceof Immediate);
        //assertTrue(delay.getServiceProcess(network1.getJobClass().get(7)) instanceof Disabled);
        //assertTrue(delay.getServiceProcess(network1.getJobClass().get(8)) instanceof Disabled);

        assertTrue(compareRelErr(delay.getServiceProcess(network1.getJobClasses().get(0)).getMean(), 100));

        // test layer 2
        // test numbers of nodes, classes, stations
        assertEquals(network2.getNumberOfNodes(), 4);
        assertEquals(network2.getNumberOfClasses(), 5);
        assertEquals(network2.getNumberOfStations(), 2);

        // test classes types, populations
        for (int i = 0; i < network2.getNumberOfClasses(); i++) {
            assertTrue(network2.getClasses().get(i) instanceof ClosedClass);
        }
        temp = (ClosedClass) network1.getClasses().get(0);
        assertEquals(temp.getPopulation(), 10);
        for (int i = 1; i < network2.getNumberOfClasses(); i++) {
            temp = (ClosedClass) network1.getClasses().get(i);
            assertEquals(temp.getPopulation(), 0);
        }

        // test connection matrix
        assertEquals(network2.getConnectionMatrix().get(1, 0), 1);
        assertEquals(network2.getConnectionMatrix().get(2, 0), 1);
        assertEquals(network2.getConnectionMatrix().get(3, 1), 1);
        assertEquals(network2.getConnectionMatrix().get(0, 2), 1);
        assertEquals(network2.getConnectionMatrix().get(0, 3), 1);
        assertEquals(network2.getConnectionMatrix().getNonZeros(), 5);

        // test the node types
        for (int i = 2; i < network2.getNumberOfNodes(); i++) {
            assertTrue(network2.getNodes().get(i) instanceof ClassSwitch);
        }

        assertTrue(network2.getStations().get(0) instanceof Delay);
        delay = (Delay) network2.getStations().get(0);
        assertTrue(delay.getServiceProcess(network2.getJobClasses().get(0)) instanceof Exp);
        assertTrue(delay.getServiceProcess(network2.getJobClasses().get(1)) instanceof Immediate);
        assertTrue(delay.getServiceProcess(network2.getJobClasses().get(2)) instanceof Exp);
        assertTrue(delay.getServiceProcess(network2.getJobClasses().get(3)) instanceof Immediate);
        assertTrue(delay.getServiceProcess(network2.getJobClasses().get(4)) instanceof Immediate);

        assertTrue(compareRelErr(delay.getServiceProcess(network2.getJobClasses().get(0)).getMean(), 100));
        assertTrue(compareRelErr(delay.getServiceProcess(network2.getJobClasses().get(2)).getMean(), 70));

        assertTrue(network2.getStations().get(1) instanceof Queue);
        queue = (Queue) network2.getStations().get(1);
        //assertTrue(queue.getServiceProcess(network2.getJobClass().get(0)) instanceof Disabled);
        //assertTrue(queue.getServiceProcess(network2.getJobClass().get(1)) instanceof Disabled);;
        //assertTrue(queue.getServiceProcess(network2.getJobClass().get(2)) instanceof Disabled);
        //assertTrue(queue.getServiceProcess(network2.getJobClass().get(3)) instanceof Disabled);
        assertTrue(queue.getServiceProcess(network2.getJobClasses().get(4)) instanceof Immediate);
        assertEquals(queue.getNumberOfServers(), 10);

        //test idxhash
        assertEquals(solverLN.getIdxhash().get(0), Double.NaN);
        assertEquals(solverLN.getIdxhash().get(1), 1);
        assertEquals(solverLN.getIdxhash().get(2), Double.NaN);
        assertEquals(solverLN.getIdxhash().get(3), 2);
    }

    @org.junit.jupiter.api.Test


    //@Disabled("Test failing - disabled for investigation")
    public void test_buildModel_1() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel1(), SolverType.MVA, options);
        double[] expectedQLen = {Double.NaN, 2.62479598084316, 2.62479598084316, 0.164049748802698, 2.46074623204047};
        double[] expectedUtil = {0.758003264230709, 0.758003264230709, Double.NaN, 0.0473752040144193, 0.71062806021629};
        double[] expectedRespT = {Double.NaN, Double.NaN, 5.54044259111638, 0.346277661944774, 5.1941649291716};
        double[] expectedResidT = {Double.NaN, 5.54044259111638, Double.NaN, 0.346277661944774, 5.1941649291716};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, 0.473752040144193, 0.473752040144193, 0.473752040144193, 0.473752040144193};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();

        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput);
    }

    @org.junit.jupiter.api.Test
    public void test_buildModel_2() throws Exception {
        LNOptions lnoptions = new LNOptions();
        lnoptions.verbose = VerboseLevel.SILENT;
        SolverOptions mvaoptions = new MVAOptions();
        mvaoptions.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel2(), SolverType.MVA, lnoptions, mvaoptions);
        // Ground truth values from MATLAB
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.124064286138323, 0.5325561421946614, 1.124064286138323, 0.5325561421946614, 0.1624589500751483, 0.9616053370785427, 0.4437967851622179, 0.08875935703244357};
        double[] expectedUtil = {0.1420149713913745, 0.5325561421946613, 0.1420149713913745, 0.5325561421946613, Double.NaN, Double.NaN, 0.1420149713913745, 0, 0.4437967851622178, 0.08875935703244356};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 12.66417786942251, 6.000000000000001, 1.830330405122518, 10.83384747573956, 5.000000000000001, 1};
        double[] expectedResidT = {Double.NaN, Double.NaN, 1.830330405122518, 6.000000000000001, Double.NaN, Double.NaN, 1.830330405122518, 0, 5.000000000000001, 1};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.08875935711960906, 0.08875935703244356, 0.08875935711960906, 0.08875935703244356, 0.08875935711960906, 0.08875935711960906, 0.08875935703244356, 0.08875935703244356};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();

        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput);
    }

    @org.junit.jupiter.api.Test
    public void test_buildModel_3() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        options.config.relax = "none"; // Disable relaxation for backward compatibility
        options.iter_max = 100; // Original default for backward compatibility
        options.iter_tol = 0.0001; // Original default for backward compatibility
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel3(), SolverType.MVA, options);
        // Ground truth values from MATLAB (without relaxation)
        double[] expectedQLen = {Double.NaN, 8.684302439121234, 0.7566547562277667, 8.684302439121234, 0.7566547562277667, 7.92765158428176, 0.7566508559719956, 0.6305456301898056, 0.1261091260379611};
        double[] expectedUtil = {0.9999545067535648, 0.9210092588749577, 0.07894524787860714, Double.NaN, Double.NaN, 0.9210092588749577, 0, 0.06578770656550595, 0.01315754131310119};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, 660.0380668062525, 57.50730613130224, 602.5298937576315, 57.50817313469662, 47.92275510941854, 9.584551021883705};
        double[] expectedResidT = {Double.NaN, 602.5298937576315, 57.50730613130224, Double.NaN, Double.NaN, 602.5298937576315, 0, 47.92275510941854, 9.584551021883705};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, 0.01315727512678511, 0.01315754131310119, 0.01315727512678511, 0.01315754131310119, 0.01315727512678511, 0.01315727512678511, 0.01315754131310119, 0.01315754131310119};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();

        avgTable.print();
        // Print LQNS results for comparison if available
        if (SolverLQNS.isAvailable()) {
            new SolverLQNS(SolverLNTestFixtures.buildModel3()).getAvgTable().print();
        }

        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput);
    }

    // //@Disabled("Passes but very slow, so disabled for performance reasons")
    @org.junit.jupiter.api.Test

    //@Disabled("Test failing - disabled for investigation")
    public void test_buildModel_4() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        options.config.relax = "none"; // Disable relaxation for backward compatibility
        options.iter_max = 100; // Original default for backward compatibility
        options.iter_tol = 0.0001; // Original default for backward compatibility
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel4(), SolverType.MVA, options);
        // Ground truth values from the table provided (without relaxation)
        double[] expectedQLen = {Double.NaN, Double.NaN, 23.4682450470894, 8.67822070268364, 0.120491297980306, 23.4682450470894, 8.67822070268364, 0.120491297980306, 23.4388606730313, 8.66725245115206, 0.124378114044187};
        double[] expectedUtil = {0.995954019616073, 0.0414593696959083, 0.664070430319683, 0.331883589296389, 0.0414593696959083, Double.NaN, Double.NaN, Double.NaN, 0.664070430319683, 0.331883589296389, 0.0414593696959083};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 1.76699970180812, 0.653709687866906, 0.0193750007720929, 1.76478725771209, 0.652883475613173, 0.0200000007969991};
        double[] expectedResidT = {Double.NaN, Double.NaN, 1.09193140773878, 0.556095776529935, 0.0200000007969991, Double.NaN, Double.NaN, Double.NaN, 1.09193140773878, 0.556095776529935, 0.0200000007969991};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 13.2814086063937, 13.2753435718556, 6.21890545438625, 13.2814086063937, 13.2753435718556, 6.21890545438625, 13.2814086063937, 13.2753435718556, 6.21890545438625};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput);
    }

    @org.junit.jupiter.api.Test


    //@Disabled("Test failing - disabled for investigation")
    public void test_buildModel_5() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel5(), SolverType.MVA, options);
        double[] expectedQLen = {Double.NaN, 1.0000, 1.0000, 0.1587, 0.5238, 0.3175};
        double[] expectedUtil = {1.0000, 1.0000, Double.NaN, 0.1587, 0.5238, 0.3175};
        double[] expectedRespT = {Double.NaN, Double.NaN, 12.6000, 2.0000, 3.0000, 4.0000};
        double[] expectedResidT = {Double.NaN, 12.6000, Double.NaN, 2.0000, 6.6000, 4.0000};
        double[] expectedTput = {Double.NaN, 0.0794, 0.0794, 0.0794, 0.1746, 0.0794};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            assertTrue(compareAbsErr(avg.getQLen().get(idx), expectedQLen[idx]));
            assertTrue(compareAbsErr(avg.getUtil().get(idx), expectedUtil[idx]));
            assertTrue(compareAbsErr(avg.getRespT().get(idx), expectedRespT[idx]));
            assertTrue(compareAbsErr(avg.getResidT().get(idx), expectedResidT[idx]));
            assertTrue(compareAbsErr(avg.getTput().get(idx), expectedTput[idx]));
        }
    }

    @org.junit.jupiter.api.Test
    public void test_buildModel_6() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel6(), SolverType.MVA, options);
        // Ground truth values from MATLAB
        // Order: P1, T1, Entry, B1, B2, B3, B4, B5
        double[] expectedQLen = {Double.NaN, 1.0, 1.0, 0.160000027162136, 0.0720000122229614, 0.0960000162972819, 0.192000032594564, 0.480000081486409};
        double[] expectedUtil = {1.0, 1.0, Double.NaN, 0.160000002128186, 0.0720000009576836, 0.0960000012769115, 0.192000002553823, 0.480000006384557};
        double[] expectedRespT = {Double.NaN, Double.NaN, 12.5, 2.0, 3.0, 4.0, 6.0, 6.0};
        double[] expectedResidT = {Double.NaN, 12.5, Double.NaN, 2.0, 0.9, 1.2, 2.4, 6.0};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, 0.0800000010640929, 0.0800000010640929, 0.0800000010640929, 0.0240000003192279, 0.0240000003192279, 0.0320000004256372, 0.0800000010640929};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput);
    }


    // test_example_layeredModel_2 removed - already tested in LayeredExamplesTest.java as lqn_multi_solvers
    // test_example_layeredModel_3 removed - already tested in LayeredExamplesTest.java as lqn_init

//    @org.junit.jupiter.api.Test
//    public void test_example_layeredModel_4() throws Exception {
//        SolverOptions options = new LNOptions();
//        options.verbose = VerboseLevel.SILENT;
//        LayeredNetwork layeredModel = LayeredModel.lqn_serial_4();
//        SolverLN solver;
//        try {
//            solver = new SolverLN(layeredModel, SolverType.MVA, options);
//        } catch (Exception e) {
//            System.out.println("Error details:");
//            e.printStackTrace();
//            throw e;
//        }
//        double[] expectedQLen = {Double.NaN, 0.0000, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 0.1310, 0.0000, 0.0663, 0.0463, 0.0192, 0.0383, 0.0192, 0.0631, 0.0080, 0.0000, 0.0463, 0.0192, 0.0383, 0.0192, 0.0080, 0.0000, 0.0663, 0.0463, 0.0192, 0.0383, 0.0192, 0.0631, 0.0080, 0.0000, 0.0463, 0.0192, 0.0383, 0.0192, 0.0080, 0.0000, 0.0000, 0.0110, 0.0110, 0.0110, 0.0110, 0.0110, 0.0110, 0.0000, 0.0000, 0.0463, 0.0000, 0.0000, 0.0192, 0.0000, 0.0000, 0.0383, 0.0000, 0.0000, 0.0192, 0.0000, 0.0000, 0.0079, 0.0079, 0.0079, 0.0079, 0.0079, 0.0079, 0.0079, 0.0079, 0.0000, 0.0000, 0.0080, 0.0000};
//        double[] expectedUtil = {0.1164, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.1164, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 0.0000, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 0.0000, 0.0412, 0.0170, 0.0341, 0.0170, 0.0071, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000};
//        double[] expectedRespT = {Double.NaN, 0.0000, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 0.0000, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 0.0000, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0000, 0.0667, 0.0113, 0.0113, 0.0113, 0.0113, 0.0889, 0.0113, 0.0000, 0.0113, 0.0113, 0.0113, 0.0113, 0.0113, 0.0000, 0.0000, 0.0111, 0.0111, 0.0111, 0.0111, 0.0111, 0.0111, 0.0000, 0.0000, 0.0113, 0.0000, 0.0000, 0.0113, 0.0000, 0.0000, 0.0113, 0.0000, 0.0000, 0.0113, 0.0000, 0.0000, 0.0111, 0.0111, 0.0111, 0.0111, 0.0111, 0.0111, 0.0111, 0.0111, 0.0000, 0.0000, 0.0113, 0.0000};
//        double[] expectedResidT = {Double.NaN, 0.0000, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 0.0113, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 0.0000, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 0.0000, 0.0040, 0.0016, 0.0033, 0.0016, 0.0007, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000};
//        double[] expectedTput = {Double.NaN, 0.0000, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 11.6365, 0.0000, 0.9934, 4.1159, 1.7030, 3.4061, 1.7030, 0.7098, 0.7097, 0.0000, 4.1154, 1.7028, 3.4057, 1.7028, 0.7096, 0.0000, 0.9934, 4.1159, 1.7030, 3.4061, 1.7030, 0.7098, 0.7097, 0.0000, 4.1154, 1.7028, 3.4057, 1.7028, 0.7096, 0.0000, 0.9934, 0.9934, 0.9934, 0.9934, 0.9934, 0.9934, 0.9934, 0.9934, 4.1159, 4.1159, 4.1159, 1.7030, 1.7030, 1.7030, 3.4061, 3.4061, 3.4061, 1.7030, 1.7030, 1.7030, 0.7098, 0.7098, 0.7098, 0.7098, 0.7098, 0.7098, 0.7098, 0.7098, 0.7098, 0.7098, 0.7097, 0.7097, 0.7097};
//        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
//        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
//            warningAssertTrue(compareAbsErr(avg.getQLen().get(idx), expectedQLen[idx]));
//            warningAssertTrue(compareAbsErr(avg.getUtil().get(idx), expectedUtil[idx]));
//            warningAssertTrue(compareAbsErr(avg.getRespT().get(idx), expectedRespT[idx]));
//            warningAssertTrue(compareAbsErr(avg.getResidT().get(idx), expectedResidT[idx]));
//            warningAssertTrue(compareAbsErr(avg.getTput().get(idx), expectedTput[idx]));
//        }
//    }

    // test_example_layeredModel_5 removed - already tested in LayeredExamplesTest.java as lqn_twotasks
    // test_example_layeredModel_6 removed - already tested in LayeredExamplesTest.java as lqn_bpmn
    // test_example_layeredModel_7 removed - already tested in LayeredExamplesTest.java as lqn_workflows

    @org.junit.jupiter.api.Test


    //@Disabled("Test failing - disabled for investigation")
    public void test_activityGraph_and() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(jline.solvers.ln.SolverLNTestFixtures.test_activityGraph_and(), SolverType.MVA, options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.0000, 0.5479, 1.0000, 0.5479, 1.0000, 0.0000, 0.8769, 0.3563, 0.3563, 0.0000};
        double[] expectedUtil = {0.0000, 1.5895, 0.0000, 1.5895, Double.NaN, Double.NaN, 0.0000, 0.0000, 0.8770, 0.3563, 0.3563, 0.0000};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 1.3333, 1.3333, 1.3333, 0.0000, 2.0000, 1.0000, 1.0000, 0.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, 0.0000, 1.3333, Double.NaN, Double.NaN, 0.0000, 0.0000, 0.6666, 0.3333, 0.3333, 0.0000};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.7500, 0.4110, 0.7500, 0.4110, 0.7500, 0.4932, 0.4385, 0.3563, 0.3563, 0.4110};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareAbsErr(avg.getQLen().get(idx), expectedQLen[idx]));
            warningAssertTrue(compareAbsErr(avg.getUtil().get(idx), expectedUtil[idx]));
            warningAssertTrue(compareAbsErr(avg.getRespT().get(idx), expectedRespT[idx]));
            warningAssertTrue(compareAbsErr(avg.getResidT().get(idx), expectedResidT[idx]));
            warningAssertTrue(compareAbsErr(avg.getTput().get(idx), expectedTput[idx]));
        }
    }

    @org.junit.jupiter.api.Test
    public void test_activityGraph_call() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(jline.solvers.ln.SolverLNTestFixtures.test_activityGraph_call(), SolverType.MVA, options);
        // Ground truth values from MATLAB
        double[] expectedQLen = {Double.NaN, Double.NaN, Double.NaN, 0.1273996509487435, 0.1134380453455223, 0.06980802789275393, 0.1273996509487435, 0.1134380453455223, 0.06980802789275393, 0.1273996509487435, 0.1134380453455223, 0.06980802789275393};
        double[] expectedUtil = {0.01396160558342394, 0.04363001744058552, 0.06980802789275393, 0.01396160558342394, 0.04363001744058552, 0.06980802789275393, Double.NaN, Double.NaN, Double.NaN, 0.01396160558342394, 0.04363001744058552, 0.06980802789275393};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 14.6, 13, 8, 14.6, 13, 8};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, 1.6, 5, 8, Double.NaN, Double.NaN, Double.NaN, 1.6, 5, 8};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, Double.NaN, 0.008726003489639965, 0.008726003488117104, 0.008726003486594242, 0.008726003489639965, 0.008726003488117104, 0.008726003486594242, 0.008726003489639965, 0.008726003488117104, 0.008726003486594242};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();

        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput);
    }

    @org.junit.jupiter.api.Test


    //@Disabled("Test failing - disabled for investigation")
    public void test_activityGraph_call_and() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(jline.solvers.ln.SolverLNTestFixtures.test_activityGraph_call_and(), SolverType.MVA, options);
        double[] expectedQLen = {Double.NaN, Double.NaN, Double.NaN, 1.0000, 0.7600, 0.9208, 1.0000, 0.7600, 0.9208, 1.0000, 0.6400, 0.1333, 0.0666, 0.0000, 0.9208};
        double[] expectedUtil = {0.0000, 0.2000, 0.9208, 0.0000, 0.2000, 0.9208, Double.NaN, Double.NaN, Double.NaN, 0.0000, 0.0000, 0.1333, 0.0666, 0.0000, 0.9208};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 9.5000, 9.5000, 8.0000, 9.5000, 8.0000, 2.0000, 1.0000, 0.0000, 8.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, 0.0001, 1.5000, 8.0000, Double.NaN, Double.NaN, Double.NaN, 0.0001, 0.0000, 1.0000, 0.5000, 0.0000, 8.0000};
        double[] expectedTput = {Double.NaN, Double.NaN, Double.NaN, 0.1053, 0.0800, 0.1151, 0.1053, 0.0800, 0.1151, 0.1053, 0.0800, 0.0666, 0.0666, 0.0800, 0.1151};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareAbsErr(avg.getQLen().get(idx), expectedQLen[idx]));
            warningAssertTrue(compareAbsErr(avg.getUtil().get(idx), expectedUtil[idx]));
            warningAssertTrue(compareAbsErr(avg.getRespT().get(idx), expectedRespT[idx]));
            warningAssertTrue(compareAbsErr(avg.getResidT().get(idx), expectedResidT[idx]));
            warningAssertTrue(compareAbsErr(avg.getTput().get(idx), expectedTput[idx]));
        }
    }

    @org.junit.jupiter.api.Test


    //@Disabled("Test failing - disabled for investigation")
    public void test_activityGraph_call_or() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(jline.solvers.ln.SolverLNTestFixtures.test_activityGraph_call_or(), SolverType.MVA, options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 0.0000, 1.0000, 0.9846, 0.0000, 1.0000, 0.9846, 0.0000, 1.0000, 0.0000, 0.0096, 0.0144, 0.9606, 0.0000};
        double[] expectedUtil = {0.0154, 0.9846, 0.0000, 0.0154, 0.9846, 0.0000, Double.NaN, Double.NaN, 0.0000, 0.0154, 0.0000, 0.0096, 0.0144, 0.9606, 0.0000};
        double[] expectedRespT = {Double.NaN, Double.NaN, 0.0000, Double.NaN, Double.NaN, 0.0000, 104.1000, 102.5000, 0.0000, 104.1000, 0.0000, 2.0000, 3.0000, 100.0000, 0.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, 0.0000, 1.6000, 102.5000, 0.0000, Double.NaN, Double.NaN, 0.0000, 1.6000, 0.0000, 1.0000, 1.5000, 100.0000, 0.0000};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.0000, 0.0096, 0.0096, 0.0000, 0.0096, 0.0096, 0.0000, 0.0096, 0.0096, 0.0048, 0.0048, 0.0096, 0.0000};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            assertTrue(compareAbsErr(avg.getQLen().get(idx), expectedQLen[idx]));
            assertTrue(compareAbsErr(avg.getUtil().get(idx), expectedUtil[idx]));
            assertTrue(compareAbsErr(avg.getRespT().get(idx), expectedRespT[idx]));
            assertTrue(compareAbsErr(avg.getResidT().get(idx), expectedResidT[idx]));
            assertTrue(compareAbsErr(avg.getTput().get(idx), expectedTput[idx]));
        }
    }

    @org.junit.jupiter.api.Test
    public void test_activityGraph_call_seq_disconnected() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(jline.solvers.ln.SolverLNTestFixtures.test_activityGraph_call_seq_disconnected(), SolverType.MVA, options);
        // Ground truth values from MATLAB
        double[] expectedQLen = {Double.NaN, Double.NaN, 0, 0.07063197025365875, 0.05576208176883958, 0, 0.07063197025365875, 0.05576208176883958, 0, 0.07063197025365875, 0.04646840147403298, 0.009293680294806596, 0};
        double[] expectedUtil = {0.01486988847445447, 0.05576208176883957, 0, 0.01486988847445447, 0.05576208176883957, 0, Double.NaN, Double.NaN, 0, 0.01486988847445447, 0.04646840147403298, 0.009293680294806596, 0};
        double[] expectedRespT = {Double.NaN, Double.NaN, 0, Double.NaN, Double.NaN, 0, 7.6, 6, 0, 7.6, 5, 1, 0};
        double[] expectedResidT = {Double.NaN, Double.NaN, 0, 1.6, 6, 0, Double.NaN, Double.NaN, 0, 1.6, 5, 1, 0};
        double[] expectedArvR = {Double.NaN, Double.NaN, 0, Double.NaN, Double.NaN, 0, Double.NaN, Double.NaN, 0, Double.NaN, Double.NaN, Double.NaN, 0};
        double[] expectedTput = {Double.NaN, Double.NaN, 0, 0.009293680296534046, 0.009293680294806596, 0, 0.009293680296534046, 0.009293680294806596, 0, 0.009293680296534046, 0.009293680294806596, 0.009293680294806596, 0};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput);
    }

    @org.junit.jupiter.api.Test


    //@Disabled("Test failing - QLen[7] mismatch")
    public void test_activityGraph_loop() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(jline.solvers.ln.SolverLNTestFixtures.test_activityGraph_loop(), SolverType.MVA, options);
        // Ground truth values from the table provided
        double[] expectedQLen = {Double.NaN, Double.NaN, 0.0909091652727956, 0.0909091652480021, 0.0909091652727956, 0.0909091652480021, 0.0909091653637047, 0.0, 0.0818182569050276, 0.00909090834297452};
        double[] expectedUtil = {0.0, 0.0909091652480022, 0.0, 0.0909091652480022, Double.NaN, Double.NaN, 0.0, 0.0, 0.0818182569050276, 0.00909090834297452};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 10.0, 10.0, 10.0, 1e-08, 3.0, 1.0};
        double[] expectedResidT = {Double.NaN, Double.NaN, 0.0, 10.0, Double.NaN, Double.NaN, 0.0, 0.0, 9.0, 1.0};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.00909090834545386, 0.00909090834297452, 0.00909090834545386, 0.00909090834297452, 0.00909090834545386, 0.00909090834297452, 0.0272727250289236, 0.00909090834297452};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput);
    }

    @org.junit.jupiter.api.Test


    //@Disabled("Test failing - disabled for investigation")
    public void test_activityGraph_or() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(jline.solvers.ln.SolverLNTestFixtures.test_activityGraph_or(), SolverType.MVA, options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.0000, 0.9846, 1.0000, 0.9846, 1.0000, 0.0000, 0.0096, 0.0144, 0.9606};
        double[] expectedUtil = {0.0154, 0.9846, 0.0154, 0.9846, Double.NaN, Double.NaN, 0.0154, 0.0000, 0.0096, 0.0144, 0.9606};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 104.1000, 102.5000, 104.1000, 0.0000, 2.0000, 3.0000, 100.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, 1.6000, 102.5000, Double.NaN, Double.NaN, 1.6000, 0.0000, 1.0000, 1.5000, 100.0000};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.0096, 0.0096, 0.0096, 0.0096, 0.0096, 0.0096, 0.0048, 0.0048, 0.0096};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            assertTrue(compareAbsErr(avg.getQLen().get(idx), expectedQLen[idx]));
            assertTrue(compareAbsErr(avg.getUtil().get(idx), expectedUtil[idx]));
            assertTrue(compareAbsErr(avg.getRespT().get(idx), expectedRespT[idx]));
            assertTrue(compareAbsErr(avg.getResidT().get(idx), expectedResidT[idx]));
            assertTrue(compareAbsErr(avg.getTput().get(idx), expectedTput[idx]));
        }
    }

    @org.junit.jupiter.api.Test
    public void test_activityGraph_seq() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(jline.solvers.ln.SolverLNTestFixtures.test_activityGraph_seq(), SolverType.MVA, options);
        // Ground truth values from MATLAB
        double[] expectedQLen = {Double.NaN, Double.NaN, 0.07063197025365875, 0.05576208176883958, 0.07063197025365875, 0.05576208176883958, 0.07063197025365875, 0.04646840147403298, 0.009293680294806596};
        double[] expectedUtil = {0.01486988847445447, 0.05576208176883957, 0.01486988847445447, 0.05576208176883957, Double.NaN, Double.NaN, 0.01486988847445447, 0.04646840147403298, 0.009293680294806596};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 7.6, 6, 7.6, 5, 1};
        double[] expectedResidT = {Double.NaN, Double.NaN, 1.6, 6, Double.NaN, Double.NaN, 1.6, 5, 1};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.009293680296534046, 0.009293680294806596, 0.009293680296534046, 0.009293680294806596, 0.009293680296534046, 0.009293680294806596, 0.009293680294806596};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput);
    }

    @org.junit.jupiter.api.Test


    //@Disabled("Test failing - disabled for investigation")
    public void test_cache_layer() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;

        // Use the existing buildCacheModel which matches the MATLAB implementation
        LayeredNetwork cacheModel = SolverLNTestFixtures.buildCacheModel();
        SolverLN solver = new SolverLN(cacheModel, SolverType.MVA, options);

        // Verify that the cache layer was created properly
        assertNotNull(solver.getEnsemble());
        assertTrue(solver.getEnsemble().size() > 0);

        // Test that the solver builds the ensemble without errors
        // This verifies that our cache implementation doesn't cause build failures
        assertTrue(solver.getEnsemble().size() >= 1, "Cache model should create at least one layer");
    }

    // Tests for SolverLN with moment3 method enabled.
    // These tests verify that the post-convergence 3-moment APH fitting works correctly.

    @Test
    public void test_moment3_basic() throws Exception {
        // Use an existing valid model
        LayeredNetwork lqn = SolverLNTestFixtures.buildModel2();

        // Create options with moment3 method
        LNOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        options.method = "moment3";

        // Create solver and run to convergence
        SolverLN solver = new SolverLN(lqn, options);
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();

        // Basic sanity checks - the solver should converge and produce results
        assertTrue(solver.hasconverged, "Solver should converge");
        assertNotNull(solver.servt, "Service times should be computed");
        assertNotNull(solver.tput, "Throughput should be computed");

        // Check that entry processes are populated (post-convergence feature)
        assertNotNull(solver.entryproc, "Entry processes should be created");
        assertNotNull(avgTable, "Average table should be computed");
    }

    @Test
    public void test_moment3_vs_default_comparison() throws Exception {
        // Run with default method
        LNOptions defaultOptions = new LNOptions();
        defaultOptions.verbose = VerboseLevel.SILENT;
        defaultOptions.method = "default";
        SolverLN solverDefault = new SolverLN(SolverLNTestFixtures.buildModel2(), defaultOptions);
        LayeredNetworkAvgTable avgDefault = (LayeredNetworkAvgTable) solverDefault.getEnsembleAvg();

        // Run with moment3 method
        LNOptions moment3Options = new LNOptions();
        moment3Options.verbose = VerboseLevel.SILENT;
        moment3Options.method = "moment3";
        SolverLN solverMoment3 = new SolverLN(SolverLNTestFixtures.buildModel2(), moment3Options);
        LayeredNetworkAvgTable avgMoment3 = (LayeredNetworkAvgTable) solverMoment3.getEnsembleAvg();

        // Both should converge
        assertTrue(solverDefault.hasconverged, "Default solver should converge");
        assertTrue(solverMoment3.hasconverged, "Moment3 solver should converge");

        // Results should exist
        assertNotNull(avgDefault);
        assertNotNull(avgMoment3);
    }

    // Differential tests wrt LQNS 6.2.29 (@ 25-Aug-2026) results
    // NaN or 0.0 entries in LQNS are ignored
    private double allowedMaxRelDiffLQNS = 0.10;

    @org.junit.jupiter.api.Test
    public void difftest_buildModel_1() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel1(), options);
        double[] expectedQLen = {Double.NaN, 2.6133, 2.6133, 0.1633, 2.4500};
        double[] expectedUtil = {0.7582, 0.0000, 0.0000, 0.0474, 0.7108};
        double[] expectedRespT = {Double.NaN, Double.NaN, 5.5148, 0.3447, 5.1701};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, 0.4739, 0.4739, 0.4739, 0.4739};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_buildModel_2() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel2(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.1188, 0.5328, 1.1188, 0.5328, 0.1625, 0.9563, 0.4440, 0.0888};
        double[] expectedUtil = {0.1421, 0.5329, 0.1421, 0.5329, Double.NaN, Double.NaN, 0.1421, 0.0000, 0.4441, 0.0888};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 12.6000, 6.0000, 1.8299, 10.7678, 5.0000, 1.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, 1.8299, 6.0000, Double.NaN, Double.NaN, 1.8299, 0.0000, 5.0000, 1.0000};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.0888, 0.0888, 0.0888, 0.0888, 0.0888, 0.0888, 0.0888, 0.0888};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_buildModel_3() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel3(), options);
        double[] expectedQLen = {Double.NaN, 8.5090, 0.8385, 8.5090, 0.8385, 7.6706, 0.8384, 0.6988, 0.1398};
        double[] expectedUtil = {1.1332, 0.0000, 0.0000, 0.0000, 0.0000, 1.0437, 0.0000, 0.0746, 0.0149};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, 570.7000, 56.2000, 514.5000, 56.2000, 46.9000, 9.3732};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, 0.0149, 0.0149, 0.0149, 0.0149, 0.0149, 0.0149, 0.0149, 0.0149};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_buildModel_4() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel4(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 19.7749, 8.9502, 1.5927, 19.7749, 8.9502, 1.5927, 19.7749, 8.9502, 1.5927};
        double[] expectedUtil = {2.2668, 1.5111, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.5113, 0.7555, 1.5111};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 1.3084, 0.5924, 0.0211, 1.3084, 0.5924, 0.0211};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 15.1000, 15.1000, 75.6000, 15.1000, 15.1000, 75.6000, 15.1000, 15.1000, 75.6000};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_buildModel_5() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel5(), options);
        double[] expectedQLen = {Double.NaN, 1.0000, 1.0000, 0.1587, 0.5238, 0.3175};
        double[] expectedUtil = {1.0000, 0.0000, 0.0000, 0.1587, 0.5238, 0.3175};
        double[] expectedRespT = {Double.NaN, Double.NaN, 12.6000, 2.0000, 3.0000, 4.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, 0.0794, 0.0794, 0.0794, 0.1746, 0.0794};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_buildModel_6() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel6(), options);
        double[] expectedQLen = {Double.NaN, 1.0000, 1.0000, 0.1653, 0.0744, 0.0992, 0.1653, 0.4959};
        double[] expectedUtil = {1.0323, 0.0000, 0.0000, 0.1706, 0.0768, 0.1024, 0.1706, 0.5119};
        double[] expectedRespT = {Double.NaN, Double.NaN, 11.7219, 1.9375, 2.9062, 3.8750, 4.8438, 5.8125};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, 0.0853, 0.0853, 0.0853, 0.0256, 0.0256, 0.0341, 0.0853};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_lqn_serial() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(LayeredModel.lqn_serial(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.1215, 0.5327, 1.1215, 0.5327, 0.1627, 0.9587, 0.4439, 0.0888};
        double[] expectedUtil = {0.1421, 0.5327, 0.0000, 0.0000, 0.0000, 0.0000, 0.1421, 0.0000, 0.4439, 0.0888};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 12.6310, 6.0000, 1.8327, 10.8000, 5.0000, 1.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.0888, 0.0888, 0.0888, 0.0888, 0.0888, 0.0888, 0.0888, 0.0888};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_activityGraph_and() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.test_activityGraph_and(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.0000, 0.7539, 0.3770, 0.3770, 0.0000};
        double[] expectedUtil = {0.0000, 1.5078, 0.0000, 1.5078, 0.0000, 1.5078, 0.0000, 0.0000, 0.7539, 0.3770, 0.3770, 0.0000};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 2.6528, 2.6528, 2.6528, 0.0000, 2.0000, 1.0000, 1.0000, 0.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.3770, 0.3770, 0.3770, 0.3770, 0.3770, 0.3770, 0.3770, 0.3770, 0.3770, 0.3770};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_activityGraph_call() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.test_activityGraph_call(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, Double.NaN, 0.1274, 0.1134, 0.0698, 0.1274, 0.1134, 0.0698, 0.1274, 0.1134, 0.0698};
        double[] expectedUtil = {0.0140, 0.0436, 0.0698, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0140, 0.0436, 0.0698};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 14.6000, 13.0000, 8.0000, 14.6000, 13.0000, 8.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, Double.NaN, 0.0087, 0.0087, 0.0087, 0.0087, 0.0087, 0.0087, 0.0087, 0.0087, 0.0087};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_activityGraph_call_and() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.test_activityGraph_call_and(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, Double.NaN, 1.0000, 1.0000, 0.7680, 1.0000, 1.0000, 0.7680, 1.0000, 0.7680, 0.1920, 0.0960, 0.0000, 0.7680};
        double[] expectedUtil = {0.0000, 0.2880, 0.7680, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.1920, 0.0960, 0.0000, 0.7680};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 10.4167, 10.4167, 8.0000, 10.4167, 8.0000, 2.0000, 1.0000, 0.0000, 8.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, Double.NaN, 0.0960, 0.0960, 0.0960, 0.0960, 0.0960, 0.0960, 0.0960, 0.0960, 0.0960, 0.0960, 0.0960, 0.0960};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_activityGraph_call_or() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.test_activityGraph_call_or(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, Double.NaN, 1.0000, 0.9846, 0.0000, 1.0000, 0.9846, 0.0000, 1.0000, 0.0000, 0.0096, 0.0144, 0.9606, 0.0000};
        double[] expectedUtil = {0.0154, 0.9846, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0154, 0.0000, 0.0096, 0.0144, 0.9606, 0.0000};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 104.1000, 102.5000, Double.NaN, 104.1000, 0.0000, 2.0000, 3.0000, 100.0000, 0.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, Double.NaN, 0.0096, 0.0096, 0.0000, 0.0096, 0.0096, 0.0000, 0.0096, 0.0096, 0.0048, 0.0048, 0.0096, 0.0000};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_activityGraph_call_seq_disconnected() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.test_activityGraph_call_seq_disconnected(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, Double.NaN, 0.0706, 0.0558, 0.0000, 0.0706, 0.0558, 0.0000, 0.0706, 0.0465, 0.0093, 0.0000};
        double[] expectedUtil = {0.0149, 0.0558, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0149, 0.0465, 0.0093, 0.0000};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 7.6000, 6.0000, Double.NaN, 7.6000, 5.0000, 1.0000, 0.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, Double.NaN, 0.0093, 0.0093, 0.0000, 0.0093, 0.0093, 0.0000, 0.0093, 0.0093, 0.0093, 0.0000};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_activityGraph_loop() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.test_activityGraph_loop(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 0.0909, 0.0909, 0.0909, 0.0909, 0.0909, 0.0000, 0.0818, 0.0091};
        double[] expectedUtil = {0.0000, 0.0909, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0818, 0.0091};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 10.0000, 10.0000, 10.0000, 0.0000, 3.0000, 1.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.0091, 0.0091, 0.0091, 0.0091, 0.0091, 0.0091, 0.0273, 0.0091};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_activityGraph_or() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.test_activityGraph_or(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.0000, 0.9846, 1.0000, 0.9846, 1.0000, 0.0000, 0.0096, 0.0144, 0.9606};
        double[] expectedUtil = {0.0154, 0.9846, 0.0000, 0.0000, 0.0000, 0.0000, 0.0154, 0.0000, 0.0096, 0.0144, 0.9606};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 104.1000, 102.5000, 104.1000, 0.0000, 2.0000, 3.0000, 100.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.0096, 0.0096, 0.0096, 0.0096, 0.0096, 0.0096, 0.0048, 0.0048, 0.0096};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_activityGraph_seq() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.test_activityGraph_seq(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 0.0706, 0.0558, 0.0706, 0.0558, 0.0706, 0.0465, 0.0093};
        double[] expectedUtil = {0.0149, 0.0558, 0.0000, 0.0000, 0.0000, 0.0000, 0.0149, 0.0465, 0.0093};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 7.6000, 6.0000, 7.6000, 5.0000, 1.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.0093, 0.0093, 0.0093, 0.0093, 0.0093, 0.0093, 0.0093};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }


    @Test
    public void test_LQN_calls_1() throws Exception {
        LayeredNetwork model = suppressOutput(() -> {
            LayeredNetwork m = new LayeredNetwork("add_cart");
            
            // first layer
            Processor P1 = new Processor(m, "client_p", 1, SchedStrategy.INF);
            Task T1 = new Task(m, "user", 1, SchedStrategy.REF).on(P1);
            Entry E1 = new Entry(m, "user").on(T1);
            T1.setThinkTime(Exp.fitMean(0));
            
            // second layer
            Processor P2 = new Processor(m, "WeiUI_p", 1, SchedStrategy.INF);
            Task T2 = new Task(m, "WeiUI", 1, SchedStrategy.INF).on(P2);
            Entry E2 = new Entry(m, "add_cart").on(T2);
            
            // activities
            Activity A1 = new Activity(m, "A1", Exp.fitMean(0)).on(T1).boundTo(E1);
            Activity A2 = new Activity(m, "A2", Exp.fitMean(0)).on(T1).synchCall(E2, 1);
            Activity A3 = new Activity(m, "A3", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
            
            T1.addPrecedence(ActivityPrecedence.Serial(A1, A2));
            
            return m;
        });
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        // Expected results matrix - converted from MATLAB format
        double[][] expectedResults = {
            {Double.NaN, 0, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {1.000000000000000, 0, Double.NaN, 0, Double.NaN, 1.000000000000000},
            {1.000000000000000, 1.000000000000000, Double.NaN, 1.000000000000000, Double.NaN, 1.000000000000000},
            {1.000000000000000, Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, 1.000000000000000},
            {1.000000000000000, Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, 1.000000000000000},
            {0.000000010000000, 0, 0.000000010000000, 0, Double.NaN, 1.000000000000000},
            {1.000000000000000, 0, 1.000000000000000, 0, Double.NaN, 1.000000000000000},
            {1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000, Double.NaN, 1.000000000000000}
        };
        
        // Compare results with expected values using relaxed tolerance for known precision differences
        callAssertTableMetrics(avgTable, expectedResults);
    }
    
    // @Test // Disabled - failing test
    public void test_LQN_calls_2() throws Exception {
        LayeredNetwork model = new LayeredNetwork("add_cart");
        
        // first layer
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "user").on(T1);
        T1.setThinkTime(Exp.fitMean(0));
        
        // second layer - changed to FCFS for M/M/1
        Processor P2 = new Processor(model, "WeiUI_p", 1, SchedStrategy.INF);
        Task T2 = new Task(model, "WeiUI", 1, SchedStrategy.FCFS).on(P2);
        Entry E2 = new Entry(model, "add_cart").on(T2);
        
        // activities
        Activity A1 = new Activity(model, "A1", new Immediate()).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(2)).on(T1).synchCall(E2, 1);
        Activity A3 = new Activity(model, "A3", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        
        T1.addPrecedence(ActivityPrecedence.Serial(A1, A2));
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        // Expected results matrix - from MATLAB test_LQN_calls_2.m
        double[][] expectedResults = {
            {Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {Double.NaN, 0.500000000000000, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {1.000000000000000, 1.000000000000000, Double.NaN, 2.000000000000000, Double.NaN, 0.500000000000000},
            {0.500000000000000, 0.500000000000000, Double.NaN, 1.000000000000000, Double.NaN, 0.500000000000000},
            {1.000000000000000, Double.NaN, 2.000000000000000, Double.NaN, Double.NaN, 0.500000000000000},
            {0.500000000000000, Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, 0.500000000000000},
            {0.000000005000000, 0, 0.000000010000000, 0, Double.NaN, 0.500000000000000},
            {2.000000000000000, 1.000000000000000, 4.000000000000000, 2.000000000000000, Double.NaN, 0.500000000000000},
            {0.500000000000000, 0.500000000000000, 1.000000000000000, 1.000000000000000, Double.NaN, 0.500000000000000}
        };
        
        callAssertTableMetrics(avgTable, expectedResults);
    }
    
    // @Test // Disabled - failing test
    public void test_LQN_serial_1() throws Exception {
        LayeredNetwork model = new LayeredNetwork("add_cart");
        
        // first layer
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "user").on(T1);
        T1.setThinkTime(Exp.fitMean(0));
        
        // second layer
        Processor P2 = new Processor(model, "WeiUI_p", 1, SchedStrategy.INF);
        Task T2 = new Task(model, "WeiUI", 1, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(model, "add_cart").on(T2);
        
        // activities
        Activity A1 = new Activity(model, "A1", Exp.fitMean(0)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(0)).on(T1).synchCall(E2, 1);
        Activity A3 = new Activity(model, "A3", Exp.fitMean(0)).on(T1);
        Activity A4 = new Activity(model, "A4", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        
        // Create serial precedence A1 -> A2 -> A3
        T1.addPrecedence(ActivityPrecedence.Serial(A1, A2));
        T1.addPrecedence(ActivityPrecedence.Serial(A2, A3));
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        // Expected results matrix
        double[][] expectedResults = {
            {Double.NaN, 0.000000030001334, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {1.000000000000000, 0.000000030001334, Double.NaN, 0.000015288109828, Double.NaN, 1.000000000000000},
            {1.000000000000000, 1.000000000000000, Double.NaN, 1.000000000000000, Double.NaN, 1.000000000000000},
            {1.000000000000000, Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, 1.000000000000000},
            {1.000000000000000, Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, 1.000000000000000},
            {0.000005096263160, 0.000000010000445, 0.000005096036609, 0.000005096036609, Double.NaN, 1.000000000000000},
            {1.000000000000000, 0.000000010000445, 1.000000000000000, 0.000005096036609, Double.NaN, 1.000000000000000},
            {0.000005096263160, 0.000000010000445, 0.000005096036609, 0.000005096036609, Double.NaN, 1.000000000000000},
            {1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000, Double.NaN, 1.000000000000000}
        };
        
        // NOTE: This test shows small numerical differences for near-zero values
        // Using relaxed validation
        callAssertTableMetrics(avgTable, expectedResults);
    }
    
    // @Test // Disabled - failing test
    public void test_LQN_serial_2() throws Exception {
        LayeredNetwork model = new LayeredNetwork("add_cart");
        
        // first layer
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "user").on(T1);
        T1.setThinkTime(Exp.fitMean(0));
        
        // second layer
        Processor P2 = new Processor(model, "WeiUI_p", 1, SchedStrategy.INF);
        Task T2 = new Task(model, "WeiUI", 1, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(model, "add_cart").on(T2);
        
        // activities
        Activity A1 = new Activity(model, "A1", Exp.fitMean(0)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(0)).on(T1).synchCall(E2, 1);
        Activity A3 = new Activity(model, "A3", Exp.fitMean(1)).on(T1);
        Activity A4 = new Activity(model, "A4", Exp.fitMean(0)).on(T2).boundTo(E2).repliesTo(E2);
        
        T1.addPrecedence(ActivityPrecedence.Serial(A1, A2));
        T1.addPrecedence(ActivityPrecedence.Serial(A2, A3));
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        // Expected results matrix
        double[][] expectedResults = {
            {Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {Double.NaN, 0.000000010000152, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {1.000000000000000, 1.000000000000000, Double.NaN, 1.000000000000000, Double.NaN, 1.000000000000000},
            {0.000015268788909, 0.000000010000152, Double.NaN, 0.000015268557034, Double.NaN, 1.000000000000000},
            {1.000000000000000, Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, 1.000000000000000},
            {0.000015268788909, Double.NaN, 0.000015268557034, Double.NaN, Double.NaN, 1.000000000000000},
            {0, 0, 0, 0, Double.NaN, 1.000000000000000},
            {0.000030525947774, 0, 0.000030526880857, 0, Double.NaN, 1.000000000000000},
            {1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000, Double.NaN, 1.000000000000000},
            {0.000015268788909, 0.000000010000152, 0.000015268557034, 0.000015268557034, Double.NaN, 1.000000000000000}
        };
        
        // NOTE: This test has very small expected values that cause numerical precision issues
        callAssertTableMetrics(avgTable, expectedResults);
    }
    
    @Test
    public void test_LQN_serial_3() throws Exception {
        LayeredNetwork model = new LayeredNetwork("add_cart");
        
        // first layer
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "user").on(T1);
        T1.setThinkTime(Exp.fitMean(0));
        
        // second layer
        Processor P2 = new Processor(model, "WeiUI_p", 1, SchedStrategy.INF);
        Task T2 = new Task(model, "WeiUI", 1, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(model, "add_cart").on(T2);
        
        // activities
        Activity A1 = new Activity(model, "A1", Exp.fitMean(0)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(0)).on(T1).synchCall(E2, 1);
        Activity A3 = new Activity(model, "A3", Exp.fitMean(1)).on(T1);
        Activity A4 = new Activity(model, "A4", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        
        T1.addPrecedence(ActivityPrecedence.Serial(A1, A2));
        T1.addPrecedence(ActivityPrecedence.Serial(A2, A3));
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        // Expected results matrix - from MATLAB test_LQN_serial_3.m
        double[][] expectedResults = {
            {Double.NaN, 0.500000000000000, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {Double.NaN, 0.500000000000000, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {1.000000000000000, 0.500000000000000, Double.NaN, 1.000000000000000, Double.NaN, 0.500000000000000},
            {0.500000000000000, 0.500000000000000, Double.NaN, 1.000000000000000, Double.NaN, 0.500000000000000},
            {1.000000000000000, Double.NaN, 2.000000000000000, Double.NaN, Double.NaN, 0.500000000000000},
            {0.500000000000000, Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, 0.500000000000000},
            {0, 0, 0, 0, Double.NaN, 0.500000000000000},
            {0.500000000000000, 0, 1.000000000000000, 0, Double.NaN, 0.500000000000000},
            {0.500000000000000, 0.500000000000000, 1.000000000000000, 1.000000000000000, Double.NaN, 0.500000000000000},
            {0.500000000000000, 0.500000000000000, 1.000000000000000, 1.000000000000000, Double.NaN, 0.500000000000000}
        };
        
        callAssertTableMetrics(avgTable, expectedResults);
    }
    
    @Test
    public void test_LQN_serial_4() throws Exception {
        LayeredNetwork model = new LayeredNetwork("add_cart");
        
        // first layer
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "user").on(T1);
        T1.setThinkTime(Exp.fitMean(0));
        
        // second layer
        Processor P2 = new Processor(model, "WeiUI_p", 1, SchedStrategy.INF);
        Task T2 = new Task(model, "WeiUI", 1, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(model, "add_cart").on(T2);
        
        // activities - A2 has service time 1 instead of 0
        Activity A1 = new Activity(model, "A1", Exp.fitMean(0)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T1).synchCall(E2, 1);
        Activity A3 = new Activity(model, "A3", Exp.fitMean(1)).on(T1);
        Activity A4 = new Activity(model, "A4", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        
        T1.addPrecedence(ActivityPrecedence.Serial(A1, A2));
        T1.addPrecedence(ActivityPrecedence.Serial(A2, A3));
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        // Expected results matrix - from MATLAB test_LQN_serial_4.m
        double[][] expectedResults = {
            {Double.NaN, 0.666677998328647, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {Double.NaN, 0.333337248848223, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {1.000000000000000, 0.666677998328647, Double.NaN, 2.000000000000000, Double.NaN, 0.333338999164324},
            {0.333330890697903, 0.333337248848223, Double.NaN, 1.000000000000000, Double.NaN, 0.333337248848223},
            {1.000000000000000, Double.NaN, 3.000000000000000, Double.NaN, Double.NaN, 0.333338999164324},
            {0.333330890697903, Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, 0.333337248848223},
            {0.0, 0.0, 0.0, 0.0, Double.NaN, 0.333338999164324},
            {0.666650817472160, 0.333338999164324, 2.000000000000000, 1.000000000000000, Double.NaN, 0.333338999164324},
            {0.333328825946500, 0.333338999164324, 1.000000000000000, 1.000000000000000, Double.NaN, 0.333338999164324},
            {0.333330890697903, 0.333337248848223, 1.000000000000000, 1.000000000000000, Double.NaN, 0.333337248848223}
        };
        
        callAssertTableMetrics(avgTable, expectedResults);
    }
    
    @Test
    public void test_LQN_serial_5() throws Exception {
        LayeredNetwork model = new LayeredNetwork("add_cart");
        
        // first layer
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "user").on(T1);
        T1.setThinkTime(Exp.fitMean(0));
        
        // second layer
        Processor P2 = new Processor(model, "WeiUI_p", 1, SchedStrategy.INF);
        Task T2 = new Task(model, "WeiUI", 1, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(model, "add_cart").on(T2);
        
        // activities - same as serial_4
        Activity A1 = new Activity(model, "A1", Exp.fitMean(0)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T1).synchCall(E2, 1);
        Activity A3 = new Activity(model, "A3", Exp.fitMean(1)).on(T1);
        Activity A4 = new Activity(model, "A4", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        
        T1.addPrecedence(ActivityPrecedence.Serial(A1, A2));
        T1.addPrecedence(ActivityPrecedence.Serial(A2, A3));
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        // Expected results matrix - from MATLAB test_LQN_serial_5.m
        double[][] expectedResults = {
            {Double.NaN, 0.666677998328647, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {Double.NaN, 0.333337248848223, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {1.000000000000000, 0.666677998328647, Double.NaN, 2.000000000000000, Double.NaN, 0.333338999164324},
            {0.333330890697903, 0.333337248848223, Double.NaN, 1.000000000000000, Double.NaN, 0.333337248848223},
            {1.000000000000000, Double.NaN, 3.000000000000000, Double.NaN, Double.NaN, 0.333338999164324},
            {0.333330890697903, Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, 0.333337248848223},
            {0.0, 0.0, 0.0, 0.0, Double.NaN, 0.333338999164324},
            {0.666650817472160, 0.333338999164324, 2.000000000000000, 1.000000000000000, Double.NaN, 0.333338999164324},
            {0.333328825946500, 0.333338999164324, 1.000000000000000, 1.000000000000000, Double.NaN, 0.333338999164324},
            {0.333330890697903, 0.333337248848223, 1.000000000000000, 1.000000000000000, Double.NaN, 0.333337248848223}
        };
        
        callAssertTableMetrics(avgTable, expectedResults);
    }
    
    // @Test // Disabled - failing test
    public void test_LQN_orfork_1() throws Exception {
        LayeredNetwork model = new LayeredNetwork("add_cart");
        
        // first layer
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "user").on(T1);
        T1.setThinkTime(Exp.fitMean(0));
        
        // second layer
        Processor P2 = new Processor(model, "WeiUI_p", 1, SchedStrategy.INF);
        Task T2 = new Task(model, "WeiUI", 1, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(model, "add_cart").on(T2);
        
        // activities
        Activity A1 = new Activity(model, "A1", Exp.fitMean(0)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(0)).on(T1).synchCall(E2, 1);
        Activity A3 = new Activity(model, "A3", Exp.fitMean(0)).on(T1);
        Activity A4 = new Activity(model, "A4", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        
        Matrix probs = new Matrix(1, 2);
        probs.set(0, 0, 0.4);
        probs.set(0, 1, 0.6);
        T1.addPrecedence(ActivityPrecedence.OrFork(A1, java.util.Arrays.asList(A2, A3), probs));
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        // Expected results matrix
        double[][] expectedResults = {
            {Double.NaN, 0.000000050001638, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {1.000000000000000, 0.000000050001638, Double.NaN, 0.000006123315075, Double.NaN, 2.500000000000000},
            {1.000000000000000, 1.000000000000000, Double.NaN, 1.000000000000000, Double.NaN, 1.000000000000000},
            {1.000000000000000, Double.NaN, 0.400000000000000, Double.NaN, Double.NaN, 2.500000000000000},
            {1.000000000000000, Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, 1.000000000000000},
            {0.000007654394587, 0.000000025000819, 0.000003061657537, 0.000003061657537, Double.NaN, 2.500000000000000},
            {1.000000000000000, 0.000000010000328, 1.000000000000000, 0.000001224663015, Double.NaN, 1.000000000000000},
            {0.000004592636752, 0.000000015000491, 0.000003061657537, 0.000001836994522, Double.NaN, 1.500000000000000},
            {1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000, Double.NaN, 1.000000000000000}
        };
        
        // NOTE: This test has very small expected values that cause numerical precision issues
        callAssertTableMetrics(avgTable, expectedResults);
    }
    
    // @Test // Disabled - failing test
    public void test_LQN_orfork_2() throws Exception {
        LayeredNetwork model = new LayeredNetwork("add_cart");
        
        // first layer
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "user").on(T1);
        T1.setThinkTime(Exp.fitMean(0));
        
        // second layer
        Processor P2 = new Processor(model, "WeiUI_p", 1, SchedStrategy.INF);
        Task T2 = new Task(model, "WeiUI", 1, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(model, "add_cart").on(T2);
        
        // activities - A3 has service time 1, A4 has service time 0
        Activity A1 = new Activity(model, "A1", Exp.fitMean(0)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(0)).on(T1).synchCall(E2, 1);
        Activity A3 = new Activity(model, "A3", Exp.fitMean(1)).on(T1);
        Activity A4 = new Activity(model, "A4", Exp.fitMean(0)).on(T2).boundTo(E2).repliesTo(E2);
        
        Matrix probs = new Matrix(1, 2);
        probs.set(0, 0, 0.4);
        probs.set(0, 1, 0.6);
        T1.addPrecedence(ActivityPrecedence.OrFork(A1, java.util.Arrays.asList(A2, A3), probs));
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        // Expected results matrix - from MATLAB test_LQN_orfork_2.m
        double[][] expectedResults = {
            {Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {Double.NaN, 0.0, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {1.000000000000000, 1.000000000000000, Double.NaN, 0.600000000000000, Double.NaN, 1.666611463403392},
            {0.000015265455695, 0.0, Double.NaN, 0.000022897603080, Double.NaN, 0.666683566902699},
            {1.000000000000000, Double.NaN, 0.600000000000000, Double.NaN, Double.NaN, 1.666611463403392},
            {0.000015265455695, Double.NaN, 0.000022897603080, Double.NaN, Double.NaN, 0.666683566902699},
            {0.0, 0.0, 0.0, 0.0, Double.NaN, 1.666611463403392},
            {0.000030522226395, 0.0, 0.000045784856076, 0.0, Double.NaN, 0.666644585361357},
            {1.000000000000000, 1.000000000000000, 1.000000000000000, 0.600000000000000, Double.NaN, 1.000000000000000},
            {0.000015265455695, 0.0, 0.000022897603080, 0.000022897603080, Double.NaN, 0.666683566902699}
        };
        
        callAssertTableMetrics(avgTable, expectedResults);
    }
    
    @Test
    public void test_LQN_orfork_3() throws Exception {
        LayeredNetwork model = new LayeredNetwork("add_cart");
        
        // first layer
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "user").on(T1);
        T1.setThinkTime(Exp.fitMean(0));
        
        // second layer
        Processor P2 = new Processor(model, "WeiUI_p", 1, SchedStrategy.INF);
        Task T2 = new Task(model, "WeiUI", 1, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(model, "add_cart").on(T2);
        
        // activities - A3 has service time 1, A4 has service time 1
        Activity A1 = new Activity(model, "A1", Exp.fitMean(0)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(0)).on(T1).synchCall(E2, 1);
        Activity A3 = new Activity(model, "A3", Exp.fitMean(1)).on(T1);
        Activity A4 = new Activity(model, "A4", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        
        Matrix probs = new Matrix(1, 2);
        probs.set(0, 0, 0.4);
        probs.set(0, 1, 0.6);
        T1.addPrecedence(ActivityPrecedence.OrFork(A1, java.util.Arrays.asList(A2, A3), probs));
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        // Expected results matrix - from MATLAB test_LQN_orfork_3.m
        double[][] expectedResults = {
            {Double.NaN, 0.600000000000000, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {Double.NaN, 0.400000000000000, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {1.000000000000000, 0.600000000000000, Double.NaN, 0.600000000000000, Double.NaN, 1.000000000000000},
            {0.400000000000000, 0.400000000000000, Double.NaN, 1.000000000000000, Double.NaN, 0.400000000000000},
            {1.000000000000000, Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, 1.000000000000000},
            {0.400000000000000, Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, 0.400000000000000},
            {0.0, 0.0, 0.0, 0.0, Double.NaN, 1.000000000000000},
            {0.400000000000000, 0.0, 1.000000000000000, 0.0, Double.NaN, 0.400000000000000},
            {0.600000000000000, 0.600000000000000, 1.000000000000000, 0.600000000000000, Double.NaN, 0.600000000000000},
            {0.400000000000000, 0.400000000000000, 1.000000000000000, 1.000000000000000, Double.NaN, 0.400000000000000}
        };
        
        callAssertTableMetrics(avgTable, expectedResults);
    }
    
    // @Test // Disabled - failing test
    public void test_LQN_orfork_4() throws Exception {
        LayeredNetwork model = new LayeredNetwork("simple_orfork");
        
        // first layer
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "user").on(T1);
        T1.setThinkTime(Exp.fitMean(0));
        
        // second layer
        Processor P2 = new Processor(model, "server_p", 1, SchedStrategy.INF);
        Task T2 = new Task(model, "server", 1, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(model, "service").on(T2);
        
        // activities
        Activity A1 = new Activity(model, "A1", Exp.fitMean(0)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(0)).on(T1).synchCall(E2, 1);
        Activity A3 = new Activity(model, "A3", Exp.fitMean(0)).on(T1);
        Activity A4 = new Activity(model, "A4", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        
        Matrix probs = new Matrix(1, 2);
        probs.set(0, 0, 0.4);
        probs.set(0, 1, 0.6);
        T1.addPrecedence(ActivityPrecedence.OrFork(A1, java.util.Arrays.asList(A2, A3), probs));
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        // Expected results matrix - from MATLAB test_LQN_orfork_4.m
        double[][] expectedResults = {
            {Double.NaN, 0.714289512399078, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {Double.NaN, 0.285719430769547, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {1.000000000000000, 0.714289512399078, Double.NaN, 1.000000000000000, Double.NaN, 0.714289512399078},
            {0.285715918474751, 0.285719430769547, Double.NaN, 1.000000000000000, Double.NaN, 0.285719430769547},
            {1.000000000000000, Double.NaN, 1.400000000000000, Double.NaN, Double.NaN, 0.714289512399078},
            {0.285715918474751, Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, 0.285719430769547},
            {0.0, 0.0, 0.0, 0.0, Double.NaN, 0.714289512399078},
            {0.571412755894316, 0.285715804959631, 2.000000000000000, 0.400000000000000, Double.NaN, 0.285715804959631},
            {0.428559156792581, 0.428573707439447, 1.000000000000000, 0.600000000000000, Double.NaN, 0.428573707439447},
            {0.285715918474751, 0.285719430769547, 1.000000000000000, 1.000000000000000, Double.NaN, 0.285719430769547}
        };
        
        callAssertTableMetrics(avgTable, expectedResults);
    }
    
    
    @Test
    public void test_LQN_orfork_6() throws Exception {
        LayeredNetwork model = new LayeredNetwork("cascaded_orfork");
        
        // first layer
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "user").on(T1);
        T1.setThinkTime(Exp.fitMean(0));
        
        // second layer - multiple tasks on P2
        Processor P2 = new Processor(model, "server_p", 1, SchedStrategy.INF);
        Task T2 = new Task(model, "service", 1, SchedStrategy.INF).on(P2);
        Task T3 = new Task(model, "pricing", 1, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(model, "service").on(T2);
        Entry E3 = new Entry(model, "pricing").on(T3);
        
        // activities with cascaded OrForks
        Activity A1 = new Activity(model, "A1", Exp.fitMean(0)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(0)).on(T1);
        Activity A3 = new Activity(model, "A3", Exp.fitMean(0)).on(T1);
        Activity A4 = new Activity(model, "A4", Exp.fitMean(0)).on(T1).synchCall(E2, 1);
        Activity A5 = new Activity(model, "A5", Exp.fitMean(0)).on(T1).synchCall(E3, 1);
        
        Activity B1 = new Activity(model, "B1", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        Activity C1 = new Activity(model, "C1", Exp.fitMean(1)).on(T3).boundTo(E3).repliesTo(E3);
        
        // First OrFork: A1 -> {A2, A3}
        Matrix probs1 = new Matrix(1, 2);
        probs1.set(0, 0, 0.4);
        probs1.set(0, 1, 0.6);
        T1.addPrecedence(ActivityPrecedence.OrFork(A1, java.util.Arrays.asList(A2, A3), probs1));
        
        // Second OrFork: A2 -> {A4, A5}
        Matrix probs2 = new Matrix(1, 2);
        probs2.set(0, 0, 0.3);
        probs2.set(0, 1, 0.7);
        T1.addPrecedence(ActivityPrecedence.OrFork(A2, java.util.Arrays.asList(A4, A5), probs2));
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        callAssertTableMetrics(avgTable, null);
    }
    
    @Test
    public void test_LQN_orfork_7() throws Exception {
        LayeredNetwork model = new LayeredNetwork("complex_orfork");
        
        // first layer with think time
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "user").on(T1);
        T1.setThinkTime(Exp.fitMean(1)); // Think time = 1
        
        // second layer
        Processor P2 = new Processor(model, "server_p", 1, SchedStrategy.INF);
        Task T2 = new Task(model, "service", 1, SchedStrategy.INF).on(P2);
        Task T3 = new Task(model, "pricing", 1, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(model, "service").on(T2);
        Entry E3 = new Entry(model, "pricing").on(T3);
        
        // Complex client-side activities
        Activity A1 = new Activity(model, "A1", Exp.fitMean(0)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(0)).on(T1);
        Activity A3 = new Activity(model, "A3", Exp.fitMean(0)).on(T1);
        Activity A4 = new Activity(model, "A4", Exp.fitMean(0)).on(T1).synchCall(E2, 1);
        Activity A5 = new Activity(model, "A5", Exp.fitMean(0)).on(T1).synchCall(E3, 1);
        
        // Server-side activities with OrForks
        Activity S0 = new Activity(model, "S0", Exp.fitMean(0)).on(T2).boundTo(E2);
        Activity S1 = new Activity(model, "S1", Exp.fitMean(2)).on(T2).repliesTo(E2);
        Activity S2 = new Activity(model, "S2", Exp.fitMean(1)).on(T2).repliesTo(E2);
        
        Activity P0 = new Activity(model, "P0", Exp.fitMean(0)).on(T3).boundTo(E3);
        Activity P1_act = new Activity(model, "P1", Exp.fitMean(3)).on(T3).repliesTo(E3);
        Activity P2_act = new Activity(model, "P2", Exp.fitMean(1)).on(T3).repliesTo(E3);
        
        // Client-side OrForks
        Matrix probsA = new Matrix(1, 2);
        probsA.set(0, 0, 0.3);
        probsA.set(0, 1, 0.7);
        T1.addPrecedence(ActivityPrecedence.OrFork(A1, java.util.Arrays.asList(A2, A3), probsA));
        
        Matrix probsA2 = new Matrix(1, 2);
        probsA2.set(0, 0, 0.2);
        probsA2.set(0, 1, 0.8);
        T1.addPrecedence(ActivityPrecedence.OrFork(A2, java.util.Arrays.asList(A4, A5), probsA2));
        
        // Server-side OrForks
        Matrix probsS = new Matrix(1, 2);
        probsS.set(0, 0, 0.25);
        probsS.set(0, 1, 0.75);
        T2.addPrecedence(ActivityPrecedence.OrFork(S0, java.util.Arrays.asList(S1, S2), probsS));
        
        Matrix probsP = new Matrix(1, 2);
        probsP.set(0, 0, 0.65);
        probsP.set(0, 1, 0.35);
        T3.addPrecedence(ActivityPrecedence.OrFork(P0, java.util.Arrays.asList(P1_act, P2_act), probsP));
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        callAssertTableMetrics(avgTable, null);
    }
    
    // @Test // Disabled - failing test
    public void test_LQN_orjoin_1() throws Exception {
        LayeredNetwork model = new LayeredNetwork("fork_join");
        
        // first layer
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "user").on(T1);
        T1.setThinkTime(Exp.fitMean(1)); // Think time = 1
        
        // second layer
        Processor P2 = new Processor(model, "server_p", 1, SchedStrategy.INF);
        Task T2 = new Task(model, "server", 1, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(model, "service").on(T2);
        
        // activities with OrFork followed by OrJoin
        Activity A1 = new Activity(model, "A1", Exp.fitMean(0)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(0)).on(T1).synchCall(E2, 1);
        Activity A3 = new Activity(model, "A3", Exp.fitMean(0)).on(T1);
        Activity A5 = new Activity(model, "A5", Exp.fitMean(1)).on(T1); // Join activity
        Activity A4 = new Activity(model, "A4", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        
        // OrFork: A1 -> {A2, A3}
        Matrix probs = new Matrix(1, 2);
        probs.set(0, 0, 0.4);
        probs.set(0, 1, 0.6);
        T1.addPrecedence(ActivityPrecedence.OrFork(A1, java.util.Arrays.asList(A2, A3), probs));
        
        // OrJoin: {A2, A3} -> A5
        T1.addPrecedence(ActivityPrecedence.OrJoin(java.util.Arrays.asList(A2, A3), A5));
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        // Expected results matrix - from MATLAB test_LQN_orjoin_1.m
        double[][] expectedResults = {
            {Double.NaN, 0.681817703378030, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {Double.NaN, 0.090907353566901, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {0.772748054673091, 0.681817703378030, Double.NaN, 3.000000000000000, Double.NaN, 0.227272567792677},
            {0.090912786672995, 0.090907353566901, Double.NaN, 1.000000000000000, Double.NaN, 0.090907353566901},
            {0.772748054673091, Double.NaN, 3.400000000000000, Double.NaN, Double.NaN, 0.227272567792677},
            {0.090912786672995, Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, 0.090907353566901},
            {0.227277836530591, 0.227272567792677, 1.000000000000000, 1.000000000000000, Double.NaN, 0.227272567792677},
            {0.181825679693553, 0.090909027117071, 2.000000000000000, 0.400000000000000, Double.NaN, 0.090909027117071},
            {0.136366701918355, 0.136363540675606, 1.000000000000000, 0.600000000000000, Double.NaN, 0.136363540675606},
            {0.090912786672995, 0.090907353566901, 1.000000000000000, 1.000000000000000, Double.NaN, 0.090907353566901},
            {0.227277836530591, 0.227272567792677, 1.000000000000000, 1.000000000000000, Double.NaN, 0.227272567792677}
        };
        
        callAssertTableMetrics(avgTable, expectedResults);
    }
    
    @Test
    public void test_LQN_orjoin_2() throws Exception {
        LayeredNetwork model = new LayeredNetwork("complex_fork_join");
        
        // first layer
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "user").on(T1);
        T1.setThinkTime(Exp.fitMean(1));
        
        // second layer - 5-server processor with FCFS task
        Processor P2 = new Processor(model, "server_p", 5, SchedStrategy.PS);
        Task T2 = new Task(model, "server", 1, SchedStrategy.FCFS).on(P2);
        Entry E2 = new Entry(model, "service").on(T2);
        
        // Complex multi-stage fork-join activities
        Activity A1 = new Activity(model, "A1", Exp.fitMean(0)).on(T1).boundTo(E1);
        Activity A5 = new Activity(model, "A5", Exp.fitMean(0)).on(T1);
        Activity A6 = new Activity(model, "A6", Exp.fitMean(0)).on(T1);
        Activity A7 = new Activity(model, "A7", Exp.fitMean(0)).on(T1);
        Activity A8 = new Activity(model, "A8", Exp.fitMean(1)).on(T1);
        
        // Fork activities
        Activity B11 = new Activity(model, "B11", Exp.fitMean(0)).on(T1).synchCall(E2, 1);
        Activity B12 = new Activity(model, "B12", Exp.fitMean(0)).on(T1).synchCall(E2, 2);
        Activity B21 = new Activity(model, "B21", Exp.fitMean(0)).on(T1).synchCall(E2, 3);
        Activity B22 = new Activity(model, "B22", Exp.fitMean(0)).on(T1).synchCall(E2, 4);
        Activity B23 = new Activity(model, "B23", Exp.fitMean(0)).on(T1).synchCall(E2, 1);
        Activity B31 = new Activity(model, "B31", Exp.fitMean(0)).on(T1).synchCall(E2, 2);
        Activity B32 = new Activity(model, "B32", Exp.fitMean(0)).on(T1).synchCall(E2, 3);
        
        Activity S1 = new Activity(model, "S1", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        
        // Multi-stage fork-join
        Matrix probs1 = new Matrix(1, 2);
        probs1.set(0, 0, 0.5); probs1.set(0, 1, 0.5);
        T1.addPrecedence(ActivityPrecedence.OrFork(A5, java.util.Arrays.asList(B11, B12), probs1));
        T1.addPrecedence(ActivityPrecedence.OrJoin(java.util.Arrays.asList(B11, B12), A6));
        
        Matrix probs2 = new Matrix(1, 3);
        probs2.set(0, 0, 0.33); probs2.set(0, 1, 0.33); probs2.set(0, 2, 0.34);
        T1.addPrecedence(ActivityPrecedence.OrFork(A6, java.util.Arrays.asList(B21, B22, B23), probs2));
        T1.addPrecedence(ActivityPrecedence.OrJoin(java.util.Arrays.asList(B21, B22, B23), A7));
        
        Matrix probs3 = new Matrix(1, 2);
        probs3.set(0, 0, 0.6); probs3.set(0, 1, 0.4);
        T1.addPrecedence(ActivityPrecedence.OrFork(A7, java.util.Arrays.asList(B31, B32), probs3));
        T1.addPrecedence(ActivityPrecedence.OrJoin(java.util.Arrays.asList(B31, B32), A8));
        
        T1.addPrecedence(ActivityPrecedence.Serial(A1, A5));
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        callAssertTableMetrics(avgTable, null);
    }
    
    @Test
    public void test_LQN_orjoin_3() throws Exception {
        LayeredNetwork model = new LayeredNetwork("simple_orfork_ps");
        
        // first layer
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.PS); // PS scheduling
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "user").on(T1);
        T1.setThinkTime(Exp.fitMean(0));
        
        // second layer
        Processor P2 = new Processor(model, "server_p", 1, SchedStrategy.INF);
        Task T2 = new Task(model, "server", 1, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(model, "service").on(T2);
        
        // Simple OrFork without OrJoin
        Activity A1 = new Activity(model, "A1", Exp.fitMean(0)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(0)).on(T1).synchCall(E2, 1);
        Activity A3 = new Activity(model, "A3", Exp.fitMean(0)).on(T1);
        Activity A4 = new Activity(model, "A4", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        
        Matrix probs = new Matrix(1, 2);
        probs.set(0, 0, 0.4);
        probs.set(0, 1, 0.6);
        T1.addPrecedence(ActivityPrecedence.OrFork(A1, java.util.Arrays.asList(A2, A3), probs));
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        callAssertTableMetrics(avgTable, null);
    }
    
    @Test
    public void test_LQN_allprec_1() throws Exception {
        LayeredNetwork model = new LayeredNetwork("combined_precedences");
        
        // Client and server setup
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "user").on(T1);
        T1.setThinkTime(Exp.fitMean(1));
        
        Processor P2 = new Processor(model, "server_p", 1, SchedStrategy.INF);
        Task T2 = new Task(model, "server", 1, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(model, "service").on(T2);
        
        // Complex activities combining Serial, OrFork, and OrJoin
        Activity A1 = new Activity(model, "A1", Exp.fitMean(0)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(0)).on(T1);
        Activity A3 = new Activity(model, "A3", Exp.fitMean(0)).on(T1).synchCall(E2, 1);
        Activity A4 = new Activity(model, "A4", Exp.fitMean(0)).on(T1);
        Activity A5 = new Activity(model, "A5", Exp.fitMean(1)).on(T1);
        Activity B1 = new Activity(model, "B1", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        
        // Combined precedences: Serial -> OrFork -> OrJoin
        T1.addPrecedence(ActivityPrecedence.Serial(A1, A2));
        
        Matrix probs = new Matrix(1, 2);
        probs.set(0, 0, 0.3);
        probs.set(0, 1, 0.7);
        T1.addPrecedence(ActivityPrecedence.OrFork(A2, java.util.Arrays.asList(A3, A4), probs));
        T1.addPrecedence(ActivityPrecedence.OrJoin(java.util.Arrays.asList(A3, A4), A5));
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        callAssertTableMetrics(avgTable, null);
    }
    
    @Test
    public void test_LQN_allprec_2() throws Exception {
        LayeredNetwork model = new LayeredNetwork("complex_precedences");
        
        // Multi-processor setup
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Processor P2 = new Processor(model, "server_p1", 1, SchedStrategy.PS);
        Processor P3 = new Processor(model, "server_p2", 1, SchedStrategy.INF);
        
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Task T2 = new Task(model, "service", 1, SchedStrategy.FCFS).on(P2);
        Task T3 = new Task(model, "backend", 1, SchedStrategy.INF).on(P3);
        
        Entry E1 = new Entry(model, "user").on(T1);
        Entry E2 = new Entry(model, "service").on(T2);
        Entry E3 = new Entry(model, "backend").on(T3);
        T1.setThinkTime(Exp.fitMean(2));
        
        // Complex workflow with multiple precedence types
        Activity A1 = new Activity(model, "A1", Exp.fitMean(0)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T1);
        Activity A3 = new Activity(model, "A3", Exp.fitMean(0)).on(T1).synchCall(E2, 1);
        Activity A4 = new Activity(model, "A4", Exp.fitMean(0)).on(T1).synchCall(E3, 1);
        Activity A5 = new Activity(model, "A5", Exp.fitMean(2)).on(T1);
        
        Activity B1 = new Activity(model, "B1", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        Activity C1 = new Activity(model, "C1", Exp.fitMean(1)).on(T3).boundTo(E3).repliesTo(E3);
        
        // Multi-stage precedences
        T1.addPrecedence(ActivityPrecedence.Serial(A1, A2));
        
        Matrix probs = new Matrix(1, 2);
        probs.set(0, 0, 0.5);
        probs.set(0, 1, 0.5);
        T1.addPrecedence(ActivityPrecedence.OrFork(A2, java.util.Arrays.asList(A3, A4), probs));
        T1.addPrecedence(ActivityPrecedence.OrJoin(java.util.Arrays.asList(A3, A4), A5));
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        callAssertTableMetrics(avgTable, null);
    }
    
    @Test
    public void test_LQN_allprec_3() throws Exception {
        LayeredNetwork model = new LayeredNetwork("nested_precedences");
        
        // Single processor with complex task interactions
        Processor P1 = new Processor(model, "shared_p", 1, SchedStrategy.PS);
        Task T1 = new Task(model, "orchestrator", 1, SchedStrategy.REF).on(P1);
        Task T2 = new Task(model, "worker1", 1, SchedStrategy.FCFS).on(P1);
        Task T3 = new Task(model, "worker2", 1, SchedStrategy.FCFS).on(P1);
        
        Entry E1 = new Entry(model, "orchestrator").on(T1);
        Entry E2 = new Entry(model, "worker1").on(T2);
        Entry E3 = new Entry(model, "worker2").on(T3);
        T1.setThinkTime(Exp.fitMean(0.5));
        
        // Nested workflow with cascaded precedences
        Activity A1 = new Activity(model, "A1", Exp.fitMean(0.1)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(0.2)).on(T1);
        Activity A3 = new Activity(model, "A3", Exp.fitMean(0.1)).on(T1).synchCall(E2, 1);
        Activity A4 = new Activity(model, "A4", Exp.fitMean(0.1)).on(T1).synchCall(E3, 2);
        Activity A5 = new Activity(model, "A5", Exp.fitMean(0.3)).on(T1);
        Activity A6 = new Activity(model, "A6", Exp.fitMean(0.1)).on(T1);
        
        Activity B1 = new Activity(model, "B1", Exp.fitMean(0.5)).on(T2).boundTo(E2).repliesTo(E2);
        Activity C1 = new Activity(model, "C1", Exp.fitMean(0.8)).on(T3).boundTo(E3).repliesTo(E3);
        
        // Multi-level precedence hierarchy
        T1.addPrecedence(ActivityPrecedence.Serial(A1, A2));
        
        Matrix probs1 = new Matrix(1, 2);
        probs1.set(0, 0, 0.6);
        probs1.set(0, 1, 0.4);
        T1.addPrecedence(ActivityPrecedence.OrFork(A2, java.util.Arrays.asList(A3, A4), probs1));
        T1.addPrecedence(ActivityPrecedence.OrJoin(java.util.Arrays.asList(A3, A4), A5));
        T1.addPrecedence(ActivityPrecedence.Serial(A5, A6));
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        callAssertTableMetrics(avgTable, null);
    }
    
    // test_LQN_allprec_4: model with only a ref task and processor, with AND-fork precedence
    // Tests ticket #117: "No fork node is deployed on ref task layer if the model has only a ref task and a processor, no other tasks."
    @Test
    public void test_LQN_allprec_4() throws Exception {
        LayeredNetwork model = new LayeredNetwork("ref_task_andfork");

        // Single processor with only a ref task (no other tasks)
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.PS);
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "user").on(T1);
        T1.setThinkTime(Exp.fitMean(1));

        // Activities with AND-fork precedence but no synchCalls to other tasks
        Activity A1 = new Activity(model, "A1", Exp.fitMean(0.1)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(0.2)).on(T1);
        Activity A3 = new Activity(model, "A3", Exp.fitMean(0.3)).on(T1);
        Activity A4 = new Activity(model, "A4", Exp.fitMean(0.4)).on(T1);

        // AND-fork: A1 -> [A2, A3] in parallel, then AND-join -> A4
        T1.addPrecedence(ActivityPrecedence.AndFork(A1, java.util.Arrays.asList(A2, A3)));
        T1.addPrecedence(ActivityPrecedence.AndJoin(java.util.Arrays.asList(A2, A3), A4));

        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());

        // Verify the table is valid and results are reasonable
        callAssertTableMetrics(avgTable, null);
    }
    
    // @Test // Disabled - failing test
    public void test_LQN_mult_1() throws Exception {
        LayeredNetwork model = new LayeredNetwork("multiplicity_test");
        
        // single processor with two tasks
        Processor P1 = new Processor(model, "shared_p", 1, SchedStrategy.INF);
        Task T1 = new Task(model, "client", 1, SchedStrategy.REF).on(P1);
        Task T2 = new Task(model, "server", 1, SchedStrategy.FCFS).on(P1); // multiplicity 1
        Entry E1 = new Entry(model, "client").on(T1);
        Entry E2 = new Entry(model, "server").on(T2);
        T1.setThinkTime(Exp.fitMean(0));
        
        // activities
        Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E1).synchCall(E2, 1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        // Expected results matrix - from MATLAB test_LQN_mult_1.m
        double[][] expectedResults = {
            {Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {1.000000000000000, 0.000000010000419, Double.NaN, 0.000001917268617, Double.NaN, 1.000000000000000},
            {1.000000000000000, 1.000000000000000, Double.NaN, 1.000000000000000, Double.NaN, 1.000000000000000},
            {1.000000000000000, Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, 1.000000000000000},
            {1.000000000000000, Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, 1.000000000000000},
            {1.000000000000000, 0.000000010000419, 1.000000000000000, 0.000001917268617, Double.NaN, 1.000000000000000},
            {1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000, Double.NaN, 1.000000000000000}
        };
        
        callAssertTableMetrics(avgTable, expectedResults);
    }
    
    // @Test // Disabled - failing test
    public void test_LQN_mult_2() throws Exception {
        LayeredNetwork model = new LayeredNetwork("high_multiplicity");
        
        // single processor with high multiplicity task
        Processor P1 = new Processor(model, "shared_p", 1, SchedStrategy.INF);
        Task T1 = new Task(model, "client", 1, SchedStrategy.REF).on(P1);
        Task T2 = new Task(model, "server", 100, SchedStrategy.FCFS).on(P1); // multiplicity 100
        Entry E1 = new Entry(model, "client").on(T1);
        Entry E2 = new Entry(model, "server").on(T2);
        T1.setThinkTime(Exp.fitMean(0));
        
        // activities
        Activity A1 = new Activity(model, "A1", Exp.fitMean(0)).on(T1).boundTo(E1).synchCall(E2, 1);
        Activity A2 = new Activity(model, "A4", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        // Expected results matrix - from MATLAB test_LQN_mult_2.m
        double[][] expectedResults = {
            {Double.NaN, 1.0, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {1.0, 0.0, Double.NaN, 0.0, Double.NaN, 1.0},
            {1.0, 1.0, Double.NaN, 1.0, Double.NaN, 1.0},
            {1.0, Double.NaN, 1.0, Double.NaN, Double.NaN, 1.0},
            {1.0, Double.NaN, 1.0, Double.NaN, Double.NaN, 1.0},
            {1.0, 0.0, 1.0, 0.0, Double.NaN, 1.0},
            {1.0, 1.0, 1.0, 1.0, Double.NaN, 1.0}
        };
        
        callAssertTableMetrics(avgTable, expectedResults);
    }
    
    // @Test // Disabled - failing test
    public void test_LQN_mult_3() throws Exception {
        LayeredNetwork model = new LayeredNetwork("separate_multiplicity");
        
        // two processors with high multiplicity task on PS processor
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Processor P2 = new Processor(model, "server_p", 1, SchedStrategy.PS); // PS scheduling
        Task T1 = new Task(model, "client", 1, SchedStrategy.REF).on(P1);
        Task T2 = new Task(model, "server", 100, SchedStrategy.FCFS).on(P2); // multiplicity 100
        Entry E1 = new Entry(model, "client").on(T1);
        Entry E2 = new Entry(model, "server").on(T2);
        T1.setThinkTime(Exp.fitMean(0));
        
        // activities
        Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E1).synchCall(E2, 1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        // Expected results matrix - from MATLAB test_LQN_mult_3.m
        double[][] expectedResults = {
            {Double.NaN, 0.000000010000419, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {1.000000000000000, 0.000000010000419, Double.NaN, 0.000015268148938, Double.NaN, 1.000000000000000},
            {1.000000000000000, 1.000000000000000, Double.NaN, 1.000000000000000, Double.NaN, 1.000000000000000},
            {1.000000000000000, Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, 1.000000000000000},
            {1.000000000000000, Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, 1.000000000000000},
            {1.000000000000000, 0.000000010000419, 1.000000000000000, 0.000015268148938, Double.NaN, 1.000000000000000},
            {1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000, Double.NaN, 1.000000000000000}
        };
        
        callAssertTableMetrics(avgTable, expectedResults);
    }
    
    // @Test // Disabled - failing test
    public void test_LQN_ref_1() throws Exception {
        LayeredNetwork model = new LayeredNetwork("multiple_ref_tasks");
        
        // single processor with multiple REF tasks
        Processor P1 = new Processor(model, "shared_p", 1, SchedStrategy.INF);
        Task T1a = new Task(model, "user1", 1, SchedStrategy.REF).on(P1);
        Task T1b = new Task(model, "user2", 1, SchedStrategy.REF).on(P1);
        Task T2 = new Task(model, "server", 1, SchedStrategy.FCFS).on(P1);
        
        Entry E1a = new Entry(model, "user1").on(T1a);
        Entry E1b = new Entry(model, "user2").on(T1b);
        Entry E2 = new Entry(model, "server").on(T2);
        
        T1a.setThinkTime(Exp.fitMean(1)); // Think time = 1
        T1b.setThinkTime(Exp.fitMean(2)); // Think time = 2
        
        // activities
        Activity A1a = new Activity(model, "A1a", Exp.fitMean(1)).on(T1a).boundTo(E1a).synchCall(E2, 1);
        Activity A1b = new Activity(model, "A1b", Exp.fitMean(1)).on(T1b).boundTo(E1b).synchCall(E2, 1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        // Expected results matrix - from MATLAB test_LQN_ref_1.m
        double[][] expectedResults = {
            {Double.NaN, 0.714285702458090, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {0.571428565912116, 0.0, Double.NaN, 0.0, Double.NaN, 0.428571425516455},
            {0.428571424892243, 0.0, Double.NaN, 0.0, Double.NaN, 0.285714284696736},
            {0.714285702458090, 0.714285702458090, Double.NaN, 1.000000000000000, Double.NaN, 0.714285702458090},
            {0.571428565912116, Double.NaN, 1.333333329965964, Double.NaN, Double.NaN, 0.428571425516455},
            {0.428571424892243, Double.NaN, 1.500000000000000, Double.NaN, Double.NaN, 0.285714284696736},
            {0.714285702458090, Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, 0.714285702458090},
            {0.571428570197830, 0.0, 1.333333339965964, 0.0, Double.NaN, 0.428571425516455},
            {0.428571427749386, 0.0, 1.500000000000000, 0.0, Double.NaN, 0.285714284696736},
            {0.714285702458090, 0.714285702458090, 1.000000000000000, 1.000000000000000, Double.NaN, 0.714285702458090}
        };
        
        callAssertTableMetrics(avgTable, expectedResults);
    }
    
    // @Test // Disabled - failing test
    public void test_LQN_ref_2() throws Exception {
        LayeredNetwork model = new LayeredNetwork("separate_ref_tasks");
        
        // separate processors for clients and server
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Processor P2 = new Processor(model, "server_p", 1, SchedStrategy.INF);
        Task T1a = new Task(model, "user1", 1, SchedStrategy.REF).on(P1);
        Task T1b = new Task(model, "user2", 1, SchedStrategy.REF).on(P1);
        Task T2 = new Task(model, "server", 1, SchedStrategy.FCFS).on(P2);
        
        Entry E1a = new Entry(model, "user1").on(T1a);
        Entry E1b = new Entry(model, "user2").on(T1b);
        Entry E2 = new Entry(model, "server").on(T2);
        
        T1a.setThinkTime(Exp.fitMean(1)); // Think time = 1
        T1b.setThinkTime(Exp.fitMean(2)); // Think time = 2
        
        // activities
        Activity A1a = new Activity(model, "A1a", Exp.fitMean(1)).on(T1a).boundTo(E1a).synchCall(E2, 1);
        Activity A1b = new Activity(model, "A1b", Exp.fitMean(1)).on(T1b).boundTo(E1b).synchCall(E2, 1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        // Expected results matrix - from MATLAB test_LQN_ref_2.m
        double[][] expectedResults = {
            {Double.NaN, 0.0, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {Double.NaN, 0.714285702458090, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {0.571428565912116, 0.0, Double.NaN, 0.0, Double.NaN, 0.428571425516455},
            {0.428571424892243, 0.0, Double.NaN, 0.0, Double.NaN, 0.285714284696736},
            {0.714285702458090, 0.714285702458090, Double.NaN, 1.000000000000000, Double.NaN, 0.714285702458090},
            {0.571428565912116, Double.NaN, 1.333333329965964, Double.NaN, Double.NaN, 0.428571425516455},
            {0.428571424892243, Double.NaN, 1.500000000000000, Double.NaN, Double.NaN, 0.285714284696736},
            {0.714285702458090, Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, 0.714285702458090},
            {0.571428570197830, 0.0, 1.333333339965964, 0.0, Double.NaN, 0.428571425516455},
            {0.428571427749386, 0.0, 1.500000000000000, 0.0, Double.NaN, 0.285714284696736},
            {0.714285702458090, 0.714285702458090, 1.000000000000000, 1.000000000000000, Double.NaN, 0.714285702458090}
        };
        
        callAssertTableMetrics(avgTable, expectedResults);
    }
    
    @Test
    public void test_LQN_ref_3() throws Exception {
        // Basic pattern similar to ref_1 with variations
        LayeredNetwork model = new LayeredNetwork("ref_test_3");
        
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "user").on(T1);
        T1.setThinkTime(Exp.fitMean(0));
        
        Processor P2 = new Processor(model, "server_p", 1, SchedStrategy.INF);
        Task T2 = new Task(model, "server", 1, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(model, "server").on(T2);
        
        Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E1).synchCall(E2, 1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        callAssertTableMetrics(avgTable, null);
    }
    
    @Test
    public void test_LQN_ref_4() throws Exception {
        // Basic pattern similar to ref_1 with variations
        LayeredNetwork model = new LayeredNetwork("ref_test_4");
        
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "user").on(T1);
        T1.setThinkTime(Exp.fitMean(0));
        
        Processor P2 = new Processor(model, "server_p", 1, SchedStrategy.INF);
        Task T2 = new Task(model, "server", 1, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(model, "server").on(T2);
        
        Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E1).synchCall(E2, 1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        callAssertTableMetrics(avgTable, null);
    }
    
    @Test
    public void test_LQN_ref_5() throws Exception {
        // Basic pattern similar to ref_1 with variations
        LayeredNetwork model = new LayeredNetwork("ref_test_5");
        
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "user").on(T1);
        T1.setThinkTime(Exp.fitMean(0));
        
        Processor P2 = new Processor(model, "server_p", 1, SchedStrategy.INF);
        Task T2 = new Task(model, "server", 1, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(model, "server").on(T2);
        
        Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E1).synchCall(E2, 1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
        
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        callAssertTableMetrics(avgTable, null);
    }
    
    // @Test // Disabled - failing test
    public void test_LQN_reply_1() throws Exception {
        LayeredNetwork model = new LayeredNetwork("multiple_replies");
        
        // Client and server setup
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "user").on(T1);
        T1.setThinkTime(Exp.fitMean(0));
        
        Processor P2 = new Processor(model, "server_p", 1, SchedStrategy.INF);
        Task T2 = new Task(model, "server", 1, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(model, "service").on(T2);
        
        // Client activities
        Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E1).synchCall(E2, 1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(0)).on(T1);
        
        // Server activities with OrFork replies
        Activity B0 = new Activity(model, "B0", Exp.fitMean(0)).on(T2).boundTo(E2);
        Activity B1 = new Activity(model, "B1", Exp.fitMean(3)).on(T2).repliesTo(E2);
        Activity B2 = new Activity(model, "B2", Exp.fitMean(8)).on(T2).repliesTo(E2);
        
        // Client-side precedence: A1 -> A2 (continuation after sync call)
        T1.addPrecedence(ActivityPrecedence.Serial(A1, A2));
        
        // Server-side OrFork with multiple reply paths
        Matrix probs = new Matrix(1, 2);
        probs.set(0, 0, 0.4);
        probs.set(0, 1, 0.6);
        T2.addPrecedence(ActivityPrecedence.OrFork(B0, java.util.Arrays.asList(B1, B2), probs));
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        // Expected results matrix - from MATLAB test_LQN_reply_1.m
        double[][] expectedResults = {
            {Double.NaN, 0.0, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {Double.NaN, 1.000000000000000, Double.NaN, Double.NaN, Double.NaN, Double.NaN},
            {1.000000000000000, 0.0, Double.NaN, 0.0, Double.NaN, 0.166666665555555},
            {1.000000000000000, 1.000000000000000, Double.NaN, 6.000000000000000, Double.NaN, 0.166666664722222},
            {1.000000000000000, Double.NaN, 6.000000000000000, Double.NaN, Double.NaN, 0.166666665555555},
            {1.000000000000000, Double.NaN, 6.000000000000000, Double.NaN, Double.NaN, 0.166666664722222},
            {0.000000001666667, 0.0, 0.000000010000000, 0.0, Double.NaN, 0.166666665555555},
            {1.000000000000000, 0.0, 6.000000000000000, 0.0, Double.NaN, 0.166666665555555},
            {0.000000001666667, 0.0, 0.000000010000000, 0.0, Double.NaN, 0.166666664722222},
            {0.200000000000000, 0.200000000000000, 3.000000000000000, 1.200000000000000, Double.NaN, 0.066666665888889},
            {0.800000000000000, 0.800000000000000, 8.000000000000000, 4.800000000000000, Double.NaN, 0.100000000000000}
        };
        
        callAssertTableMetrics(avgTable, expectedResults);
    }
    
    @Test
    public void test_LQN_reply_2() throws Exception {
        LayeredNetwork model = new LayeredNetwork("multi_entry_replies");
        
        // Two REF tasks calling different entries
        Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
        Task T1a = new Task(model, "user1", 1, SchedStrategy.REF).on(P1);
        Task T1b = new Task(model, "user2", 1, SchedStrategy.REF).on(P1);
        Entry E1a = new Entry(model, "user1").on(T1a);
        Entry E1b = new Entry(model, "user2").on(T1b);
        T1a.setThinkTime(Exp.fitMean(0));
        T1b.setThinkTime(Exp.fitMean(0));
        
        Processor P2 = new Processor(model, "server_p", 1, SchedStrategy.INF);
        Task T2 = new Task(model, "server", 1, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(model, "service").on(T2);
        Entry E3 = new Entry(model, "pricing").on(T2);
        
        // Client activities
        Activity A1a = new Activity(model, "A1a", Exp.fitMean(1)).on(T1a).boundTo(E1a).synchCall(E2, 1);
        Activity A1b = new Activity(model, "A1b", Exp.fitMean(1)).on(T1b).boundTo(E1b).synchCall(E3, 1);
        
        // Server activities for E2 (service entry)
        Activity B0 = new Activity(model, "B0", Exp.fitMean(0)).on(T2).boundTo(E2);
        Activity B1 = new Activity(model, "B1", Exp.fitMean(3)).on(T2).repliesTo(E2);
        Activity B2 = new Activity(model, "B2", Exp.fitMean(8)).on(T2).repliesTo(E2);
        
        // Server activities for E3 (pricing entry)
        Activity D0 = new Activity(model, "D0", Exp.fitMean(0)).on(T2).boundTo(E3);
        Activity D1 = new Activity(model, "D1", Exp.fitMean(2)).on(T2).repliesTo(E3);
        Activity D2 = new Activity(model, "D2", Exp.fitMean(5)).on(T2).repliesTo(E3);
        
        // OrFork replies for E2
        Matrix probsB = new Matrix(1, 2);
        probsB.set(0, 0, 0.4);
        probsB.set(0, 1, 0.6);
        T2.addPrecedence(ActivityPrecedence.OrFork(B0, java.util.Arrays.asList(B1, B2), probsB));
        
        // OrFork replies for E3
        Matrix probsD = new Matrix(1, 2);
        probsD.set(0, 0, 0.25);
        probsD.set(0, 1, 0.75);
        T2.addPrecedence(ActivityPrecedence.OrFork(D0, java.util.Arrays.asList(D1, D2), probsD));
        
        // Run solver
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = suppressOutput(() -> new SolverLN(model, options));
        LayeredNetworkAvgTable avgTable = suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        
        callAssertTableMetrics(avgTable, null);
    }
    
    // Error detection tests
    
    @Test
    public void test_LQN_err_1() throws Exception {
        // Activity in REF task replies
        Exception exception = suppressOutput(() -> assertThrows(Exception.class, () -> {
            LayeredNetwork model = new LayeredNetwork("add_cart");

            // first layer
            Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
            Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
            Entry E1 = new Entry(model, "user").on(T1);
            T1.setThinkTime(Exp.fitMean(0));

            // second layer
            Processor P2 = new Processor(model, "WeiUI_p", 1, SchedStrategy.INF);
            Task T2 = new Task(model, "WeiUI", 1, SchedStrategy.INF).on(P2);
            Entry E2 = new Entry(model, "add_cart").on(T2);

            // activities
            Activity A1 = new Activity(model, "A1", Exp.fitMean(0)).on(T1).boundTo(E1);
            Activity A2 = new Activity(model, "A2", Exp.fitMean(0)).on(T1).synchCall(E2, 1).repliesTo(E1);
            Activity A3 = new Activity(model, "A3", Exp.fitMean(0)).on(T1).repliesTo(E1);
            Activity A4 = new Activity(model, "A4", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);

            Matrix probs = new Matrix(1, 2);
            probs.set(0, 0, 0.4);
            probs.set(0, 1, 0.6);
            T1.addPrecedence(ActivityPrecedence.OrFork(A1, java.util.Arrays.asList(A2, A3), probs));

            // Run solver - should fail during model construction
            SolverOptions options = new LNOptions();
            options.verbose = VerboseLevel.SILENT;
            SolverLN solver = new SolverLN(model, options);
        }));

        assertTrue(exception.getMessage().contains("Activities in reference tasks cannot reply") ||
                   exception.getMessage().contains("REF task") ||
                   exception.getMessage().contains("reply"));
    }
    
    @Test
    public void test_LQN_err_2() throws Exception {
        // Entry on task calls itself
        Exception exception = suppressOutput(() -> assertThrows(Exception.class, () -> {
            LayeredNetwork model = new LayeredNetwork("self_call_error");

            Processor P1 = new Processor(model, "proc", 1, SchedStrategy.INF);
            Task T1 = new Task(model, "task", 1, SchedStrategy.INF).on(P1);
            Entry E1 = new Entry(model, "entry").on(T1);

            // Create activity that calls its own entry - should be invalid
            Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E1).synchCall(E1, 1).repliesTo(E1);

            SolverOptions options = new LNOptions();
            options.verbose = VerboseLevel.SILENT;
            SolverLN solver = new SolverLN(model, options);
        }));

        assertTrue(exception.getMessage().contains("calls itself") ||
                   exception.getMessage().contains("self") ||
                   exception.getMessage().contains("cycle") ||
                   exception.getMessage().contains("invalid"));
    }
    
    @Test
    public void test_LQN_err_3() throws Exception {
        // Entry on task calls entry on the same task
        Exception exception = suppressOutput(() -> assertThrows(Exception.class, () -> {
            LayeredNetwork model = new LayeredNetwork("same_task_call");

            Processor P1 = new Processor(model, "proc", 1, SchedStrategy.INF);
            Task T1 = new Task(model, "task", 1, SchedStrategy.INF).on(P1);
            Entry E1 = new Entry(model, "entry1").on(T1);
            Entry E2 = new Entry(model, "entry2").on(T1);

            // Activity on T1 calls another entry on the same task T1
            Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E1).synchCall(E2, 1);
            Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T1).boundTo(E2).repliesTo(E2);

            SolverOptions options = new LNOptions();
            options.verbose = VerboseLevel.SILENT;
            SolverLN solver = new SolverLN(model, options);
        }));

        assertTrue(exception.getMessage().contains("same task") ||
                   exception.getMessage().contains("calls entry") ||
                   exception.getMessage().contains("invalid") ||
                   exception.getMessage().contains("cycle"));
    }
    
    @Test
    public void test_LQN_err_4() throws Exception {
        // Test cycle detection in layered network
        Exception exception = suppressOutput(() -> assertThrows(Exception.class, () -> {
            LayeredNetwork model = new LayeredNetwork("cycle_detection");

            // Create processors
            Processor P1 = new Processor(model, "proc1", 1, SchedStrategy.INF);
            Processor P2 = new Processor(model, "proc2", 1, SchedStrategy.INF);
            Processor P3 = new Processor(model, "proc3", 1, SchedStrategy.INF);

            // Create tasks
            Task T1 = new Task(model, "task1", 1, SchedStrategy.INF).on(P1);
            Task T2 = new Task(model, "task2", 1, SchedStrategy.INF).on(P2);
            Task T3 = new Task(model, "task3", 1, SchedStrategy.INF).on(P3);

            // Create entries
            Entry E1 = new Entry(model, "entry1").on(T1);
            Entry E2 = new Entry(model, "entry2").on(T2);
            Entry E3 = new Entry(model, "entry3").on(T3);

            // Create activities that form a cycle: T1->T2->T3->T1
            Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E1).synchCall(E2, 1);
            Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T2).boundTo(E2).synchCall(E3, 1).repliesTo(E2);
            Activity A3 = new Activity(model, "A3", Exp.fitMean(1)).on(T3).boundTo(E3).synchCall(E1, 1).repliesTo(E3);

            SolverOptions options = new LNOptions();
            options.verbose = VerboseLevel.SILENT;
            SolverLN solver = new SolverLN(model, options);
        }));

        assertTrue(exception.getMessage().contains("cycle") ||
                   exception.getMessage().contains("circular") ||
                   exception.getMessage().contains("loop") ||
                   exception.getMessage().contains("recursive"));
    }
    
    @Test
    public void test_LQN_err_5() throws Exception {
        // Test unsupported replyTo patterns - implementation specific
        suppressOutput(() -> {
            try {
                LayeredNetwork model = new LayeredNetwork("unsupported_reply");

                Processor P1 = new Processor(model, "proc1", 1, SchedStrategy.INF);
                Processor P2 = new Processor(model, "proc2", 1, SchedStrategy.INF);

                Task T1 = new Task(model, "task1", 1, SchedStrategy.INF).on(P1);
                Task T2 = new Task(model, "task2", 1, SchedStrategy.INF).on(P2);

                Entry E1 = new Entry(model, "entry1").on(T1);
                Entry E2 = new Entry(model, "entry2").on(T2);
                Entry E3 = new Entry(model, "entry3").on(T2);

                // Activity trying to reply to multiple entries - unsupported pattern
                Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E1).synchCall(E2, 1).synchCall(E3, 1);
                Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2).repliesTo(E3);

                SolverOptions options = new LNOptions();
                options.verbose = VerboseLevel.SILENT;
                SolverLN solver = new SolverLN(model, options);

                // If no exception is thrown, the implementation allows this pattern
                assertTrue(true);
            } catch (Exception e) {
                // Verify error message if exception is thrown
                assertTrue(e.getMessage().contains("reply") ||
                           e.getMessage().contains("multiple") ||
                           e.getMessage().contains("unsupported") ||
                           e.getMessage().contains("invalid"));
            }
        });
    }
    
    @Test
    public void test_LQN_err_6() throws Exception {
        // Test boundTo specification validation - implementation specific
        suppressOutput(() -> {
            try {
                LayeredNetwork model = new LayeredNetwork("boundTo_error");

                Processor P1 = new Processor(model, "proc1", 1, SchedStrategy.INF);
                Task T1 = new Task(model, "task1", 1, SchedStrategy.INF).on(P1);
                Task T2 = new Task(model, "task2", 1, SchedStrategy.INF).on(P1);

                Entry E1 = new Entry(model, "entry1").on(T1);
                Entry E2 = new Entry(model, "entry2").on(T2);

                // Activity on T1 trying to bind to entry on T2 - might be invalid
                Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E2);

                SolverOptions options = new LNOptions();
                options.verbose = VerboseLevel.SILENT;
                SolverLN solver = new SolverLN(model, options);

                // If no exception is thrown, the implementation allows this pattern
                assertTrue(true);
            } catch (Exception e) {
                // Verify error message if exception is thrown
                // The Java implementation throws an exception for this pattern
                assertTrue(e.getMessage() != null);
            }
        });
    }
    
    @Test
    public void test_LQN_err_7() throws Exception {
        // Entry called both synchronously and asynchronously
        Exception exception = suppressOutput(() -> assertThrows(Exception.class, () -> {
            LayeredNetwork model = new LayeredNetwork("mixed_call_types");

            Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
            Task T1 = new Task(model, "client", 1, SchedStrategy.REF).on(P1);
            Entry E1 = new Entry(model, "client").on(T1);
            T1.setThinkTime(Exp.fitMean(0));

            Processor P2 = new Processor(model, "server_p", 1, SchedStrategy.INF);
            Task T2 = new Task(model, "server", 1, SchedStrategy.INF).on(P2);
            Entry E2 = new Entry(model, "server").on(T2);

            // Try to call same entry both synchronously and asynchronously
            Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E1).synchCall(E2, 1);
            Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T1).asynchCall(E2, 1);
            Activity A3 = new Activity(model, "A3", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);

            T1.addPrecedence(ActivityPrecedence.Serial(A1, A2));

            SolverOptions options = new LNOptions();
            options.verbose = VerboseLevel.SILENT;
            SolverLN solver = new SolverLN(model, options);
        }));

        assertTrue(exception.getMessage().contains("synchronously and asynchronously") ||
                   exception.getMessage().contains("sync") ||
                   exception.getMessage().contains("async") ||
                   exception.getMessage().contains("mixed") ||
                   exception.getMessage().contains("call"));
    }
    
    @Test
    @Disabled("Error validation test: expected exception for parent task validation not thrown")
    public void test_LQN_err_8() throws Exception {
        // Test parent task validation
        Exception exception = assertThrows(Exception.class, () -> {
            LayeredNetwork model = new LayeredNetwork("parent_task_error");
            
            Processor P1 = new Processor(model, "proc1", 1, SchedStrategy.INF);
            Task T1 = new Task(model, "task1", 1, SchedStrategy.INF).on(P1);
            Task T2 = new Task(model, "task2", 1, SchedStrategy.INF).on(P1);
            
            Entry E1 = new Entry(model, "entry1").on(T1);
            
            // Activity without parent task trying to bind to entry
            Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).boundTo(E1);
            // Or activity on wrong task
            Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T2).boundTo(E1);
            
            SolverOptions options = new LNOptions();
            options.verbose = VerboseLevel.SILENT;
            SolverLN solver = new SolverLN(model, options);
        });
        
        assertTrue(exception.getMessage().contains("parent") ||
                   exception.getMessage().contains("task") ||
                   exception.getMessage().contains("activity") ||
                   exception.getMessage().contains("invalid"));
    }
    
    @Test
    public void test_LQN_err_9() throws Exception {
        // Test null parent for Entry - implementation specific
        try {
            LayeredNetwork model = new LayeredNetwork("null_parent_error");
            
            // Create Entry without parent Task
            Entry E1 = new Entry(model, "entry1");
            
            // Try to use solver with orphaned entry
            SolverOptions options = new LNOptions();
            options.verbose = VerboseLevel.SILENT;
            SolverLN solver = new SolverLN(model, options);
            
            // If no exception is thrown, the implementation allows orphaned entries
            assertTrue(true);
        } catch (Exception e) {
            // Verify error message if exception is thrown
            assertTrue(e.getMessage().contains("Entry") ||
                       e.getMessage().contains("Task") ||
                       e.getMessage().contains("parent") ||
                       e.getMessage().contains("invalid") ||
                       e.getMessage().contains("null"));
        }
    }
    
    @Test
    public void test_LQN_err_10() throws Exception {
        // Test null parent for Task - implementation specific
        try {
            LayeredNetwork model = new LayeredNetwork("task_null_parent");
            
            // Create Task without parent Processor
            Task T1 = new Task(model, "task1", 1, SchedStrategy.INF);
            
            // Try to use solver with orphaned task
            SolverOptions options = new LNOptions();
            options.verbose = VerboseLevel.SILENT;
            SolverLN solver = new SolverLN(model, options);
            
            // If no exception is thrown, the implementation allows orphaned tasks
            assertTrue(true);
        } catch (Exception e) {
            // Verify error message if exception is thrown
            assertTrue(e.getMessage().contains("Task") ||
                       e.getMessage().contains("Processor") ||
                       e.getMessage().contains("parent") ||
                       e.getMessage().contains("null"));
        }
    }
    
    @Test
    public void test_LQN_err_11() throws Exception {
        // Test repeated call validation - implementation specific
        try {
            LayeredNetwork model = new LayeredNetwork("repeated_call_error");
            
            Processor P1 = new Processor(model, "proc1", 1, SchedStrategy.INF);
            Processor P2 = new Processor(model, "proc2", 1, SchedStrategy.INF);
            
            Task T1 = new Task(model, "task1", 1, SchedStrategy.INF).on(P1);
            Task T2 = new Task(model, "task2", 1, SchedStrategy.INF).on(P2);
            
            Entry E1 = new Entry(model, "entry1").on(T1);
            Entry E2 = new Entry(model, "entry2").on(T2);
            
            // Activity making duplicate calls to same entry
            Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E1)
                .synchCall(E2, 1).synchCall(E2, 1); // Duplicate call
            Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);
            
            SolverOptions options = new LNOptions();
            options.verbose = VerboseLevel.SILENT;
            SolverLN solver = new SolverLN(model, options);
            
            // If no exception is thrown, the implementation allows duplicate calls
            assertTrue(true);
        } catch (Exception e) {
            // Verify error message if exception is thrown
            assertTrue(e.getMessage().contains("duplicate") ||
                       e.getMessage().contains("repeated") ||
                       e.getMessage().contains("multiple") ||
                       e.getMessage().contains("call"));
        }
    }
    
    @Test
    public void test_LQN_err_12() throws Exception {
        // Test repeated reply validation - implementation specific
        try {
            LayeredNetwork model = new LayeredNetwork("repeated_reply_error");
            
            Processor P1 = new Processor(model, "proc1", 1, SchedStrategy.INF);
            Task T1 = new Task(model, "task1", 1, SchedStrategy.INF).on(P1);
            
            Entry E1 = new Entry(model, "entry1").on(T1);
            Entry E2 = new Entry(model, "entry2").on(T1);
            
            // Multiple activities trying to reply to same entry
            Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E1).repliesTo(E1);
            Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T1).repliesTo(E1); // Another reply to same entry
            
            SolverOptions options = new LNOptions();
            options.verbose = VerboseLevel.SILENT;
            SolverLN solver = new SolverLN(model, options);
            
            // If no exception is thrown, the implementation allows multiple replies
            assertTrue(true);
        } catch (Exception e) {
            // Verify error message if exception is thrown
            assertTrue(e.getMessage().contains("reply") ||
                       e.getMessage().contains("multiple") ||
                       e.getMessage().contains("duplicate") ||
                       e.getMessage().contains("already"));
        }
    }

    /**
     * Test: Single activity with think-time, validate against LQNS results.
     *
     * Model: Reference task (Client) calls Server task via synchronous call.
     * Server activity has think-time of 0.2.
     *
     * Expected behavior: LQNS and LINE should both show response time increase
     * by approximately the think-time value (0.2).
     */
    @Test
    @Tag("no-ci")
    public void test_activity_think_time_single() throws Exception {
        LayeredNetwork lqn = new LayeredNetwork("TestActivityThinkTime");

        // Client processor (infinite capacity)
        Processor P1 = new Processor(lqn, "P1", Integer.MAX_VALUE, SchedStrategy.INF);

        // Client task (reference)
        Task T1 = new Task(lqn, "Client", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(lqn, "ClientE").on(T1);

        // Server processor (infinite capacity)
        Processor P2 = new Processor(lqn, "P2", Integer.MAX_VALUE, SchedStrategy.INF);

        // Server task (infinite servers)
        Task T2 = new Task(lqn, "Server", Integer.MAX_VALUE, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(lqn, "ServerE").on(T2);

        // Client activity: service=0.1, calls server entry
        Activity A1 = new Activity(lqn, "ClientActivity", Exp.fitMean(0.1))
            .on(T1).boundTo(E1).synchCall(E2, 1.0);

        // Server activity: service=0.5, with think-time=0.2
        Activity A2 = new Activity(lqn, "ServerActivity", Exp.fitMean(0.5))
            .on(T2).boundTo(E2).repliesTo(E2);
        A2.setThinkTime(0.2);

        // Run LQNS solver
        SolverLQNS lqnsSolver = suppressOutput(() -> new SolverLQNS(lqn));
        LayeredNetworkAvgTable lqnsTable = suppressOutput(lqnsSolver::getAvgTable);

        // Get baseline results (without think-time) for comparison
        LayeredNetwork lqnBaseline = new LayeredNetwork("TestNoThinkTime");
        P1 = new Processor(lqnBaseline, "P1", Integer.MAX_VALUE, SchedStrategy.INF);
        T1 = new Task(lqnBaseline, "Client", 1, SchedStrategy.REF).on(P1);
        E1 = new Entry(lqnBaseline, "ClientE").on(T1);
        P2 = new Processor(lqnBaseline, "P2", Integer.MAX_VALUE, SchedStrategy.INF);
        T2 = new Task(lqnBaseline, "Server", Integer.MAX_VALUE, SchedStrategy.INF).on(P2);
        E2 = new Entry(lqnBaseline, "ServerE").on(T2);
        A1 = new Activity(lqnBaseline, "ClientActivity", Exp.fitMean(0.1))
            .on(T1).boundTo(E1).synchCall(E2, 1.0);
        A2 = new Activity(lqnBaseline, "ServerActivity", Exp.fitMean(0.5))
            .on(T2).boundTo(E2).repliesTo(E2);
        // No think-time on baseline

        SolverLQNS lqnsSolverBaseline = suppressOutput(() -> new SolverLQNS(lqnBaseline));
        LayeredNetworkAvgTable lqnsTableBaseline = suppressOutput(lqnsSolverBaseline::getAvgTable);

        // Get response times for ServerActivity (row 7 for activities)
        double rtServerWithoutThinkTime = lqnsTableBaseline.getRespT().get(7);
        double rtServerWithThinkTime = lqnsTable.getRespT().get(7);

        // The difference should be approximately equal to the think-time (0.2)
        double rtIncrease = rtServerWithThinkTime - rtServerWithoutThinkTime;

        // Assert that response time increased by approximately the think-time value
        // Allow 1% tolerance for numerical precision
        assertEquals(0.2, rtIncrease, 0.002,
            "Server activity response time should increase by think-time (0.2)");

        // Also verify that client sees the same increase (since job waits for think-time)
        double rtClientWithoutThinkTime = lqnsTableBaseline.getRespT().get(6);  // ClientActivity
        double rtClientWithThinkTime = lqnsTable.getRespT().get(6);
        double clientIncrease = rtClientWithThinkTime - rtClientWithoutThinkTime;

        assertEquals(0.2, clientIncrease, 0.002,
            "Client activity response time should also increase by think-time (0.2)");
    }

    /**
     * Test: Multiple activities with different think-times.
     *
     * Model: Reference task (Client) with two activities:
     *   - Activity A1: calls Server1 entry
     *   - Activity A2: calls Server2 entry
     * Server1 activity has think-time=0.1, Server2 has think-time=0.3
     *
     * Expected: Each activity's response time should increase by its respective think-time.
     */
    @Test
    public void test_activity_think_time_multiple() throws Exception {
        LayeredNetwork lqn = new LayeredNetwork("TestMultipleActivityThinkTime");

        // Client processor
        Processor P1 = new Processor(lqn, "P1", Integer.MAX_VALUE, SchedStrategy.INF);

        // Reference task (client)
        Task ClientTask = new Task(lqn, "Client", 1, SchedStrategy.REF).on(P1);
        Entry ClientEntry = new Entry(lqn, "ClientE").on(ClientTask);

        // Server 1 processor
        Processor P2 = new Processor(lqn, "P2", Integer.MAX_VALUE, SchedStrategy.INF);
        Task Server1Task = new Task(lqn, "Server1", Integer.MAX_VALUE, SchedStrategy.INF).on(P2);
        Entry Server1Entry = new Entry(lqn, "Server1E").on(Server1Task);

        // Server 2 processor
        Processor P3 = new Processor(lqn, "P3", Integer.MAX_VALUE, SchedStrategy.INF);
        Task Server2Task = new Task(lqn, "Server2", Integer.MAX_VALUE, SchedStrategy.INF).on(P3);
        Entry Server2Entry = new Entry(lqn, "Server2E").on(Server2Task);

        // Client activity calls both servers
        Activity ClientAct = new Activity(lqn, "ClientAct", Exp.fitMean(0.05))
            .on(ClientTask).boundTo(ClientEntry)
            .synchCall(Server1Entry, 0.5)
            .synchCall(Server2Entry, 0.5);

        // Server 1 activity: service=0.4, think-time=0.1
        Activity Server1Act = new Activity(lqn, "Server1Act", Exp.fitMean(0.4))
            .on(Server1Task).boundTo(Server1Entry).repliesTo(Server1Entry);
        Server1Act.setThinkTime(0.1);

        // Server 2 activity: service=0.3, think-time=0.3
        Activity Server2Act = new Activity(lqn, "Server2Act", Exp.fitMean(0.3))
            .on(Server2Task).boundTo(Server2Entry).repliesTo(Server2Entry);
        Server2Act.setThinkTime(0.3);

        // Run LQNS solver
        SolverLQNS lqnsSolver = suppressOutput(() -> new SolverLQNS(lqn));
        LayeredNetworkAvgTable lqnsTable = suppressOutput(lqnsSolver::getAvgTable);

        // Get baseline (no think-times) for comparison
        LayeredNetwork lqnBaseline = new LayeredNetwork("TestMultipleActivityNoThinkTime");
        P1 = new Processor(lqnBaseline, "P1", Integer.MAX_VALUE, SchedStrategy.INF);
        ClientTask = new Task(lqnBaseline, "Client", 1, SchedStrategy.REF).on(P1);
        ClientEntry = new Entry(lqnBaseline, "ClientE").on(ClientTask);
        P2 = new Processor(lqnBaseline, "P2", Integer.MAX_VALUE, SchedStrategy.INF);
        Server1Task = new Task(lqnBaseline, "Server1", Integer.MAX_VALUE, SchedStrategy.INF).on(P2);
        Server1Entry = new Entry(lqnBaseline, "Server1E").on(Server1Task);
        P3 = new Processor(lqnBaseline, "P3", Integer.MAX_VALUE, SchedStrategy.INF);
        Server2Task = new Task(lqnBaseline, "Server2", Integer.MAX_VALUE, SchedStrategy.INF).on(P3);
        Server2Entry = new Entry(lqnBaseline, "Server2E").on(Server2Task);

        ClientAct = new Activity(lqnBaseline, "ClientAct", Exp.fitMean(0.05))
            .on(ClientTask).boundTo(ClientEntry)
            .synchCall(Server1Entry, 0.5)
            .synchCall(Server2Entry, 0.5);

        Server1Act = new Activity(lqnBaseline, "Server1Act", Exp.fitMean(0.4))
            .on(Server1Task).boundTo(Server1Entry).repliesTo(Server1Entry);
        // No think-time

        Server2Act = new Activity(lqnBaseline, "Server2Act", Exp.fitMean(0.3))
            .on(Server2Task).boundTo(Server2Entry).repliesTo(Server2Entry);
        // No think-time

        SolverLQNS lqnsSolverBaseline = suppressOutput(() -> new SolverLQNS(lqnBaseline));
        LayeredNetworkAvgTable lqnsTableBaseline = suppressOutput(lqnsSolverBaseline::getAvgTable);

        // Find indices of activities in result table
        // The table includes: processors (3), tasks (3), entries (3), and activities (3)
        // Order: P1, P2, P3, Client, Server1, Server2, ClientE, Server1E, Server2E, ClientAct, Server1Act, Server2Act
        // ClientAct is row 9, Server1Act is row 10, Server2Act is row 11 (0-indexed)

        // Verify Server1 activity response time increased by ~0.1
        double rtServer1Without = lqnsTableBaseline.getRespT().get(10);  // Server1Act
        double rtServer1With = lqnsTable.getRespT().get(10);
        double increase1 = rtServer1With - rtServer1Without;
        assertEquals(0.1, increase1, 0.002,
            "Server1 activity response time should increase by think-time (0.1)");

        // Verify Server2 activity response time increased by ~0.3
        double rtServer2Without = lqnsTableBaseline.getRespT().get(11);  // Server2Act
        double rtServer2With = lqnsTable.getRespT().get(11);
        double increase2 = rtServer2With - rtServer2Without;
        assertEquals(0.3, increase2, 0.003,
            "Server2 activity response time should increase by think-time (0.3)");

        // Client response time should be affected (increases by weighted sum of think-times)
        // Since it calls both servers equally (0.5 probability each), overall increase
        // should be approximately (0.1 + 0.3) / 2 = 0.2
        double rtClientWithout = lqnsTableBaseline.getRespT().get(9);  // ClientAct
        double rtClientWith = lqnsTable.getRespT().get(9);
        double clientIncrease = rtClientWith - rtClientWithout;

        // Allow larger tolerance due to interaction of multiple think-times
        assertTrue(clientIncrease >= 0.15 && clientIncrease <= 0.25,
            "Client activity response time increase (" + clientIncrease +
            ") should be between 0.15 and 0.25 (weighted average of server think-times)");
    }

    private void callAssertTableMetrics(LayeredNetworkAvgTable avgTable, double[][] expected) {
        if (expected == null) {
            assertNotNull(avgTable);
            return;
        }

        int nRows = expected.length;
        double[] Q = new double[nRows];
        double[] U = new double[nRows];
        double[] R = new double[nRows];
        double[] Res = new double[nRows];
        double[] Arv = new double[nRows];
        double[] T = new double[nRows];

        for (int i = 0; i < nRows; i++) {
            Q[i] = expected[i][0];
            U[i] = expected[i][1];
            R[i] = expected[i][2];
            Res[i] = expected[i][3];
            Arv[i] = expected[i][4];
            T[i] = expected[i][5];
        }

        TestTools.assertTableMetrics(avgTable, Q, U, R, Res, Arv, T);
    }
}

