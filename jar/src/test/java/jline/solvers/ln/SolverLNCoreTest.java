package jline.solvers.ln;

import jline.TestTools;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.constant.CallType;
import jline.lang.constant.SchedStrategy;
import jline.VerboseLevel;
import jline.lang.layered.*;
import jline.lang.nodes.ClassSwitch;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVA;

import static jline.TestTools.*;
import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests verifying core LQN decomposition and internal data structures.
 */
class SolverLNCoreTest extends SolverLNTestBase {

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
        assertEquals(lqn.mult.getNumCols(), 5);
        assertEquals(lqn.mult.getNonZeros(), 4);
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
}
