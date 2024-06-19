package jline.solvers.ln;
import jline.lang.*;
import jline.lang.constant.CallType;
import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.Exp;
import jline.lang.distributions.Immediate;
import jline.lang.nodes.ClassSwitch;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.layered.*;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.mva.SolverMVA;
import jline.util.Matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

class SolverLNTest {

    LayeredNetwork buildModel1() throws Exception{
        LayeredNetwork model = new LayeredNetwork("test_LQN_1");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 50, SchedStrategy.REF);
        T1.setThinkTime(new Exp(1.0/100.0));
        T1.on(P1);

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);

        Activity A1 = new Activity(model, "A1", new Exp(10));
        A1.on(T1);
        A1.boundTo(E1);

        Activity A2 = new Activity(model, "A2", new Exp(1.0/1.5));
        A2.on(T1);

        T1.addPrecedence(ActivityPrecedence.Sequence("A1","A2"));
        return model;
    }

    LayeredNetwork buildModel2() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_LQN_2");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 10, SchedStrategy.REF);
        T1.on(P1);
        T1.setThinkTime(new Exp(1.0/100));

        Task T2 = new Task(model, "T2", 1, SchedStrategy.FCFS);
        T2.on(P2);
        T2.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);

        Entry E2 = new Entry(model, "E2");
        E2.on(T2);

        Activity A1 = new Activity(model, "A1", new Exp(1/1.6));
        A1.on(T1);
        A1.boundTo(E1);

        Activity A2 = new Activity(model, "A2", new Immediate());
        A2.on(T1);
        A2.synchCall(E2);

        Activity A3 = new Activity(model, "A3", new Exp(1.0/5));
        A3.on(T2);
        A3.boundTo(E2);

        Activity A4 = new Activity(model, "A4", new Exp(1));
        A4.on(T2);
        A4.repliesTo(E2);

        T1.addPrecedence(ActivityPrecedence.Sequence("A1","A2"));
        T2.addPrecedence(ActivityPrecedence.Sequence("A3","A4"));
        return model;
    }

    // This is the former TestSolverLN2
    LayeredNetwork buildModel3() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_LQN_3");
        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 10, SchedStrategy.REF);
        T1.on(P1);
        T1.setThinkTime(new Exp(1.0/100));

        Task T2 = new Task(model, "T2", 10, SchedStrategy.FCFS);
        T2.on(P1);
        T2.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);

        Entry E2 = new Entry(model, "E2");
        E2.on(T2);

        Activity A1 = new Activity(model, "A1", new Exp(1/70.0));
        A1.on(T1);
        A1.boundTo(E1);

        Activity A2 = new Activity(model, "A2", new Immediate());
        A2.on(T1);
        A2.synchCall(E2);

        Activity A3 = new Activity(model, "A3", new Exp(1.0/5));
        A3.on(T2);
        A3.boundTo(E2);

        Activity A4 = new Activity(model, "A4", new Exp(1));
        A4.on(T2);
        A4.repliesTo(E2);

        T1.addPrecedence(ActivityPrecedence.Sequence("A1","A2"));
        T2.addPrecedence(ActivityPrecedence.Sequence("A3","A4"));
        return model;
    }

    /* Model for calls - seems to run very very slowly with LN, to be checked */
    LayeredNetwork buildModel4() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_LQN_4");

        Processor P1 = new Processor(model, "P1", 2, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 3, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 50, SchedStrategy.REF); T1.on(P1);
        T1.setThinkTime(new Exp((double) 1/ 2));
        Task T2 = new Task(model, "T2", 50, SchedStrategy.FCFS); T2.on(P1);
        T2.setThinkTime(new Exp((double) 1/ 3));
        Task T3 = new Task(model, "T3", 25, SchedStrategy.FCFS); T3.on(P2);
        T3.setThinkTime(new Exp((double) 1 /4));

        Entry E1 = new Entry(model, "E1"); E1.on(T1);
        Entry E2 = new Entry(model, "E2"); E2.on(T2);
        Entry E3 = new Entry(model, "E3"); E3.on(T3);

        Activity A1 = new Activity(model, "AS1", new Exp(10));
        A1.on(T1); A1.boundTo(E1); A1.synchCall(E2,1);
        Activity A2 = new Activity(model, "AS2", new Exp(20));
        A2.on(T2); A2.boundTo(E2); A2.synchCall(E3,5); A2.repliesTo(E2);
        Activity A3 = new Activity(model, "AS3", new Exp(50));
        A3.on(T3); A3.boundTo(E3); A3.repliesTo(E3);

        return model;
    }

    /* Model for Loop Precedences */
    LayeredNetwork buildModel5() throws Exception {
        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", 10, SchedStrategy.INF);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF); T1.on(P1);
        T1.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "Entry"); E1.on(T1);

        Activity A1 = new Activity(model, "A1", new Exp(0.5)); A1.on(T1); A1.boundTo(E1);
        Activity A2 = new Activity(model, "A2", new Exp(0.333333)); A2.on(T1);
        Activity A3 = new Activity(model, "A3", new Exp(0.25)); A3.on(T1);


        // Loop Activity Precedence
        ArrayList<String> precActs = new ArrayList<String>();
        precActs.add("A2");
        precActs.add("A3");
        T1.addPrecedence(ActivityPrecedence.Loop("A1", precActs, new Matrix(2.2)));
        return model;
    }

    LayeredNetwork buildModel6() throws Exception {
        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", 10, SchedStrategy.FCFS);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF); T1.on(P1);
        T1.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "Entry"); E1.on(T1);

        Activity A1 = new Activity(model, "B1", new Exp(0.5)); A1.on(T1); A1.boundTo(E1);
        Activity A2 = new Activity(model, "B2", new Exp(0.333333)); A2.on(T1);
        Activity A3 = new Activity(model, "B3", new Exp(0.25)); A3.on(T1);
        Activity A4 = new Activity(model, "B4", new Exp(0.2)); A4.on(T1);
        Activity A5 = new Activity(model, "B5", new Exp(0.166667)); A5.on(T1);


        // OrFork Activity Precedence
        ArrayList<String> precActs = new ArrayList<String>();
        Matrix probs = new Matrix(1,3);
        precActs.add("B2");
        precActs.add("B3");
        precActs.add("B4");
        probs.set(0,0,0.3);
        probs.set(0,1,0.3);
        probs.set(0,2,0.4);
        T1.addPrecedence(ActivityPrecedence.OrFork("B1", precActs, probs));

        // OrJoin Activity Precedence
        precActs.clear();
        precActs.add("B2");
        precActs.add("B3");
        precActs.add("B4");
        T1.addPrecedence(ActivityPrecedence.OrJoin(precActs, "B5"));

        return model;
    }

    final double errorRate = 0.02;

    boolean compare(double actual, double expected){
        return Math.abs((actual - expected) / expected) <= errorRate;
    }

    boolean compareOutput(double actual, double expected){
        if (Double.isNaN(actual) && Double.isNaN(expected)) {
            return true;
        }
        return Math.abs(actual - expected) <= errorRate;
    }

    @org.junit.jupiter.api.Test
    public void testSingleLayerResult() throws Exception {
        SolverLN solverLN = new SolverLN(buildModel1());
        Network network = solverLN.getEnsemble().get(0);
        SolverMVA layersolver = new SolverMVA(network);
        layersolver.runAnalyzer();

        //test queue length
        assertTrue(compare(layersolver.result.QN.get(0,0),47.38));
        assertTrue(compare(layersolver.result.QN.get(1,2),0.1633));
        assertTrue(compare(layersolver.result.QN.get(1,3),2.45));

        //test through put
        assertTrue(compare(layersolver.result.TN.get(0,0),0.4738));
        assertTrue(compare(layersolver.result.TN.get(1,2),0.4738));
        assertTrue(compare(layersolver.result.TN.get(1,3),0.4738));

        //test utilization
        assertTrue(compare(layersolver.result.UN.get(0,0),47.38));
        assertTrue(compare(layersolver.result.UN.get(1,2),0.04738));
        assertTrue(compare(layersolver.result.UN.get(1,3),0.7108));

        //test response time
        assertTrue(compare(layersolver.result.RN.get(0,0),100));
        assertTrue(compare(layersolver.result.RN.get(1,2),0.3447));
        assertTrue(compare(layersolver.result.RN.get(1,3),5.178));

        //test idxhash
        assertEquals(solverLN.getIdxhash().get(0), Double.NaN);
        assertEquals(solverLN.getIdxhash().get(1), 1);
        assertEquals(solverLN.getIdxhash().get(2), Double.NaN);
    }

    @org.junit.jupiter.api.Test
    public void testMultipleLayerAndMatrix() throws Exception {
        SolverLN solverLN = new SolverLN(buildModel2());

        //Layer 1
        Network network1 = solverLN.getEnsemble().get(0);
        SolverMVA layersolver1 = new SolverMVA(network1);
        layersolver1.runAnalyzer();

        //Layer 2
        Network network2 = solverLN.getEnsemble().get(1);
        SolverMVA layersolver2 = new SolverMVA(network2);
        layersolver2.runAnalyzer();

        //Layer 3
        Network network3 = solverLN.getEnsemble().get(2);
        SolverMVA layersolver3 = new SolverMVA(network3);
        layersolver3.runAnalyzer();

        // test Servt_classes_updmap
        assertEquals(solverLN.getServt_classes_updmap().get(1,1),1);
        assertEquals(solverLN.getServt_classes_updmap().get(1,2),7);
        assertEquals(solverLN.getServt_classes_updmap().get(1,3),2);
        assertEquals(solverLN.getServt_classes_updmap().get(1,4),3);
        assertEquals(solverLN.getServt_classes_updmap().get(2,1),1);
        assertEquals(solverLN.getServt_classes_updmap().get(2,2),8);
        assertEquals(solverLN.getServt_classes_updmap().get(2,3),2);
        assertEquals(solverLN.getServt_classes_updmap().get(2,4),4);
        assertEquals(solverLN.getServt_classes_updmap().get(3,1),2);
        assertEquals(solverLN.getServt_classes_updmap().get(3,2),9);
        assertEquals(solverLN.getServt_classes_updmap().get(3,3),2);
        assertEquals(solverLN.getServt_classes_updmap().get(3,4),3);
        assertEquals(solverLN.getServt_classes_updmap().get(4,1),2);
        assertEquals(solverLN.getServt_classes_updmap().get(4,2),10);
        assertEquals(solverLN.getServt_classes_updmap().get(4,3),2);
        assertEquals(solverLN.getServt_classes_updmap().get(4,4),4);

        // test Call_classes_updmap
        assertEquals(solverLN.getCall_classes_updmap().get(1,1),1);
        assertEquals(solverLN.getCall_classes_updmap().get(1,2),1);
        assertEquals(solverLN.getCall_classes_updmap().get(1,3),1);
        assertEquals(solverLN.getCall_classes_updmap().get(1,4),5);
        assertEquals(solverLN.getCall_classes_updmap().get(2,1),4);
        assertEquals(solverLN.getCall_classes_updmap().get(2,2),1);
        assertEquals(solverLN.getCall_classes_updmap().get(2,3),2);
        assertEquals(solverLN.getCall_classes_updmap().get(2,4),5);

        // test getThinkt_classes_updmap
        assertEquals(solverLN.getThinkt_classes_updmap().get(1,1),2);
        assertEquals(solverLN.getThinkt_classes_updmap().get(1,2),4);
        assertEquals(solverLN.getThinkt_classes_updmap().get(1,3),1);
        assertEquals(solverLN.getThinkt_classes_updmap().get(1,4),1);
        assertEquals(solverLN.getThinkt_classes_updmap().get(2,1),4);
        assertEquals(solverLN.getThinkt_classes_updmap().get(2,2),7);
        assertEquals(solverLN.getThinkt_classes_updmap().get(2,3),1);
        assertEquals(solverLN.getThinkt_classes_updmap().get(2,4),3);
        assertEquals(solverLN.getThinkt_classes_updmap().get(3,1),4);
        assertEquals(solverLN.getThinkt_classes_updmap().get(3,2),8);
        assertEquals(solverLN.getThinkt_classes_updmap().get(3,3),1);
        assertEquals(solverLN.getThinkt_classes_updmap().get(3,4),4);


        //test Layer1
        ////test queue length
        assertTrue(compare(layersolver1.result.QN.get(0,0),9.8189));
        assertTrue(compare(layersolver1.result.QN.get(1,2),0.1811));

        //test through put
        assertTrue(compare(layersolver1.result.TN.get(0,0),0.09819));
        assertTrue(compare(layersolver1.result.TN.get(1,3),0.09819));

        //test utilization
        assertTrue(compare(layersolver1.result.UN.get(0,0),9.819));
        assertTrue(compare(layersolver1.result.UN.get(1,2),0.1571));

        //test response time
        assertTrue(compare(layersolver1.result.RN.get(0,0),100));
        assertTrue(compare(layersolver1.result.RN.get(1,2),1.844));


        //test layer2
        //test queue length
        assertTrue(compare(layersolver2.result.QN.get(1,2),0.8333));
        assertTrue(compare(layersolver2.result.QN.get(1,3),0.1667));

        //test through put
        assertTrue(compare(layersolver2.result.TN.get(1,2),0.1667));
        assertTrue(compare(layersolver2.result.TN.get(1,3),0.1667));

        //test utilization
        assertTrue(compare(layersolver2.result.UN.get(1,2),0.8333));
        assertTrue(compare(layersolver2.result.UN.get(1,3),0.1667));

        //test response time
        assertTrue(compare(layersolver2.result.RN.get(1,2),5));
        assertTrue(compare(layersolver2.result.RN.get(1,3),1));


        //test layer3
        //test queue length
        assertTrue(compare(layersolver3.result.QN.get(0,0),9.8425));
        assertTrue(compare(layersolver3.result.QN.get(0,2),0.1575));

        //test through put
        assertTrue(compare(layersolver3.result.TN.get(0,0),0.09843));
        assertTrue(compare(layersolver3.result.TN.get(0,2),0.09843));

        //test utilization
        assertTrue(compare(layersolver3.result.UN.get(0,0),9.8425));
        assertTrue(compare(layersolver3.result.UN.get(0,2),0.1575));

        //test response time
        assertTrue(compare(layersolver3.result.RN.get(0,0),100));
        assertTrue(compare(layersolver3.result.RN.get(0,2),1.6));

        //test idxhash
        assertEquals(solverLN.getIdxhash().get(0), Double.NaN);
        assertEquals(solverLN.getIdxhash().get(1), 1);
        assertEquals(solverLN.getIdxhash().get(2), 2);
        assertEquals(solverLN.getIdxhash().get(3), Double.NaN);
        assertEquals(solverLN.getIdxhash().get(4), 3);
    }

    @org.junit.jupiter.api.Test
    public void testgetStruct() throws Exception{
        LayeredNetworkStruct lqn = buildModel2().getStruct();
        assertEquals(lqn.nidx,10);
        assertEquals(lqn.nhosts,2);
        assertEquals(lqn.ntasks,2);
        assertEquals(lqn.nentries,2);
        assertEquals(lqn.nacts,4);
        assertEquals(lqn.ncalls,1);
        assertEquals(lqn.hshift,0);
        assertEquals(lqn.tshift,2);
        assertEquals(lqn.eshift,4);
        assertEquals(lqn.ashift,6);
        assertEquals(lqn.cshift,0);

        assertEquals(lqn.tasksof.size(),2);
        assertEquals(lqn.tasksof.get(1).size(),1);
        assertEquals(lqn.tasksof.get(2).size(),1);
        assertEquals(lqn.tasksof.get(1).get(0),3);
        assertEquals(lqn.tasksof.get(2).get(0),4);

        assertEquals(lqn.entriesof.size(),2);
        assertEquals(lqn.entriesof.get(4).size(),1);
        assertEquals(lqn.entriesof.get(3).size(),1);
        assertEquals(lqn.entriesof.get(4).get(0),6);
        assertEquals(lqn.entriesof.get(3).get(0),5);

        assertEquals(lqn.actsof.size(),2);
        assertEquals(lqn.actsof.get(4).size(),2);
        assertEquals(lqn.actsof.get(3).size(),2);
        assertEquals(lqn.actsof.get(4).get(0),9);
        assertEquals(lqn.actsof.get(4).get(1),10);
        assertEquals(lqn.actsof.get(3).get(0),7);
        assertEquals(lqn.actsof.get(3).get(1),8);

        assertEquals(lqn.callsof.get(8).size(),1);
        assertEquals(lqn.callsof.get(8).get(0),1);

        assertTrue(lqn.hostdem.get(3).isImmediate());
        assertTrue(lqn.hostdem.get(4).isImmediate());
        assertTrue(lqn.hostdem.get(5).isImmediate());
        assertTrue(lqn.hostdem.get(6).isImmediate());
        assertTrue(lqn.hostdem.get(8).isImmediate());
        assertTrue(lqn.hostdem.get(7) instanceof Exp);
        assertTrue(lqn.hostdem.get(9) instanceof Exp);
        assertTrue(lqn.hostdem.get(10) instanceof Exp);
        assertEquals(lqn.hostdem.get(7).getMean(),1/0.625);
        assertEquals(lqn.hostdem.get(9).getMean(),1/0.2);
        assertEquals(lqn.hostdem.get(10).getMean(),1.0);

        assertEquals(lqn.think.size(),2);
        assertTrue(lqn.think.get(3) instanceof Exp);
        assertEquals(lqn.think.get(3).getMean(),1/0.01);
        assertTrue(lqn.think.get(4).isImmediate());

        assertEquals(lqn.sched.size(),4);
        assertSame(lqn.sched.get(1), SchedStrategy.PS);
        assertSame(lqn.sched.get(2), SchedStrategy.PS);
        assertSame(lqn.sched.get(3),SchedStrategy.REF);
        assertSame(lqn.sched.get(4),SchedStrategy.FCFS);

        assertEquals(lqn.schedid.getNumCols(),5);
        assertEquals(lqn.schedid.getNumRows(),1);
        assertEquals(lqn.schedid.getNonZeroLength(),4);
        assertEquals(lqn.schedid.get(0,1),6);
        assertEquals(lqn.schedid.get(0,2),6);
        assertEquals(lqn.schedid.get(0,3),14);
        assertEquals(lqn.schedid.get(0,4),1);

        assertEquals(lqn.names.size(),10);
        assertEquals(lqn.names.get(1),"P1");
        assertEquals(lqn.names.get(2),"P2");
        assertEquals(lqn.names.get(3),"T1");
        assertEquals(lqn.names.get(4),"T2");
        assertEquals(lqn.names.get(5),"E1");
        assertEquals(lqn.names.get(6),"E2");
        assertEquals(lqn.names.get(7),"A1");
        assertEquals(lqn.names.get(8),"A2");
        assertEquals(lqn.names.get(9),"A3");
        assertEquals(lqn.names.get(10),"A4");

        assertEquals(lqn.hashnames.size(),10);
        assertEquals(lqn.hashnames.get(1),"P:P1");
        assertEquals(lqn.hashnames.get(2),"P:P2");
        assertEquals(lqn.hashnames.get(3),"R:T1");
        assertEquals(lqn.hashnames.get(4),"T:T2");
        assertEquals(lqn.hashnames.get(5),"E:E1");
        assertEquals(lqn.hashnames.get(6),"E:E2");
        assertEquals(lqn.hashnames.get(7),"A:A1");
        assertEquals(lqn.hashnames.get(8),"A:A2");
        assertEquals(lqn.hashnames.get(9),"A:A3");
        assertEquals(lqn.hashnames.get(10),"A:A4");

        assertEquals(lqn.mult.getNumRows(),1);
        assertEquals(lqn.mult.getNumCols(),5);
        assertEquals(lqn.mult.getNonZeroLength(),4);
        assertEquals(lqn.mult.get(1),1);
        assertEquals(lqn.mult.get(2),1);
        assertEquals(lqn.mult.get(3),10);
        assertEquals(lqn.mult.get(4),1);

        assertEquals(lqn.repl.getNumRows(),1);
        assertEquals(lqn.repl.getNumCols(),5);
        assertEquals(lqn.repl.getNonZeroLength(),4);
        assertEquals(lqn.repl.get(1),1);
        assertEquals(lqn.repl.get(2),1);
        assertEquals(lqn.repl.get(3),1);
        assertEquals(lqn.repl.get(4),1);

        assertEquals(lqn.type.getNumRows(),1);
        assertEquals(lqn.type.getNumCols(),11);
        assertEquals(lqn.type.getNonZeroLength(),8);
        assertEquals(lqn.type.get(3),1);
        assertEquals(lqn.type.get(4),1);
        assertEquals(lqn.type.get(5),2);
        assertEquals(lqn.type.get(6),2);
        assertEquals(lqn.type.get(7),3);
        assertEquals(lqn.type.get(8),3);
        assertEquals(lqn.type.get(9),3);
        assertEquals(lqn.type.get(10),3);

        assertEquals(lqn.nitems.getNonZeroLength(),0);

        assertTrue(lqn.itemcap.isEmpty());

        assertEquals(lqn.replacement.getNonZeroLength(),0);

        assertTrue(lqn.itemproc.isEmpty());

        assertEquals(lqn.calltype.size(),1);
        assertSame(lqn.calltype.get(1), CallType.SYNC);

        assertEquals(lqn.callpair.getNumCols(),3);
        assertEquals(lqn.callpair.getNumRows(),9);
        assertEquals(lqn.callpair.getNonZeroLength(),2);
        assertEquals(lqn.callpair.get(1,1),8);
        assertEquals(lqn.callpair.get(1,2),6);

        assertEquals(lqn.callproc.size(),1);
        assertNotNull(lqn.callproc.get(1));

        assertEquals(lqn.callnames.size(),1);
        assertEquals(lqn.callnames.get(1),"A2=>E2");

        assertEquals(lqn.callhashnames.size(),1);
        assertEquals(lqn.callhashnames.get(1),"A:A2=>E:E2");

        assertEquals(lqn.actpretype.getNumRows(),1);
        assertEquals(lqn.actpretype.getNumCols(),11);
        assertEquals(lqn.actpretype.getNonZeroLength(),2);
        assertEquals(lqn.actpretype.get(0,7),1);
        assertEquals(lqn.actpretype.get(0,9),1);

        assertEquals(lqn.actposttype.getNumCols(),11);
        assertEquals(lqn.actposttype.getNumRows(),1);
        assertEquals(lqn.actposttype.getNonZeroLength(),2);
        assertEquals(lqn.actposttype.get(0,8),11);
        assertEquals(lqn.actposttype.get(0,10),11);

        assertEquals(lqn.graph.getNumRows(),11);
        assertEquals(lqn.graph.getNumCols(),11);
        assertEquals(lqn.graph.getNonZeroLength(),9);
        assertEquals(lqn.graph.get(3,1),1);
        assertEquals(lqn.graph.get(4,2),1);
        assertEquals(lqn.graph.get(3,5),1);
        assertEquals(lqn.graph.get(4,6),1);
        assertEquals(lqn.graph.get(8,6),1);
        assertEquals(lqn.graph.get(5,7),1);
        assertEquals(lqn.graph.get(7,8),1);
        assertEquals(lqn.graph.get(6,9),1);
        assertEquals(lqn.graph.get(9,10),1);

        assertEquals(lqn.parent.getNonZeroLength(),8);
        assertEquals(lqn.parent.getNumRows(),1);
        assertEquals(lqn.parent.getNumCols(),11);
        assertEquals(lqn.parent.get(0,3),1);
        assertEquals(lqn.parent.get(0,4),2);
        assertEquals(lqn.parent.get(0,5),3);
        assertEquals(lqn.parent.get(0,6),4);
        assertEquals(lqn.parent.get(0,7),3);
        assertEquals(lqn.parent.get(0,8),3);
        assertEquals(lqn.parent.get(0,9),4);
        assertEquals(lqn.parent.get(0,10),4);

        assertEquals(lqn.replygraph.getNumCols(),11);
        assertEquals(lqn.replygraph.getNumRows(),11);
        assertEquals(lqn.replygraph.getNonZeroLength(),2);
        assertEquals(lqn.replygraph.get(8,5),1);
        assertEquals(lqn.replygraph.get(10,6),1);

        assertEquals(lqn.iscache.getNonZeroLength(),0);

        assertEquals(lqn.iscaller.getNonZeroLength(),4);
        assertEquals(lqn.iscaller.getNumRows(),11);
        assertEquals(lqn.iscaller.getNumCols(),11);
        assertEquals(lqn.iscaller.get(3,4),1);
        assertEquals(lqn.iscaller.get(8,4),1);
        assertEquals(lqn.iscaller.get(3,6),1);
        assertEquals(lqn.iscaller.get(8,6),1);

        assertEquals(lqn.issynccaller.getNumCols(),11);
        assertEquals(lqn.issynccaller.getNumRows(),11);
        assertEquals(lqn.issynccaller.getNonZeroLength(),4);
        assertEquals(lqn.issynccaller.get(3,4),1);
        assertEquals(lqn.issynccaller.get(8,4),1);
        assertEquals(lqn.issynccaller.get(3,6),1);
        assertEquals(lqn.issynccaller.get(8,6),1);

        assertEquals(lqn.isasynccaller.getNonZeroLength(),0);

        assertEquals(lqn.isref.getNumRows(),1);
        assertEquals(lqn.isref.getNonZeroLength(),1);
        assertEquals(lqn.isref.getNumCols(),5);
        assertEquals(lqn.isref.get(0,3),1);
    }

    @org.junit.jupiter.api.Test
    public void testLayerStruct() throws  Exception{
        // this test is to test the network's job classes, stations and connections

        SolverLN solverLN = new SolverLN(buildModel3());

        Network network1 = solverLN.getEnsemble().get(0);

        Network network2 = solverLN.getEnsemble().get(1);

        // test layer 1
        // test numbers of nodes, classes, stations
        assertEquals(network1.getNumberOfNodes(),6);
        assertEquals(network1.getNumberOfClasses(),9);
        assertEquals(network1.getNumberOfStations(),2);

        // test classes types, populations
        for(int i=0;i<network1.getNumberOfClasses();i++) {
            assertTrue(network1.getClasses().get(i) instanceof ClosedClass);
        }
        ClosedClass temp = (ClosedClass) network1.getClasses().get(0);
        assertEquals(temp.getPopulation(),10);
        temp = (ClosedClass) network1.getClasses().get(5);
        assertEquals(temp.getPopulation(),10);
        for(int i=1;i<9;i++) {
            if(i!=5) {
                temp = (ClosedClass) network1.getClasses().get(i);
                assertEquals(temp.getPopulation(), 0);
            }
        }

        // test connection matrix
        assertEquals(network1.getConnectionMatrix().get(2,0),1);
        assertEquals(network1.getConnectionMatrix().get(4,0),1);
        assertEquals(network1.getConnectionMatrix().get(3,1),1);
        assertEquals(network1.getConnectionMatrix().get(5,1),1);
        assertEquals(network1.getConnectionMatrix().get(0,2),1);
        assertEquals(network1.getConnectionMatrix().get(0,3),1);
        assertEquals(network1.getConnectionMatrix().get(1,4),1);
        assertEquals(network1.getConnectionMatrix().get(1,5),1);
        assertEquals(network1.getConnectionMatrix().getNonZeroLength(),8);

        // test the node types
        for(int i=2;i<network1.getNumberOfNodes();i++){
            assertTrue(network1.getNodes().get(i) instanceof ClassSwitch);
        }


        // for every delay and queue test the Service Process
        assertTrue(network1.getStations().get(1) instanceof Queue);
        Queue queue = (Queue) network1.getStations().get(1);
        //assertTrue(queue.getServiceProcess(network1.getJobClass().get(0)) instanceof Disabled);
        //assertTrue(queue.getServiceProcess(network1.getJobClass().get(1)) instanceof Disabled);
        assertTrue(queue.getServiceProcess(network1.getJobClass().get(2)) instanceof Exp);
        assertTrue(queue.getServiceProcess(network1.getJobClass().get(3)) instanceof Immediate);
        assertTrue(queue.getServiceProcess(network1.getJobClass().get(4)) instanceof Immediate);
        //assertTrue(queue.getServiceProcess(network1.getJobClass().get(5)) instanceof Disabled);
        //assertTrue(queue.getServiceProcess(network1.getJobClass().get(6)) instanceof Disabled);
        assertTrue(queue.getServiceProcess(network1.getJobClass().get(7)) instanceof Exp);
        assertTrue(queue.getServiceProcess(network1.getJobClass().get(8)) instanceof Exp);

        assertTrue(compare(queue.getServiceProcess(network1.getJobClass().get(2)).getMean(),70));
        assertTrue(compare(queue.getServiceProcess(network1.getJobClass().get(7)).getMean(),5));
        assertTrue(compare(queue.getServiceProcess(network1.getJobClass().get(8)).getMean(),1));
        assertEquals(queue.getNumberOfServers(),1);

        assertTrue(network1.getStations().get(0) instanceof Delay);
        Delay delay = (Delay) network1.getStations().get(0);
        assertTrue(delay.getServiceProcess(network1.getJobClass().get(0)) instanceof Exp);
        assertTrue(delay.getServiceProcess(network1.getJobClass().get(1)) instanceof Immediate);
        //assertTrue(delay.getServiceProcess(network1.getJobClass().get(2)) instanceof Disabled);
        //assertTrue(delay.getServiceProcess(network1.getJobClass().get(3)) instanceof Disabled);
        assertTrue(delay.getServiceProcess(network1.getJobClass().get(4)) instanceof Immediate);
        assertTrue(delay.getServiceProcess(network1.getJobClass().get(5)) instanceof Immediate);
        assertTrue(delay.getServiceProcess(network1.getJobClass().get(6)) instanceof Immediate);
        //assertTrue(delay.getServiceProcess(network1.getJobClass().get(7)) instanceof Disabled);
        //assertTrue(delay.getServiceProcess(network1.getJobClass().get(8)) instanceof Disabled);

        assertTrue(compare(delay.getServiceProcess(network1.getJobClass().get(0)).getMean(),100));

        // test layer 2
        // test numbers of nodes, classes, stations
        assertEquals(network2.getNumberOfNodes(),4);
        assertEquals(network2.getNumberOfClasses(),5);
        assertEquals(network2.getNumberOfStations(),2);

        // test classes types, populations
        for(int i=0;i<network2.getNumberOfClasses();i++) {
            assertTrue(network2.getClasses().get(i) instanceof ClosedClass);
        }
        temp = (ClosedClass) network1.getClasses().get(0);
        assertEquals(temp.getPopulation(),10);
        for(int i=1;i<network2.getNumberOfClasses();i++) {
            temp = (ClosedClass) network1.getClasses().get(i);
            assertEquals(temp.getPopulation(), 0);
        }

        // test connection matrix
        assertEquals(network2.getConnectionMatrix().get(1,0),1);
        assertEquals(network2.getConnectionMatrix().get(2,0),1);
        assertEquals(network2.getConnectionMatrix().get(3,1),1);
        assertEquals(network2.getConnectionMatrix().get(0,2),1);
        assertEquals(network2.getConnectionMatrix().get(0,3),1);
        assertEquals(network2.getConnectionMatrix().getNonZeroLength(),5);

        // test the node types
        for(int i=2;i<network2.getNumberOfNodes();i++){
            assertTrue(network2.getNodes().get(i) instanceof ClassSwitch);
        }

        assertTrue(network2.getStations().get(0) instanceof Delay);
        delay = (Delay) network2.getStations().get(0);
        assertTrue(delay.getServiceProcess(network2.getJobClass().get(0)) instanceof Exp);
        assertTrue(delay.getServiceProcess(network2.getJobClass().get(1)) instanceof Immediate);
        assertTrue(delay.getServiceProcess(network2.getJobClass().get(2)) instanceof Exp);
        assertTrue(delay.getServiceProcess(network2.getJobClass().get(3)) instanceof Immediate);
        assertTrue(delay.getServiceProcess(network2.getJobClass().get(4)) instanceof Immediate);

        assertTrue(compare(delay.getServiceProcess(network2.getJobClass().get(0)).getMean(),100));
        assertTrue(compare(delay.getServiceProcess(network2.getJobClass().get(2)).getMean(),70));

        assertTrue(network2.getStations().get(1) instanceof Queue);
        queue = (Queue) network2.getStations().get(1);
        //assertTrue(queue.getServiceProcess(network2.getJobClass().get(0)) instanceof Disabled);
        //assertTrue(queue.getServiceProcess(network2.getJobClass().get(1)) instanceof Disabled);;
        //assertTrue(queue.getServiceProcess(network2.getJobClass().get(2)) instanceof Disabled);
        //assertTrue(queue.getServiceProcess(network2.getJobClass().get(3)) instanceof Disabled);
        assertTrue(queue.getServiceProcess(network2.getJobClass().get(4)) instanceof Immediate);
        assertEquals(queue.getNumberOfServers(),10);

        //test idxhash
        assertEquals(solverLN.getIdxhash().get(0), Double.NaN);
        assertEquals(solverLN.getIdxhash().get(1), 1);
        assertEquals(solverLN.getIdxhash().get(2), Double.NaN);
        assertEquals(solverLN.getIdxhash().get(3), 2);
    }

    @org.junit.jupiter.api.Test
    public void testSolveCallModel() throws  Exception{
        SolverLN solver = new SolverLN(buildModel4());
        double[] expectedQlen = {Double.NaN, Double.NaN, 12.1041, 0.0039, 0.1244, 12.1041, 0.0039, 0.1244, 12.1041, 0.0039, 0.1244};
        double[] expectedUtil = {0.9473, 0.0415, 0.9470, 0.0003, 0.0415, Double.NaN, Double.NaN, Double.NaN, 0.9470, 0.0003, 0.0415};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 0.6391, 0.2887, 0.0200, 0.6391, 0.2887, 0.0200};
        double[] expectedResidT = {Double.NaN, Double.NaN, 0.3531, 0.1884, 0.0200, Double.NaN, Double.NaN, Double.NaN, 0.3531, 0.1884, 0.0200};
        double[] expectedTput = {0.0134, 6.2189, 18.9393, 0.0134, 6.2189, 18.9393, 0.0134, 6.2189, 18.9393, 0.0134, 6.2189};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            assertTrue(compareOutput(avg.getQLen().get(idx), expectedQlen[idx]));
            assertTrue(compareOutput(avg.getUtil().get(idx), expectedUtil[idx]));
            assertTrue(compareOutput(avg.getRespT().get(idx), expectedRespT[idx]));
            assertTrue(compareOutput(avg.getResidT().get(idx), expectedResidT[idx]));
            assertTrue(compareOutput(avg.getTput().get(idx), expectedTput[idx]));
        }
    }

    @org.junit.jupiter.api.Test
    public void testSolveLoopModel() throws  Exception{
        SolverLN solver = new SolverLN(buildModel5());
        double[] expectedQlen = {Double.NaN, 1, 1, 0.1587, 0.5238, 0.3174};
        double[] expectedUtil = {1, 1, Double.NaN, 0.1587, 0.5238, 0.3174};
        double[] expectedRespT = {Double.NaN, Double.NaN, 12.6, 2, 3, 4};
        double[] expectedResidT = {Double.NaN, 12.6, Double.NaN, 2, 6.6, 4};
        double[] expectedTput = {0.07936, 0.07936, 0.07936, 0.07936, 0.1746, 0.07936};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            assertTrue(compareOutput(avg.getQLen().get(idx), expectedQlen[idx]));
            assertTrue(compareOutput(avg.getUtil().get(idx), expectedUtil[idx]));
            assertTrue(compareOutput(avg.getRespT().get(idx), expectedRespT[idx]));
            assertTrue(compareOutput(avg.getResidT().get(idx), expectedResidT[idx]));
            assertTrue(compareOutput(avg.getTput().get(idx), expectedTput[idx]));
        }
    }
    @org.junit.jupiter.api.Test
    public void testSolveOrForkOrJoinModel() throws  Exception{
        SolverLN solver = new SolverLN(buildModel6());
        double[] expectedQlen = {Double.NaN, 1, 1, 0.16529, 0.07438, 0.09917, 0.1652, 0.4958};
        double[] expectedUtil = {0.1, 0.1, Double.NaN, 0.01652, 0.00743, 0.00991, 0.01652, 0.04958};
        double[] expectedRespT = {Double.NaN, Double.NaN, 12.1, 2, 3, 4, 5, 6};
        double[] expectedResidT = {Double.NaN, 12.1, Double.NaN, 2, 0.9, 1.2, 2, 6};
        double[] expectedTput = {0.08264, 0.08264, 0.08264, 0.08264, 0.02479, 0.02479, 0.03305, 0.08264};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            assertTrue(compareOutput(avg.getQLen().get(idx), expectedQlen[idx]));
            assertTrue(compareOutput(avg.getUtil().get(idx), expectedUtil[idx]));
            assertTrue(compareOutput(avg.getRespT().get(idx), expectedRespT[idx]));
            assertTrue(compareOutput(avg.getResidT().get(idx), expectedResidT[idx]));
            assertTrue(compareOutput(avg.getTput().get(idx), expectedTput[idx]));
        }
    }
}

