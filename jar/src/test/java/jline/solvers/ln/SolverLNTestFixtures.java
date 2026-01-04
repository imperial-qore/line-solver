package jline.solvers.ln;

import jline.lang.constant.SchedStrategy;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.layered.*;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.lang.processes.DiscreteSampler;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.Arrays;

public class SolverLNTestFixtures {

    public static LayeredNetwork buildModel1() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_LQN_1");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 50, SchedStrategy.REF);
        T1.setThinkTime(new Exp(1.0 / 100.0));
        T1.on(P1);

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);

        Activity A1 = new Activity(model, "A1", new Exp(10));
        A1.on(T1);
        A1.boundTo(E1);

        Activity A2 = new Activity(model, "A2", new Exp(1.0 / 1.5));
        A2.on(T1);

        T1.addPrecedence(ActivityPrecedence.Serial("A1", "A2"));
        return model;
    }

    public static LayeredNetwork buildModel2() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_LQN_2");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 10, SchedStrategy.REF);
        T1.on(P1);
        T1.setThinkTime(new Exp(1.0 / 100));

        Task T2 = new Task(model, "T2", 1, SchedStrategy.FCFS);
        T2.on(P2);
        T2.setThinkTime(Immediate.getInstance());

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);

        Entry E2 = new Entry(model, "E2");
        E2.on(T2);

        Activity A1 = new Activity(model, "A1", new Exp(1 / 1.6));
        A1.on(T1);
        A1.boundTo(E1);

        Activity A2 = new Activity(model, "A2", Immediate.getInstance());
        A2.on(T1);
        A2.synchCall(E2);

        Activity A3 = new Activity(model, "A3", new Exp(1.0 / 5));
        A3.on(T2);
        A3.boundTo(E2);

        Activity A4 = new Activity(model, "A4", new Exp(1));
        A4.on(T2);
        A4.repliesTo(E2);

        T1.addPrecedence(ActivityPrecedence.Serial("A1", "A2"));
        T2.addPrecedence(ActivityPrecedence.Serial("A3", "A4"));
        return model;
    }

    // This is the former TestSolverLN2
    public static LayeredNetwork buildModel3() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_LQN_3");
        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 10, SchedStrategy.REF);
        T1.on(P1);
        T1.setThinkTime(new Exp(1.0 / 100));

        Task T2 = new Task(model, "T2", 10, SchedStrategy.FCFS);
        T2.on(P1);
        T2.setThinkTime(Immediate.getInstance());

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);

        Entry E2 = new Entry(model, "E2");
        E2.on(T2);

        Activity A1 = new Activity(model, "A1", new Exp(1 / 70.0));
        A1.on(T1);
        A1.boundTo(E1);

        Activity A2 = new Activity(model, "A2", Immediate.getInstance());
        A2.on(T1);
        A2.synchCall(E2);

        Activity A3 = new Activity(model, "A3", new Exp(1.0 / 5));
        A3.on(T2);
        A3.boundTo(E2);

        Activity A4 = new Activity(model, "A4", new Exp(1));
        A4.on(T2);
        A4.repliesTo(E2);

        T1.addPrecedence(ActivityPrecedence.Serial("A1", "A2"));
        T2.addPrecedence(ActivityPrecedence.Serial("A3", "A4"));
        return model;
    }

    /* Model for calls - seems to run very very slowly with LN, to be checked */
    public static LayeredNetwork buildModel4() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_LQN_4");

        Processor P1 = new Processor(model, "P1", 2, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 3, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 50, SchedStrategy.REF);
        T1.on(P1);
        T1.setThinkTime(new Exp((double) 1 / 2));
        Task T2 = new Task(model, "T2", 50, SchedStrategy.FCFS);
        T2.on(P1);
        T2.setThinkTime(new Exp((double) 1 / 3));
        Task T3 = new Task(model, "T3", 25, SchedStrategy.FCFS);
        T3.on(P2);
        T3.setThinkTime(new Exp((double) 1 / 4));

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);
        Entry E2 = new Entry(model, "E2");
        E2.on(T2);
        Entry E3 = new Entry(model, "E3");
        E3.on(T3);

        Activity A1 = new Activity(model, "AS1", new Exp(10));
        A1.on(T1);
        A1.boundTo(E1);
        A1.synchCall(E2, 1);
        Activity A2 = new Activity(model, "AS2", new Exp(20));
        A2.on(T2);
        A2.boundTo(E2);
        A2.synchCall(E3, 5);
        A2.repliesTo(E2);
        Activity A3 = new Activity(model, "AS3", new Exp(50));
        A3.on(T3);
        A3.boundTo(E3);
        A3.repliesTo(E3);

        return model;
    }

    /* Model for Loop Precedences */
    public static LayeredNetwork buildModel5() throws Exception {
        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", 10, SchedStrategy.INF);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF);
        T1.on(P1);
        T1.setThinkTime(Immediate.getInstance());

        Entry E1 = new Entry(model, "Entry");
        E1.on(T1);

        Activity A1 = new Activity(model, "A1", new Exp(0.5));
        A1.on(T1);
        A1.boundTo(E1);
        Activity A2 = new Activity(model, "A2", new Exp(0.333333));
        A2.on(T1);
        Activity A3 = new Activity(model, "A3", new Exp(0.25));
        A3.on(T1);


        // Loop Activity Precedence
        ArrayList<String> precActs = new ArrayList<String>();
        precActs.add("A2");
        precActs.add("A3");
        T1.addPrecedence(ActivityPrecedence.Loop("A1", precActs, Matrix.singleton(2.2)));
        return model;
    }

    public static LayeredNetwork buildModel6() throws Exception {
        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", 10, SchedStrategy.FCFS);
        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF).on(P1).setThinkTime(Immediate.getInstance());
        Entry E1 = new Entry(model, "Entry").on(T1);

        Activity A1 = new Activity(model, "B1", Exp.fitMean(2)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "B2", Exp.fitMean(3)).on(T1);
        Activity A3 = new Activity(model, "B3", Exp.fitMean(4)).on(T1);
        Activity A4 = new Activity(model, "B4", Exp.fitMean(6)).on(T1);
        Activity A5 = new Activity(model, "B5", Exp.fitMean(6)).on(T1);

        // OrFork Activity Precedence
        ArrayList<String> precActs = new ArrayList<String>();
        Matrix probs = new Matrix(1, 3);
        precActs.add("B2");
        precActs.add("B3");
        precActs.add("B4");
        probs.set(0, 0, 0.3);
        probs.set(0, 1, 0.3);
        probs.set(0, 2, 0.4);
        T1.addPrecedence(ActivityPrecedence.OrFork("B1", precActs, probs));

        // OrJoin Activity Precedence
        ArrayList<String> precActs2 = new ArrayList<String>();
        precActs2.add("B2");
        precActs2.add("B3");
        precActs2.add("B4");
        T1.addPrecedence(ActivityPrecedence.OrJoin(precActs2, "B5"));

        return model;
    }

    public static LayeredNetwork test_activityGraph_seq() throws Exception {

        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF).on(P1).setThinkTime(new Exp(0.01));
        Task T2 = new Task(model, "T2", 1, SchedStrategy.FCFS).on(P2).setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);

        Activity A11 = new Activity(model, "A11", new Exp(0.625)).on(T1).boundTo(E1).synchCall(E2, 1.0);

        Activity A21 = new Activity(model, "A21", new Exp(0.2)).on(T2).boundTo(E2);
        Activity A22 = new Activity(model, "A22", new Exp(1)).on(T2).repliesTo(E2);

        T2.addPrecedence(ActivityPrecedence.Serial("A21", "A22"));

        return model;
    }

    public static LayeredNetwork test_activityGraph_call() throws Exception {

        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);
        Processor P3 = new Processor(model, "P3", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF).on(P1);
        T1.setThinkTime(new Exp(0.01));
        Task T2 = new Task(model, "T2", 1, SchedStrategy.FCFS).on(P2);
        T2.setThinkTime(new Immediate());
        Task T3 = new Task(model, "T3", 1, SchedStrategy.FCFS).on(P3);
        T3.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);
        Entry E3 = new Entry(model, "E3").on(T3);

        Activity A1 = new Activity(model, "A11", new Exp(0.625)).on(T1);
        A1.boundTo(E1);
        A1.synchCall(E2, 1);
        Activity A2 = new Activity(model, "A21", new Exp(0.2)).on(T2);
        A2.boundTo(E2);
        A2.synchCall(E3, 1);
        A2.repliesTo(E2);
        Activity A3 = new Activity(model, "A31", new Exp(0.125)).on(T3);
        A3.boundTo(E3);
        A3.repliesTo(E3);

        return model;
    }

    public static LayeredNetwork test_activityGraph_call_seq() throws Exception {

        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);
        Processor P3 = new Processor(model, "P3", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF).on(P1);
        T1.setThinkTime(new Exp(0.01));
        Task T2 = new Task(model, "T2", 1, SchedStrategy.FCFS).on(P2);
        T2.setThinkTime(new Immediate());
        Task T3 = new Task(model, "T3", 1, SchedStrategy.FCFS).on(P3);
        T3.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);
        Entry E3 = new Entry(model, "E3").on(T3);

        Activity A1 = new Activity(model, "A11", new Exp(0.625)).on(T1);
        A1.boundTo(E1);
        A1.synchCall(E2, 1).synchCall(E3, 1);
        Activity A2 = new Activity(model, "A21", new Exp(0.2)).on(T2);
        A2.boundTo(E2);
        Activity A3 = new Activity(model, "A22", new Exp(1)).on(T2);
        A3.synchCall(E3, 1);
        A3.repliesTo(E2);
        Activity A4 = new Activity(model, "A31", new Exp(0.125)).on(T3);
        A4.boundTo(E3);
        A4.repliesTo(E3);

        T2.addPrecedence(ActivityPrecedence.Serial("A21", "A22"));

        return model;
    }

    public static LayeredNetwork test_activityGraph_call_seq_disconnected() throws Exception {
        // this model has an unreachable connected component (P3,T3,E3,A4) isolated from the reference task
        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);
        Processor P3 = new Processor(model, "P3", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF).on(P1);
        T1.setThinkTime(new Exp(0.01));
        Task T2 = new Task(model, "T2", 1, SchedStrategy.FCFS).on(P2);
        T2.setThinkTime(new Immediate());
        Task T3 = new Task(model, "T3", 1, SchedStrategy.FCFS).on(P3);
        T3.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);
        Entry E3 = new Entry(model, "E3").on(T3);

        Activity A1 = new Activity(model, "A11", new Exp(0.625)).on(T1);
        A1.boundTo(E1);
        A1.synchCall(E2, 1);
        Activity A2 = new Activity(model, "A21", new Exp(0.2)).on(T2);
        A2.boundTo(E2);
        Activity A3 = new Activity(model, "A22", new Exp(1)).on(T2);
        A3.repliesTo(E2);
        Activity A4 = new Activity(model, "A31", new Exp(0.125)).on(T3);
        A4.boundTo(E3);
        A4.repliesTo(E3);

        T2.addPrecedence(ActivityPrecedence.Serial("A21", "A22"));

        return model;
    }

    public static LayeredNetwork test_activityGraph_call_and() throws Exception {
        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", Integer.MAX_VALUE, SchedStrategy.INF);
        Processor P3 = new Processor(model, "P3", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF).on(P1);
        T1.setThinkTime(new Immediate());
        Task T2 = new Task(model, "T2", Integer.MAX_VALUE, SchedStrategy.INF).on(P2);
        T2.setThinkTime(new Immediate());
        Task T3 = new Task(model, "T3", 1, SchedStrategy.FCFS).on(P3);
        T3.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);
        Entry E3 = new Entry(model, "E3").on(T3);

        Activity A1 = new Activity(model, "A11", new Immediate()).on(T1);
        A1.boundTo(E1);
        A1.synchCall(E2, 1);
        Activity A2 = new Activity(model, "A20", new Immediate()).on(T2);
        A2.boundTo(E2);
        A2.synchCall(E3, 1);
        Activity A3 = new Activity(model, "A21a", new Exp(0.5)).on(T2);
        Activity A4 = new Activity(model, "A21b", new Exp(1)).on(T2);
        Activity A5 = new Activity(model, "A22", new Immediate()).on(T2);
        A5.repliesTo(E2);
        Activity A6 = new Activity(model, "A31", new Exp(0.125)).on(T3);
        A6.boundTo(E3);
        A6.repliesTo(E3);

        // AndFork Activity Precedence
        ArrayList<String> postActs = new ArrayList<String>();
        postActs.add("A21a");
        postActs.add("A21b");
        T2.addPrecedence(ActivityPrecedence.AndFork("A20", postActs));

        // AndJoin Activity Precedence
        ArrayList<String> precActs = new ArrayList<String>();
        precActs.add("A21a");
        precActs.add("A21b");
        T2.addPrecedence(ActivityPrecedence.AndJoin(precActs, "A22"));

        return model;
    }

    public static LayeredNetwork test_activityGraph_call_or() throws Exception {

        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", Integer.MAX_VALUE, SchedStrategy.INF);
        Processor P3 = new Processor(model, "P3", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF).on(P1);
        T1.setThinkTime(new Immediate());
        Task T2 = new Task(model, "T2", Integer.MAX_VALUE, SchedStrategy.INF).on(P2);
        T2.setThinkTime(new Immediate());
        Task T3 = new Task(model, "T3", 1, SchedStrategy.FCFS).on(P3);
        T3.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);
        Entry E3 = new Entry(model, "E3").on(T3);

        Activity A1 = new Activity(model, "A11", new Exp(0.625)).on(T1);
        A1.boundTo(E1);
        A1.synchCall(E2, 1);
        Activity A2 = new Activity(model, "A20", new Immediate()).on(T2);
        A2.boundTo(E2);
        Activity A3 = new Activity(model, "A21a", new Exp(0.5)).on(T2);
        Activity A4 = new Activity(model, "A21b", new Exp(0.333333)).on(T2);
        Activity A5 = new Activity(model, "A22", new Exp(0.01)).on(T2);
        A5.repliesTo(E2);
        Activity A6 = new Activity(model, "A31", new Exp(0.125)).on(T3);
        A6.boundTo(E3);
        A6.repliesTo(E3);


        // OrFork Activity Precedence
        ArrayList<String> precActs = new ArrayList<String>();
        Matrix probs = new Matrix(1, 2);
        precActs.add("A21a");
        precActs.add("A21b");
        probs.set(0, 0, 0.5);
        probs.set(0, 1, 0.5);
        T2.addPrecedence(ActivityPrecedence.OrFork("A20", precActs, probs));

        // OrJoin Activity Precedence
        precActs = new ArrayList<String>();
        precActs.add("A21a");
        precActs.add("A21b");
        T2.addPrecedence(ActivityPrecedence.OrJoin(precActs, "A22"));
        return model;
    }

    public static LayeredNetwork test_activityGraph_loop() throws Exception {

        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", Integer.MAX_VALUE, SchedStrategy.INF);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF).on(P1);
        T1.setThinkTime(new Exp(0.01));
        Task T2 = new Task(model, "T2", Integer.MAX_VALUE, SchedStrategy.INF).on(P2);
        T2.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);

        Activity A1 = new Activity(model, "A11", new Immediate()).on(T1);
        A1.boundTo(E1);
        A1.synchCall(E2, 1);
        Activity A2 = new Activity(model, "A20", new Immediate()).on(T2);
        A2.boundTo(E2);
        Activity A3 = new Activity(model, "A21", new Exp(0.333333)).on(T2);
        Activity A4 = new Activity(model, "A22", new Exp(1)).on(T2);
        A4.repliesTo(E2);


        // Loop Activity Precedence
        ArrayList<String> precActs = new ArrayList<String>();
        precActs.add("A21");
        precActs.add("A22");
        T2.addPrecedence(ActivityPrecedence.Loop("A20", precActs, Matrix.singleton(3)));

        return model;
    }

    public static LayeredNetwork test_activityGraph_or() throws Exception {

        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", Integer.MAX_VALUE, SchedStrategy.INF);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF).on(P1);
        T1.setThinkTime(new Immediate());
        Task T2 = new Task(model, "T2", Integer.MAX_VALUE, SchedStrategy.INF).on(P2);
        T2.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);

        Activity A1 = new Activity(model, "A11", new Exp(0.625)).on(T1);
        A1.boundTo(E1);
        A1.synchCall(E2, 1);
        Activity A2 = new Activity(model, "A20", new Immediate()).on(T2);
        A2.boundTo(E2);
        Activity A3 = new Activity(model, "A21a", new Exp(0.5)).on(T2);
        Activity A4 = new Activity(model, "A21b", new Exp(0.333333)).on(T2);
        Activity A5 = new Activity(model, "A22", new Exp(0.01)).on(T2);
        A5.repliesTo(E2);

        // OrFork Activity Precedence
        ArrayList<String> precActs = new ArrayList<String>();
        Matrix probs = new Matrix(1, 2);
        precActs.add("A21a");
        precActs.add("A21b");
        probs.set(0, 0, 0.5);
        probs.set(0, 1, 0.5);
        T2.addPrecedence(ActivityPrecedence.OrFork("A20", precActs, probs));

        // OrJoin Activity Precedence
        precActs = new ArrayList<String>();
        precActs.add("A21a");
        precActs.add("A21b");
        T2.addPrecedence(ActivityPrecedence.OrJoin(precActs, "A22"));

        return model;
    }

    public static LayeredNetwork test_activityGraph_and() throws Exception {

        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", Integer.MAX_VALUE, SchedStrategy.INF);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF).on(P1);
        T1.setThinkTime(new Immediate());
        Task T2 = new Task(model, "T2", Integer.MAX_VALUE, SchedStrategy.INF).on(P2);
        T2.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);

        Activity A1 = new Activity(model, "A11", new Immediate()).on(T1);
        A1.boundTo(E1);
        A1.synchCall(E2, 1);
        Activity A2 = new Activity(model, "A20", new Immediate()).on(T2);
        A2.boundTo(E2);
        Activity A3 = new Activity(model, "A21a", new Exp(0.5)).on(T2);
        Activity A4 = new Activity(model, "A21b", new Exp(1)).on(T2);
        Activity A5 = new Activity(model, "A21c", new Exp(1)).on(T2);
        Activity A6 = new Activity(model, "A22", new Immediate()).on(T2);
        A6.repliesTo(E2);

        // AndFork Activity Precedence
        ArrayList<String> postActs = new ArrayList<String>();
        postActs.add("A21a");
        postActs.add("A21b");
        postActs.add("A21c");
        T2.addPrecedence(ActivityPrecedence.AndFork("A20", postActs));

        // AndJoin Activity Precedence
        ArrayList<String> precActs = new ArrayList<String>();
        precActs.add("A21a");
        precActs.add("A21b");
        precActs.add("A21c");
        T2.addPrecedence(ActivityPrecedence.AndJoin(precActs, "A22"));

        return model;
    }

    public static LayeredNetwork buildCacheModel() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_cache_LQN");

        // Create processors
        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor PC = new Processor(model, "PC", 1, SchedStrategy.PS);

        // Create client task
        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF);
        T1.on(P1);

        // Create cache task with proper parameters
        int totalitems = 4;
        int cachecapacity = 2;
        CacheTask C2 = new CacheTask(model, "C2", totalitems, cachecapacity, ReplacementStrategy.RR, 1, SchedStrategy.FCFS);
        C2.on(PC);

        // Create entries
        Entry E1 = new Entry(model, "E1");
        E1.on(T1);

        // Create ItemEntry with popularity distribution
        double[] popularityProbs = new double[totalitems];
        for (int i = 0; i < totalitems; i++) {
            popularityProbs[i] = 1.0 / totalitems;
        }
        Matrix pMatrix = new Matrix(popularityProbs);
        DiscreteSampler pAccess = new DiscreteSampler(pMatrix);
        ItemEntry I2 = new ItemEntry(model, "I2", totalitems, pAccess);
        I2.on(C2);

        // Create activities
        Activity A1 = new Activity(model, "A1", Immediate.getInstance());
        A1.on(T1);
        A1.boundTo(E1);
        A1.synchCall(I2, 1);

        Activity AC2 = new Activity(model, "AC2", Immediate.getInstance());
        AC2.on(C2);
        AC2.boundTo(I2);

        Activity AC2h = new Activity(model, "AC2h", new Exp(1.0));
        AC2h.on(C2);
        AC2h.repliesTo(I2);

        Activity AC2m = new Activity(model, "AC2m", new Exp(0.5));
        AC2m.on(C2);
        AC2m.repliesTo(I2);

        // Add cache precedence
        C2.addPrecedence(ActivityPrecedence.CacheAccess(AC2, Arrays.asList(AC2h, AC2m)));

        return model;
    }


}

