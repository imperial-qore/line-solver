/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.VerboseLevel;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.layered.*;
import jline.lang.processes.APH;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import java.net.URL;
import jline.solvers.SolverOptions;
import jline.solvers.ln.LN;
import jline.solvers.lqns.LQNS;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Examples of layered networks
 */
public class LayeredModel {

    /**
     * Basic layered network with multiple processors and synchronous calls.
     * <p>
     * Features:
     * - Two PS processors P1 (2 cores) and P2 (3 cores)
     * - Three tasks: T1 (50 jobs, REF), T2 (50 jobs, FCFS), T3 (25 jobs, FCFS)
     * - Synchronous call chain: T1 -> E2 (1 call), T2 -> E3 (5 calls)
     * - Demonstrates basic multi-tier application modeling
     *
     * @return configured basic layered network model
     * @throws Exception if model creation fails
     */
    public static LayeredNetwork lqn_basic() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_LQN_4");

        Processor P1 = new Processor(model, "P1", 2, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 3, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 50, SchedStrategy.REF).on(P1).setThinkTime(new Exp(1.0 / 2));
        Task T2 = new Task(model, "T2", 50, SchedStrategy.FCFS).on(P1).setThinkTime(new Exp(1.0 / 3));
        Task T3 = new Task(model, "T3", 25, SchedStrategy.FCFS).on(P2).setThinkTime(new Exp(1.0 / 4));

        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);
        Entry E3 = new Entry(model, "E3").on(T3);

        Activity A1 = new Activity(model, "AS1", new Exp(10.0)).on(T1);
        A1.boundTo(E1);
        A1.synchCall(E2, 1);
        Activity A2 = new Activity(model, "AS2", new Exp(20.0)).on(T2);
        A2.boundTo(E2);
        A2.synchCall(E3, 5);
        A2.repliesTo(E2);
        Activity A3 = new Activity(model, "AS3", new Exp(50.0)).on(T3);
        A3.boundTo(E3);
        A3.repliesTo(E3);

        return model;
    }

    /**
     * Basic layered network with two processors and synchronous calls.
     * <p>
     * Features:
     * - Two PS processors P1 and P2 with single server each
     * - Task T1 (10 jobs, REF) with exponential think time (0.01)
     * - Task T2 (1 job, FCFS) with immediate think time
     * - Activities with serial precedence and synchronous calls
     * - Demonstrates basic layered network modeling
     *
     * @return configured layered network model
     * @throws Exception if model creation fails
     */
    public static LayeredNetwork lqn_serial() throws Exception {
        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 10, SchedStrategy.REF).on(P1).setThinkTime(Exp.fitMean(100.0));
        Task T2 = new Task(model, "T2", 1, SchedStrategy.FCFS).on(P2).setThinkTime(Immediate.getInstance());

        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);

        Activity A1 = new Activity(model, "AS1", Exp.fitMean(1.6)).on(T1);
        A1.boundTo(E1);
        Activity A2 = new Activity(model, "AS2", Immediate.getInstance()).on(T1);
        A2.synchCall(E2, 1);
        Activity A3 = new Activity(model, "AS3", Exp.fitMean(5.0)).on(T2);
        A3.boundTo(E2);
        Activity A4 = new Activity(model, "AS4", Exp.fitMean(1.0)).on(T2);
        A4.repliesTo(E2);

        T1.addPrecedence(ActivityPrecedence.Serial("AS1", "AS2"));
        T2.addPrecedence(ActivityPrecedence.Serial("AS3", "AS4"));

        // Model solution
//        LN solver = new LN(model, SolverType.MVA);
//        solver.getEnsembleAvg();
        return model;
    }

    /**
     * Layered network with infinite servers and APH service.
     * <p>
     * Features:
     * - Two infinite capacity processors (P1, P2)
     * - Task T1 (1 job, REF) with Erlang think time
     * - Task T2 (infinite jobs, INF) with immediate think time
     * - APH service distribution with high variability (SCV=10)
     * - Synchronous call with multiplicity 3
     *
     * @return configured layered network with infinite servers
     * @throws Exception if model creation fails
     */
    public static LayeredNetwork lqn_multi_solvers() throws Exception {

        LayeredNetwork model = new LayeredNetwork("LQN1");

        Processor P1 = new Processor(model, "P1", Integer.MAX_VALUE, SchedStrategy.INF);
        Processor P2 = new Processor(model, "P2", Integer.MAX_VALUE, SchedStrategy.INF);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF).on(P1);
        T1.setThinkTime(Erlang.fitMeanAndOrder(0.0001, 2));
        Task T2 = new Task(model, "T2", Integer.MAX_VALUE, SchedStrategy.INF).on(P2);
        T2.setThinkTime(Immediate.getInstance());

        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);

        Activity A1 = new Activity(model, "A1", new Exp(1)).on(T1);
        A1.boundTo(E1);
        A1.synchCall(E2, 3);
        Activity A2 = new Activity(model, "A2", APH.fitMeanAndSCV(1, 10)).on(T2);
        A2.boundTo(E2);
        A2.repliesTo(E2);


        // Model solution
//        LN solver = new LN(model, SolverType.MVA);
//        solver.getEnsembleAvg();
        return model;
    }

    /**
     * Layered network with function task demonstrating FaaS/serverless.
     * <p>
     * Features:
     * - 2 processors (P1 infinite, P2 with 4 servers)
     * - FunctionTask demonstrating cold start and keep-alive behavior
     * - Setup time (cold start penalty) and delay-off time (warm duration)
     * - Synchronous calls between client and function layers
     * - Representative of serverless/FaaS architectures
     *
     * @return configured layered network with function task
     * @throws Exception if model creation fails
     */
    public static LayeredNetwork lqn_function() throws Exception {

        LayeredNetwork model = new LayeredNetwork("faas_test_example");

        // Definition of processors, tasks and entries
        Processor P1 = new Processor(model, "P1", Integer.MAX_VALUE, SchedStrategy.INF);
        Processor P2 = new Processor(model, "P2", 4, SchedStrategy.FCFS);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "E1").on(T1);

        //Task T2 = new Task(model, "T2", 1, SchedStrategy.FCFS).on(P2); // Commented alternative
        FunctionTask T2 = new FunctionTask(model, "F2", 6, SchedStrategy.FCFS).on(P2).setThinkTime(Exp.fitMean(8.0));
        T2.setSetupTime(new Exp(1.0));      // Cold start time
        T2.setDelayOffTime(new Exp(2.0));   // Time before function instance is removed

        Entry E2 = new Entry(model, "E2").on(T2);

        // T3 = Task(model, 'T3', 1, SchedStrategy.FCFS).on(P2); // Commented alternative
        // E3 = Entry(model, 'E3').on(T3); // Commented alternative

        // Definition of activities
        Activity A1 = new Activity(model, "A1", new Exp(1.0)).on(T1).boundTo(E1).synchCall(E2, 1.0);
        Activity A2 = new Activity(model, "A2", new Exp(3.0)).on(T2).boundTo(E2).repliesTo(E2);
        // A3 = Activity(model, 'A3', Exp(5.0)).on(T3).boundTo(E3).repliesTo(E3); // Commented alternative

        return model;
    }

    /**
     * Layered network with multiple entries and synchronous calls.
     * <p>
     * Features:
     * - Two PS processors with high task populations
     * - Task T1 (100 jobs, REF) with Erlang think time
     * - Task T2 (100 jobs, INF) with multiple entries (E2, E3)
     * - Activity A1 makes synchronous calls to both E2 and E3
     * - Serial precedence pattern in T2 activities
     * - Demonstrates concurrent service requests
     *
     * @return configured multi-entry layered network
     * @throws Exception if model creation fails
     */
    public static LayeredNetwork lqn_twotasks() throws Exception {

        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Task T1 = new Task(model, "T1", 100, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "E1").on(T1);

        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);
        Task T2 = new Task(model, "T2", Integer.MAX_VALUE, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(model, "E2").on(T2);
        Entry E3 = new Entry(model, "E3").on(T2);

        T1.setThinkTime(Erlang.fitMeanAndOrder(10, 1));

        Activity A1 = new Activity(model, "A1", new Exp(1)).on(T1).boundTo(E1).synchCall(E2).synchCall(E3, 1);
        Activity A20 = new Activity(model, "A20", new Exp(1)).on(T2).boundTo(E2);
        Activity A21 = new Activity(model, "A21", new Exp(1)).on(T2);
        Activity A22 = new Activity(model, "A22", new Exp(1)).on(T2).repliesTo(E2);

        T2.addPrecedence(ActivityPrecedence.Serial(A20, A21, A22));

        Activity A5 = new Activity(model, "A3", new Exp(1)).on(T2);
        A5.boundTo(E3);
        A5.repliesTo(E3);


        // Model solution
//        LN solver = new LN(model, SolverType.MVA);
//        solver.getEnsembleAvg();
        return model;
    }

    /**
     * Complex layered network with fork-join patterns and multiple precedence types.
     * <p>
     * Features:
     * - 7 processors with varying capacities and scheduling strategies
     * - 7 tasks with different multiplicity and reference patterns
     * - 16 entries across all tasks
     * - OrFork, AndFork, OrJoin, AndJoin activity precedence patterns
     * - Demonstrates advanced layered network control flow
     * - Complex synchronous call patterns between layers
     *
     * @return configured fork-join layered network
     * @throws Exception if model creation fails
     */
    public static LayeredNetwork lqn_bpmn() throws Exception {

        //URI fileURI = GettingStarted.class.getResource("/lqn_bpmn.xml").toURI();
        //Mode model = LayeredNetwork.parseXML(fileURI.getRawPath());
        LayeredNetwork model = new LayeredNetwork("lqn_bpmn.xml");

        Processor P1 = new Processor(model, "R1_Processor", 100, SchedStrategy.FCFS);
        Processor P2 = new Processor(model, "R2_Processor", Integer.MAX_VALUE, SchedStrategy.INF);
        Processor P3 = new Processor(model, "R3_Processor", 2, SchedStrategy.FCFS);
        Processor P4 = new Processor(model, "R1A_Processor", 7, SchedStrategy.FCFS);
        Processor P5 = new Processor(model, "R1B_Processor", 3, SchedStrategy.FCFS);
        Processor P6 = new Processor(model, "R2A_Processor", 4, SchedStrategy.FCFS);
        Processor P7 = new Processor(model, "R2B_Processor", 5, SchedStrategy.FCFS);

        Task T1 = new Task(model, "R1_Task", 100, SchedStrategy.REF).on(P1);
        T1.setThinkTime(Exp.fitMean(20.0));
        Task T2 = new Task(model, "R2_Task", Integer.MAX_VALUE, SchedStrategy.INF).on(P2);
        T2.setThinkTime(Immediate.getInstance());
        Task T3 = new Task(model, "R3_Task", 2, SchedStrategy.FCFS).on(P3);
        T3.setThinkTime(Immediate.getInstance());
        Task T4 = new Task(model, "R1A_Task", 7, SchedStrategy.FCFS).on(P4);
        T4.setThinkTime(Immediate.getInstance());
        Task T5 = new Task(model, "R1B_Task", 3, SchedStrategy.FCFS).on(P5);
        T5.setThinkTime(Immediate.getInstance());
        Task T6 = new Task(model, "R2A_Task", 4, SchedStrategy.FCFS).on(P6);
        T6.setThinkTime(Immediate.getInstance());
        Task T7 = new Task(model, "R2B_Task", 5, SchedStrategy.FCFS).on(P7);
        T7.setThinkTime(Immediate.getInstance());

        Entry E1 = new Entry(model, "R1_Ref_Entry").on(T1);

        Entry E2 = new Entry(model, "R2_Synch_A2_Entry").on(T2);
        Entry E3 = new Entry(model, "R2_Synch_A5_Entry").on(T2);

        Entry E4 = new Entry(model, "R3_Synch_A9_Entry").on(T3);

        Entry E5 = new Entry(model, "R1A_Synch_A1_Entry").on(T4);
        Entry E6 = new Entry(model, "R1A_Synch_A2_Entry").on(T4);
        Entry E7 = new Entry(model, "R1A_Synch_A3_Entry").on(T4);

        Entry E8 = new Entry(model, "R1B_Synch_A4_Entry").on(T5);
        Entry E9 = new Entry(model, "R1B_Synch_A5_Entry").on(T5);
        Entry E10 = new Entry(model, "R1B_Synch_A6_Entry").on(T5);

        Entry E11 = new Entry(model, "R2A_Synch_A7_Entry").on(T6);
        Entry E12 = new Entry(model, "R2A_Synch_A8_Entry").on(T6);
        Entry E13 = new Entry(model, "R2A_Synch_A11_Entry").on(T6);

        Entry E14 = new Entry(model, "R2B_Synch_A9_Entry").on(T7);
        Entry E15 = new Entry(model, "R2B_Synch_A10_Entry").on(T7);
        Entry E16 = new Entry(model, "R2B_Synch_A12_Entry").on(T7);

        Activity A1 = new Activity(model, "A1_Empty", Immediate.getInstance()).on(T1);
        A1.boundTo(E1);
        A1.synchCall(E5, 1);
        Activity A2 = new Activity(model, "A2_Empty", Immediate.getInstance()).on(T1);
        A2.synchCall(E6, 1);
        Activity A3 = new Activity(model, "A5_Empty", Immediate.getInstance()).on(T1);
        A3.synchCall(E9, 1);
        Activity A4 = new Activity(model, "A6_Empty", Immediate.getInstance()).on(T1);
        A4.synchCall(E10, 1);
        Activity A5 = new Activity(model, "A3_Empty", Immediate.getInstance()).on(T1);
        A5.synchCall(E7, 1);
        Activity A6 = new Activity(model, "A4_Empty", Immediate.getInstance()).on(T1);
        A6.synchCall(E8, 1);

        Activity A7 = new Activity(model, "E4_Empty", Immediate.getInstance()).on(T2);
        A7.boundTo(E2);
        Activity A8 = new Activity(model, "A7_Empty", Immediate.getInstance()).on(T2);
        A8.synchCall(E11, 1);
        Activity A9 = new Activity(model, "A8_Empty", Immediate.getInstance()).on(T2);
        A9.synchCall(E12, 1);
        Activity A10 = new Activity(model, "A9_Empty", Immediate.getInstance()).on(T2);
        A10.synchCall(E14, 1);
        Activity A11 = new Activity(model, "A11_Empty", Immediate.getInstance()).on(T2);
        A11.synchCall(E13, 1);
        A11.repliesTo(E2);
        Activity A12 = new Activity(model, "A12_Empty", Immediate.getInstance()).on(T2);
        A12.boundTo(E3);
        A12.synchCall(E16, 1);
        A12.repliesTo(E3);
        Activity A13 = new Activity(model, "A10_Empty", Immediate.getInstance()).on(T2);
        A13.synchCall(E15, 1);

        Activity A14 = new Activity(model, "A13", Exp.fitMean(10.0)).on(T3);
        A14.boundTo(E4);
        A14.repliesTo(E4);
        Activity A15 = new Activity(model, "A1", Exp.fitMean(7.0)).on(T4);
        A15.boundTo(E5);
        A15.repliesTo(E5);
        Activity A16 = new Activity(model, "A2", Exp.fitMean(4.0)).on(T4);
        A16.boundTo(E6);
        Activity A17 = new Activity(model, "A3", Exp.fitMean(5.0)).on(T4);
        A17.boundTo(E7);
        A17.repliesTo(E7);
        Activity A18 = new Activity(model, "A2_Res_Empty", Immediate.getInstance()).on(T4);
        A18.synchCall(E2, 1);
        A18.repliesTo(E6);
        Activity A19 = new Activity(model, "A4", Exp.fitMean(8.0)).on(T5);
        A19.boundTo(E8);
        A19.repliesTo(E8);
        Activity A20 = new Activity(model, "A5", Exp.fitMean(4.0)).on(T5);
        A20.boundTo(E9);
        Activity A21 = new Activity(model, "A6", Exp.fitMean(6.0)).on(T5);
        A21.boundTo(E10);
        A21.repliesTo(E10);
        Activity A22 = new Activity(model, "A5_Res_Empty", Immediate.getInstance()).on(T5);
        A22.synchCall(E3, 1);
        A22.repliesTo(E9);
        Activity A23 = new Activity(model, "A7", Exp.fitMean(6.0)).on(T6);
        A23.boundTo(E11);
        A23.repliesTo(E11);
        Activity A24 = new Activity(model, "A8", Exp.fitMean(8.0)).on(T6);
        A24.boundTo(E12);
        A24.repliesTo(E12);
        Activity A25 = new Activity(model, "A11", Exp.fitMean(4.0)).on(T6);
        A25.boundTo(E13);
        A25.repliesTo(E13);
        Activity A26 = new Activity(model, "A9", Exp.fitMean(4.0)).on(T7);
        A26.boundTo(E14);
        Activity A27 = new Activity(model, "A10", Exp.fitMean(6.0)).on(T7);
        A27.boundTo(E15);
        A27.repliesTo(E15);
        Activity A28 = new Activity(model, "A12", Exp.fitMean(8.0)).on(T7);
        A28.boundTo(E16);
        A28.repliesTo(E16);
        Activity A29 = new Activity(model, "A9_Res_Empty", Immediate.getInstance()).on(T7);
        A29.synchCall(E4, 1);
        A29.repliesTo(E14);

        T1.addPrecedence(ActivityPrecedence.Serial("A1_Empty", "A2_Empty"));
        T1.addPrecedence(ActivityPrecedence.Serial("A5_Empty", "A6_Empty"));
        T2.addPrecedence(ActivityPrecedence.Serial("E4_Empty", "A7_Empty"));
        T2.addPrecedence(ActivityPrecedence.Serial("A9_Empty", "A10_Empty"));
        T4.addPrecedence(ActivityPrecedence.Serial("A2", "A2_Res_Empty"));
        T5.addPrecedence(ActivityPrecedence.Serial("A5", "A5_Res_Empty"));
        T7.addPrecedence(ActivityPrecedence.Serial("A9", "A9_Res_Empty"));

        // OrFork Activity Precedence
        ArrayList<String> precActs = new ArrayList<String>();
        Matrix probs = new Matrix(1, 2);
        precActs.add("A3_Empty");
        precActs.add("A4_Empty");
        probs.set(0, 0, 0.6);
        probs.set(0, 1, 0.4);
        T1.addPrecedence(ActivityPrecedence.OrFork("A2_Empty", precActs, probs));

        // AndFork Activity Precedence
        ArrayList<String> postActs = new ArrayList<String>();
        postActs.add("A8_Empty");
        postActs.add("A9_Empty");
        T2.addPrecedence(ActivityPrecedence.AndFork("A7_Empty", postActs));

        // OrJoin Activity Precedence
        ArrayList<String> precActs2 = new ArrayList<String>();
        precActs2.add("A3_Empty");
        precActs2.add("A4_Empty");
        T1.addPrecedence(ActivityPrecedence.OrJoin(precActs2, "A5_Empty"));

        // AndJoin Activity Precedence
        ArrayList<String> precActs3 = new ArrayList<String>();
        precActs3.add("A8_Empty");
        precActs3.add("A10_Empty");
        T2.addPrecedence(ActivityPrecedence.AndJoin(precActs3, "A11_Empty"));

        // Model solution
        //LN solver = new LN(model, SolverType.MVA);
        //solver.getEnsembleAvg();
        return model;
    }

    /**
     * Layered network demonstrating loop and fork-join precedence patterns.
     * <p>
     * Features:
     * - 3 processors: P1 (INF), P2 (INF), P3 (5 servers, PS)
     * - Task T1 with loop precedence (3 iterations)
     * - Task T2 with AndFork/AndJoin patterns
     * - Task T3 with OrFork/OrJoin patterns (30%, 30%, 40% probabilities)
     * - Nested synchronous calls between tasks
     * - Demonstrates complex control flow in layered networks
     *
     * @return configured layered network with loops and forks
     * @throws Exception if model creation fails
     */
    public static LayeredNetwork lqn_workflows() throws Exception {
        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", Integer.MAX_VALUE, SchedStrategy.INF);
        Processor P2 = new Processor(model, "P2", Integer.MAX_VALUE, SchedStrategy.INF);
        Processor P3 = new Processor(model, "P3", 5, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF).on(P1);
        T1.setThinkTime(Immediate.getInstance());
        Task T2 = new Task(model, "T2", 1, SchedStrategy.INF).on(P2);
        T2.setThinkTime(Immediate.getInstance());
        Task T3 = new Task(model, "T3", 20, SchedStrategy.INF).on(P3);
        T3.setThinkTime(Exp.fitMean(10.0));

        Entry E1 = new Entry(model, "Entry").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);
        Entry E3 = new Entry(model, "E1").on(T3);

        Activity A1 = new Activity(model, "A1", Exp.fitMean(1.0)).on(T1);
        A1.boundTo(E1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(2.0)).on(T1);
        Activity A3 = new Activity(model, "A3", Exp.fitMean(3.0)).on(T1);
        A3.synchCall(E2, 1);

        Activity A4 = new Activity(model, "B1", Exp.fitMean(0.1)).on(T2);
        A4.boundTo(E2);
        Activity A5 = new Activity(model, "B2", Exp.fitMean(0.2)).on(T2);
        Activity A6 = new Activity(model, "B3", Exp.fitMean(0.3)).on(T2);
        Activity A7 = new Activity(model, "B4", Exp.fitMean(0.4)).on(T2);
        Activity A8 = new Activity(model, "B5", Exp.fitMean(0.5)).on(T2);
        Activity A9 = new Activity(model, "B6", Exp.fitMean(0.6)).on(T2);
        A9.synchCall(E3, 1);
        A9.repliesTo(E2);
        Activity A10 = new Activity(model, "C1", Exp.fitMean(0.1)).on(T3);
        A10.boundTo(E3);
        Activity A11 = new Activity(model, "C2", Exp.fitMean(0.2)).on(T3);
        Activity A12 = new Activity(model, "C3", Exp.fitMean(0.3)).on(T3);
        Activity A13 = new Activity(model, "C4", Exp.fitMean(0.4)).on(T3);
        Activity A14 = new Activity(model, "C5", Exp.fitMean(0.5)).on(T3);
        A14.repliesTo(E3);

        // Sequential Activity Precedence
        T2.addPrecedence(ActivityPrecedence.Serial("B4", "B5"));

        // Loop Activity Precedence
        ArrayList<String> precActs = new ArrayList<String>();
        precActs.add("A2");
        precActs.add("A3");
        T1.addPrecedence(ActivityPrecedence.Loop("A1", precActs, Matrix.singleton(3)));

        // OrFork Activity Precedence
        precActs = new ArrayList<String>();
        Matrix probs = new Matrix(1, 3);
        precActs.add("C2");
        precActs.add("C3");
        precActs.add("C4");
        probs.set(0, 0, 0.3);
        probs.set(0, 1, 0.3);
        probs.set(0, 2, 0.4);
        T3.addPrecedence(ActivityPrecedence.OrFork("C1", precActs, probs));

        // OrJoin Activity Precedence
        precActs = new ArrayList<String>();
        precActs.add("C2");
        precActs.add("C3");
        precActs.add("C4");
        T3.addPrecedence(ActivityPrecedence.OrJoin(precActs, "C5"));

        // AndFork Activity Precedence
        //precActs = new ArrayList<String>();
        ArrayList<String> postActs = new ArrayList<String>();
        postActs.add("B2");
        postActs.add("B3");
        postActs.add("B4");
        T2.addPrecedence(ActivityPrecedence.AndFork("B1", postActs));

        // AndJoin Activity Precedence
        precActs = new ArrayList<String>();
        precActs.add("B2");
        precActs.add("B3");
        precActs.add("B5");
        T2.addPrecedence(ActivityPrecedence.AndJoin(precActs, "B6"));

        return model;
    }


    /**
     * Large enterprise layered network model (ofbizExample substitute).
     * <p>
     * Features:
     * - 5 processors representing different system layers
     * - Client layer with reference task (10 jobs) 
     * - 3 service layers with FCFS tasks
     * - Database layer with FCFS task
     * - Serial activity precedence with synchronous calls
     * - Representative of Apache OFBiz-style enterprise architecture
     * - Each service layer makes database calls
     *
     * @return configured large layered network model
     * @throws Exception if model creation fails
     */
    public static LayeredNetwork lqn_ofbiz() throws Exception {
        // Load from XML resource
        try {
            // For tests, try to load from test classpath first
            ClassLoader testClassLoader = Thread.currentThread().getContextClassLoader();
            URL resource = testClassLoader != null ? testClassLoader.getResource("ofbizExample.xml") : null;
            
            // Fall back to main class loader if test classpath doesn't work
            if (resource == null) {
                resource = LayeredModel.class.getClassLoader().getResource("ofbizExample.xml");
            }
            
            if (resource != null) {
                String xmlPath = resource.getPath();
                return LayeredNetwork.parseXML(xmlPath);
            } else {
                return null;
            }
        } catch (Exception e) {
            return null;
        }
    }

    /**
     * Main method for testing and demonstrating layered model examples.
     *
     * <p>Currently contains commented code for various testing scenarios:
     * - XML parsing and export functionality
     * - Solver integration and performance analysis
     * - Model structure inspection and visualization
     * - Ensemble analysis for layered networks
     *
     * @param args command line arguments (not used)
     * @throws Exception if any example execution fails
     */
    public static void main(String[] args) throws Exception {
        LayeredNetwork model = lqn_basic();
        //new LQNS(model).getAvgTable().print();
        new LN(model, SolverType.MVA).getAvgTable().print();
    }
}
