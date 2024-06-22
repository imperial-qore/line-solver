package jline.examples;

import jline.lang.constant.*;
import jline.lang.distributions.*;
import jline.lang.layered.*;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.ln.SolverLN;
import jline.util.Matrix;

import java.util.Arrays;

/**
 * Examples of layered queueing network models
 */
public class LayeredModel {

    public static LayeredNetwork testSimple() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_simple");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 50, SchedStrategy.REF);
        T1.setThinkTime(new Exp(0.01));
        T1.on(P1);

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);

        Activity A1 = new Activity(model, "A1", new Exp(1));
        A1.on(T1);
        A1.boundTo(E1);
//        T1.addPrecedence(ActivityPrecedence.Sequence("A1","A2"));
        return model;
    }

    public static LayeredNetwork test0() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_LQN_1");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 50, SchedStrategy.REF);
        Task T2 = new Task(model, "T2", 50, SchedStrategy.FCFS);
        T1.setThinkTime(new Exp(1));
        T1.on(P1);
        T2.setThinkTime(new Exp(1));
        T2.on(P2);

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);

        Entry E2 = new Entry(model, "E2");
        E2.on(T2);


        Activity A1 = new Activity(model, "A1", new Exp(1));
        A1.on(T1);
        A1.boundTo(E1);
        A1.synchCall(E2, 1);

        Activity A2 = new Activity(model, "A2", new Exp(1));
        A2.on(T2);
        A2.boundTo(E2);
        A2.repliesTo(E2);
//        T1.addPrecedence(ActivityPrecedence.Sequence("A1","A2"));
        return model;
    }

    public static LayeredNetwork test1() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_LQN_1");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 50, SchedStrategy.REF);
        T1.setThinkTime(new Exp(1.0 / 100));
        T1.on(P1);

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);

        Activity A1 = new Activity(model, "A1", new Exp(10));
        A1.on(T1);
        A1.boundTo(E1);

        Activity A2 = new Activity(model, "A2", new Exp(1.0 / 1.5));
        A2.on(T1);

        T1.addPrecedence(ActivityPrecedence.Sequence("A1", "A2"));

        return model;
    }

    public static LayeredNetwork test2() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_LQN_2");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 10, SchedStrategy.REF);
        T1.on(P1);
        T1.setThinkTime(new Exp(1.0 / 100));

        Task T2 = new Task(model, "T2", 1, SchedStrategy.FCFS);
        T2.on(P2);
        T2.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);

        Entry E2 = new Entry(model, "E2");
        E2.on(T2);

        Activity A1 = new Activity(model, "A1", new Exp(1 / 1.6));
        A1.on(T1);
        A1.boundTo(E1);

        Activity A2 = new Activity(model, "A2", new Immediate());
        A2.on(T1);

        Activity A3 = new Activity(model, "A3", new Immediate());
        A3.on(T1);
        A3.synchCall(E2, 1);

        Activity A4 = new Activity(model, "A4", new Exp(1.0 / 5));
        A4.on(T2);
        A4.boundTo(E2);

        Activity A5 = new Activity(model, "A5", new Exp(1));
        A5.on(T2);
        A5.repliesTo(E2);

        T1.addPrecedence(ActivityPrecedence.Sequence("A1", "A2"));
        T1.addPrecedence(ActivityPrecedence.Sequence("A2", "A3"));
        T2.addPrecedence(ActivityPrecedence.Sequence("A4", "A5"));

        return model;
    }
    public static LayeredNetwork test35() throws Exception {

        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", Integer.MAX_VALUE, SchedStrategy.INF);
        Processor P2 = new Processor(model, "P2", Integer.MAX_VALUE, SchedStrategy.INF);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF); T1.on(P1);
        T1.setThinkTime(Erlang.fitMeanAndSCV(0.0001,0.5));
        Task T2 = new Task(model, "T2", Integer.MAX_VALUE, SchedStrategy.INF); T2.on(P2);
        T2.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1"); E1.on(T1);
        Entry E2 = new Entry(model, "E2"); E2.on(T2);

        Activity A1 = new Activity(model, "A1", new Exp(1)); A1.on(T1); A1.boundTo(E1); A1.synchCall(E2,3);
        Activity A2 = new Activity(model, "A2", APH.fitMeanAndSCV(1,10)); A2.on(T2); A2.boundTo(E2); A2.repliesTo(E2);

        return model;
    }
    public static LayeredNetwork test3() throws Exception {

        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", Integer.MAX_VALUE, SchedStrategy.INF);
        Processor P2 = new Processor(model, "P2", Integer.MAX_VALUE, SchedStrategy.INF);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF); T1.on(P1);
        T1.setThinkTime(Erlang.fitMeanAndSCV(0.0001,0.5));
        Task T2 = new Task(model, "T2", Integer.MAX_VALUE, SchedStrategy.INF); T2.on(P2);
        T2.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1"); E1.on(T1);
        Entry E2 = new Entry(model, "E2"); E2.on(T2);

        Activity A1 = new Activity(model, "A1", new Exp(1)); A1.on(T1); A1.boundTo(E1); A1.synchCall(E2,3);
        Activity A2 = new Activity(model, "A2", APH.fitMeanAndSCV(1.0,10.0)); A2.on(T2); A2.boundTo(E2); A2.repliesTo(E2);

        return model;
    }

    public static LayeredNetwork test4() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_LQN_4");
        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 10, SchedStrategy.REF);
        T1.on(P1);
        T1.setThinkTime(new Exp(100));

        Task T2 = new Task(model, "T2", 1, SchedStrategy.FCFS);
        T2.on(P2);
        T2.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);

        Entry E2 = new Entry(model, "E2");
        E2.on(T2);

        Activity A1 = new Activity(model, "A1", new Exp( 1.6));
        A1.on(T1);
        A1.boundTo(E1);

        Activity A2 = new Activity(model, "A2", new Immediate());
        A2.on(T1);
        A2.synchCall(E2);

        Activity A3 = new Activity(model, "A3", new Exp(5));
        A3.on(T2);
        A3.boundTo(E2);

        Activity A4 = new Activity(model, "A4", new Exp(1));
        A4.on(T2);
        A4.repliesTo(E2);

        T1.addPrecedence(ActivityPrecedence.Sequence("A1", "A2"));
        T2.addPrecedence(ActivityPrecedence.Sequence("A3", "A4"));

        return model;
    }

    public static LayeredNetwork testAndForkJoin() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_and_fork_join");
        Processor P1 = new Processor(model, "P1", 10, SchedStrategy.INF);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF);
        T1.on(P1);
        T1.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);

        Activity A1 = new Activity(model, "A1", new Exp((double) 1/2));
        A1.on(T1);
        A1.boundTo(E1);

        Activity A2 = new Activity(model, "A2", new Exp((double)1/3));
        A2.on(T1);

        Activity A3 = new Activity(model, "A3", new Exp((double)1/4));
        A3.on(T1);

        Activity A4 = new Activity(model, "A4", new Exp((double)1/5));
        A4.on(T1);

        Activity A5 = new Activity(model, "A5", new Exp((double)1/6));
        A5.on(T1);


        T1.addPrecedence(ActivityPrecedence.AndFork("A1", Arrays.asList("A2", "A3", "A4"), new Matrix(0, 0)));
        T1.addPrecedence(ActivityPrecedence.AndJoin(Arrays.asList("A2", "A3", "A4"), "A5", new Matrix(1)));

        return model;
    }

    public static LayeredNetwork testOrForkJoin() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_or_fork_join");
        Processor P1 = new Processor(model, "P1", 10, SchedStrategy.FCFS);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF);
        T1.on(P1);
        T1.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);

        Activity A1 = new Activity(model, "A1", new Exp((double) 1/2));
        A1.on(T1);
        A1.boundTo(E1);

        Activity A2 = new Activity(model, "A2", new Exp((double)1/3));
        A2.on(T1);

        Activity A3 = new Activity(model, "A3", new Exp((double)1/4));
        A3.on(T1);

        Activity A4 = new Activity(model, "A4", new Exp((double)1/5));
        A4.on(T1);

        Activity A5 = new Activity(model, "A5", new Exp((double)1/6));
        A5.on(T1);

        Matrix params = new Matrix(1, 3);
        params.set(0,0,0.3);
        params.set(0,1,0.3);
        params.set(0,2,0.4);
        T1.addPrecedence(ActivityPrecedence.OrFork("A1", Arrays.asList("A2", "A3", "A4"),  params));
        T1.addPrecedence(ActivityPrecedence.OrJoin(Arrays.asList("A2", "A3", "A4"), "A5"));

        return model;
    }


    public static LayeredNetwork testLoop() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_loop");
        Processor P1 = new Processor(model, "P1", 10, SchedStrategy.INF);
//        Processor P2 = new Processor(model, "P2", 10, SchedStrategy.INF);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF);
        T1.on(P1);
        T1.setThinkTime(new Immediate());
//
//        Task T2 = new Task(model, "T2", 1, SchedStrategy.INF);
//        T2.on(P2);
//        T2.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);

//        Entry E2 = new Entry(model, "E2");
//        E2.on(T2);

        Activity A1 = new Activity(model, "A1", new Exp( (double)1/2));
        A1.on(T1);
        A1.boundTo(E1);

        Activity A2 = new Activity(model, "A2", new Exp((double) 1/3));
        A2.on(T1);

        Activity A3 = new Activity(model, "A3", new Exp((double) 1/4));
        A3.on(T1);
//        A3.synchCall(E2);

//        Activity B1 = new Activity(model, "B1", new Exp( 1/2));
//        B1.on(T2);
//        B1.boundTo(E2);

//        Activity B2 = new Activity(model, "B2", new Exp(1/3));
//        B2.on(T2);
//        B2.repliesTo(E2);

        T1.addPrecedence(ActivityPrecedence.Loop("A1", Arrays.asList("A2", "A3"), new Matrix(2)));
//        T1.addPrecedence(ActivityPrecedence.Sequence("B1", "B2"));

        return model;
    }

    public static LayeredNetwork testAllPrecedences() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_loop_network");
        Processor P1 = new Processor(model, "P1", 10, SchedStrategy.INF);
        Processor P2 = new Processor(model, "P2", 10, SchedStrategy.INF);
        Processor P3 = new Processor(model, "P3", 5, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF);
        T1.on(P1);
        T1.setThinkTime(new Immediate());

        Task T2 = new Task(model, "T2", 1, SchedStrategy.INF);
        T2.on(P2);
        T2.setThinkTime(new Immediate());

        Task T3 = new Task(model, "T3", 20, SchedStrategy.INF);
        T3.on(P3);
        T3.setThinkTime(new Exp((double) 1 / 10));

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);

        Entry E2 = new Entry(model, "E2");
        E2.on(T2);

        Entry E3 = new Entry(model, "E3");
        E3.on(T3);

        Activity A1 = new Activity(model, "A1", new Exp( 1));
        A1.on(T1);
        A1.boundTo(E1);

        Activity A2 = new Activity(model, "A2", new Exp(1/2));
        A2.on(T1);

        Activity A3 = new Activity(model, "A3", new Exp(1/3));
        A3.on(T1);
        A3.synchCall(E2);

        Activity B1 = new Activity(model, "B1", new Exp( 10));
        B1.on(T2);
        B1.boundTo(E2);

        Activity B2 = new Activity(model, "B2", new Exp(5));
        B2.on(T2);

        Activity B3 = new Activity(model, "B3", new Exp(10/3));
        B3.on(T2);

        Activity B4 = new Activity(model, "B4", new Exp(10/4));
        B4.on(T2);

        Activity B5 = new Activity(model, "B5", new Exp(2));
        B5.on(T2);

        Activity B6 = new Activity(model, "B6", new Exp(10/6));
        B6.on(T2);
        B6.synchCall(E3);
        B6.repliesTo(E2);


        Activity C1 = new Activity(model, "C1", new Exp( 10));
        C1.on(T3);
        C1.boundTo(E3);

        Activity C2 = new Activity(model, "C2", new Exp(5));
        C2.on(T3);

        Activity C3 = new Activity(model, "C3", new Exp(10/3));
        C3.on(T3);

        Activity C4 = new Activity(model, "C4", new Exp(10/4));
        C4.on(T3);

        Activity C5 = new Activity(model, "C5", new Exp(2));
        C5.on(T3);
        C5.repliesTo(E3);

        T1.addPrecedence(ActivityPrecedence.Loop("A1", Arrays.asList("A2", "A3"), new Matrix(3)));
        T2.addPrecedence(ActivityPrecedence.Sequence("B4", "B5"));
        T2.addPrecedence(ActivityPrecedence.AndFork("B1", Arrays.asList("B2", "B3", "B4"), new Matrix(0, 0)));
        T2.addPrecedence(ActivityPrecedence.AndJoin(Arrays.asList("B2", "B3", "B5"), "B6", new Matrix(0, 0)));
        T3.addPrecedence(ActivityPrecedence.OrFork("C1", Arrays.asList("C2", "C3", "C4"), new Matrix(Arrays.asList(0.3, 0.3, 0.4))));
        T3.addPrecedence(ActivityPrecedence.OrJoin(Arrays.asList("C2", "C3", "C4"), "C5"));

        return model;
    }

    public static LayeredNetwork ex4() throws Exception {

        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "FrontEnd_CPU_Processor", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "USAGE_DELAY_Processor", 1000, SchedStrategy.INF);
        Processor P3 = new Processor(model, "UsageScenario_userType1_1_Processor", 1, SchedStrategy.FCFS);
        Processor P4 = new Processor(model, "RequestHandler_HandlerIF_main_345_Processor", 1, SchedStrategy.FCFS);
        Processor P5 = new Processor(model, "RequestHandler_HandlerIF_login_345_Processor", 1, SchedStrategy.FCFS);
        Processor P6 = new Processor(model, "RequestHandler_HandlerIF_checkLogin_345_Processor", 1, SchedStrategy.FCFS);
        Processor P7 = new Processor(model, "RequestHandler_HandlerIF_logout_345_Processor", 1, SchedStrategy.FCFS);
        Processor P8 = new Processor(model, "UsageScenario_userType2_7_Processor", 1, SchedStrategy.FCFS);
        Processor P9 = new Processor(model, "RequestHandler_HandlerIF_quickadd_345_Processor", 1, SchedStrategy.FCFS);

        Task T1 = new Task(model, "FrontEnd_CPU_Task", 1000, SchedStrategy.INF); T1.on(P1);
        T1.setThinkTime(new Immediate());
        Task T2 = new Task(model, "USAGE_DELAY_Task", 1000, SchedStrategy.INF); T2.on(P2);
        T2.setThinkTime(new Immediate());
        Task T3 = new Task(model, "UsageScenario_userType1_1_Task", 10, SchedStrategy.REF); T3.on(P3);
        T3.setThinkTime(new Exp(0.1));
        Task T4 = new Task(model, "RequestHandler_HandlerIF_main_345_Task", 1000, SchedStrategy.INF); T4.on(P4);
        T4.setThinkTime(new Immediate());
        Task T5 = new Task(model, "RequestHandler_HandlerIF_login_345_Task", 1000, SchedStrategy.INF); T5.on(P5);
        T5.setThinkTime(new Immediate());
        Task T6 = new Task(model, "RequestHandler_HandlerIF_checkLogin_345_Task", 1000, SchedStrategy.INF); T6.on(P6);
        T6.setThinkTime(new Immediate());
        Task T7 = new Task(model, "RequestHandler_HandlerIF_logout_345_Task", 1000, SchedStrategy.INF); T7.on(P7);
        T7.setThinkTime(new Immediate());
        Task T8 = new Task(model, "UsageScenario_userType2_7_Task", 10, SchedStrategy.REF); T8.on(P8);
        T8.setThinkTime(new Exp(0.0714286));
        Task T9 = new Task(model, "RequestHandler_HandlerIF_quickadd_345_Task", 1000, SchedStrategy.INF); T9.on(P9);
        T9.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "FrontEnd_CPU_Entry"); E1.on(T1);
        Entry E2 = new Entry(model, "InternalAction_main__Iu1-wMhoEeKON4DtRoKCMw_34_50_Entry"); E2.on(T1);
        Entry E3 = new Entry(model, "InternalAction_login__YgxRGw26EeSPwb7XgvxhWQ_34_50_Entry"); E3.on(T1);
        Entry E4 = new Entry(model, "InternalAction_CheckLogin__vHMZsMhoEeKON4DtRoKCMw_34_50_Entry"); E4.on(T1);
        Entry E5 = new Entry(model, "InternalAction_logout__j-2E4MhrEeKON4DtRoKCMw_34_50_Entry"); E5.on(T1);
        Entry E6 = new Entry(model, "InternalAction_addcart__5JEHQMhoEeKON4DtRoKCMw_34_50_Entry"); E6.on(T1);
        Entry E7 = new Entry(model, "USAGE_DELAY0_Entry"); E7.on(T2);
        Entry E8 = new Entry(model, "UsageScenario_userType1_1_Entry"); E8.on(T3);
        Entry E9 = new Entry(model, "RequestHandler_HandlerIF_main_345_Entry"); E9.on(T4);
        Entry E10 = new Entry(model, "RequestHandler_HandlerIF_login_345_Entry"); E10.on(T5);
        Entry E11 = new Entry(model, "RequestHandler_HandlerIF_checkLogin_345_Entry"); E11.on(T6);
        Entry E12 = new Entry(model, "RequestHandler_HandlerIF_logout_345_Entry"); E12.on(T7);
        Entry E13 = new Entry(model, "UsageScenario_userType2_7_Entry"); E13.on(T8);
        Entry E14 = new Entry(model, "RequestHandler_HandlerIF_quickadd_345_Entry"); E14.on(T9);

        Activity A1 = new Activity(model, "FrontEnd_CPU_Activity", new Immediate()); A1.on(T1); A1.boundTo(E1); A1.repliesTo(E1);
        Activity A2 = new Activity(model, "InternalAction_main__Iu1-wMhoEeKON4DtRoKCMw_34_50_Activity", new Exp(100)); A2.on(T1); A2.boundTo(E2); A2.repliesTo(E2);
        Activity A3 = new Activity(model, "InternalAction_login__YgxRGw26EeSPwb7XgvxhWQ_34_50_Activity", new Exp(100)); A3.on(T1); A3.boundTo(E3); A3.repliesTo(E3);
        Activity A4 = new Activity(model, "InternalAction_CheckLogin__vHMZsMhoEeKON4DtRoKCMw_34_50_Activity", new Exp(100)); A4.on(T1); A4.boundTo(E4); A4.repliesTo(E4);
        Activity A5 = new Activity(model, "InternalAction_logout__j-2E4MhrEeKON4DtRoKCMw_34_50_Activity", new Exp(100)); A5.on(T1); A5.boundTo(E5); A5.repliesTo(E5);
        Activity A6 = new Activity(model, "InternalAction_addcart__5JEHQMhoEeKON4DtRoKCMw_34_50_Activity", new Exp(100)); A6.on(T1); A6.boundTo(E6); A6.repliesTo(E6);
        Activity A7 = new Activity(model, "USAGE_DELAY0_Activity", new Immediate()); A7.on(T2); A7.boundTo(E7); A7.repliesTo(E7);
        Activity A8 = new Activity(model, "Start2", new Immediate()); A8.on(T3); A8.boundTo(E8);
        Activity A9 = new Activity(model, "EntryLevelSystemCallcall1main1", new Immediate()); A9.on(T3); A9.synchCall(E9,1);
        Activity A10 = new Activity(model, "EntryLevelSystemCallcall1login", new Immediate()); A10.on(T3); A10.synchCall(E10,1);
        Activity A11 = new Activity(model, "EntryLevelSystemCallcall1checkLogin1", new Immediate()); A11.on(T3); A11.synchCall(E11,1);
        Activity A12 = new Activity(model, "EntryLevelSystemCallcall1checkLogin2", new Immediate()); A12.on(T3); A12.synchCall(E11,1);
        Activity A13 = new Activity(model, "EntryLevelSystemCallcall1logout", new Immediate()); A13.on(T3); A13.synchCall(E12,1);
        Activity A14 = new Activity(model, "EntryLevelSystemCallcall1main2", new Immediate()); A14.on(T3); A14.synchCall(E9,1);
        Activity A15 = new Activity(model, "Stop6", new Immediate()); A15.on(T3);
        Activity A16 = new Activity(model, "StartAction_start__EkLVIMhoEeKON4DtRoKCMw_34_5", new Immediate()); A16.on(T4); A16.boundTo(E9);
        Activity A17 = new Activity(model, "InternalAction_main__Iu1-wMhoEeKON4DtRoKCMw_34_50", new Immediate()); A17.on(T4); A17.synchCall(E2,1);
        Activity A18 = new Activity(model, "StopAction_stop__EkL8MMhoEeKON4DtRoKCMw_34_5", new Immediate()); A18.on(T4); A18.repliesTo(E9);
        Activity A19 = new Activity(model, "StartAction_aName__XJKk8g26EeSPwb7XgvxhWQ_34_5", new Immediate()); A19.on(T5); A19.boundTo(E10);
        Activity A20 = new Activity(model, "InternalAction_login__YgxRGw26EeSPwb7XgvxhWQ_34_50", new Immediate()); A20.on(T5); A20.synchCall(E3,1);
        Activity A21 = new Activity(model, "StopAction_aName__ae0SIA26EeSPwb7XgvxhWQ_34_5", new Immediate()); A21.on(T5); A21.repliesTo(E10);
        Activity A22 = new Activity(model, "StartAction_start__EkI44MhoEeKON4DtRoKCMw_34_5", new Immediate()); A22.on(T6); A22.boundTo(E11);
        Activity A23 = new Activity(model, "InternalAction_CheckLogin__vHMZsMhoEeKON4DtRoKCMw_34_50", new Immediate()); A23.on(T6); A23.synchCall(E4,1);
        Activity A24 = new Activity(model, "StopAction_stop__EkJf8MhoEeKON4DtRoKCMw_34_5", new Immediate()); A24.on(T6); A24.repliesTo(E11);
        Activity A25 = new Activity(model, "StartAction_start__EkVtMMhoEeKON4DtRoKCMw_34_5", new Immediate()); A25.on(T7); A25.boundTo(E12);
        Activity A26 = new Activity(model, "InternalAction_logout__j-2E4MhrEeKON4DtRoKCMw_34_50", new Immediate()); A26.on(T7); A26.synchCall(E5,1);
        Activity A27 = new Activity(model, "StopAction_stop__EkVtMchoEeKON4DtRoKCMw_34_5", new Immediate()); A27.on(T7); A27.repliesTo(E12);
        Activity A28 = new Activity(model, "Start8", new Immediate()); A28.on(T8); A28.boundTo(E13);
        Activity A29 = new Activity(model, "EntryLevelSystemCallcall2main1", new Immediate()); A29.on(T8); A29.synchCall(E9,1);
        Activity A30 = new Activity(model, "EntryLevelSystemCallcall2login", new Immediate()); A30.on(T8); A30.synchCall(E10,1);
        Activity A31 = new Activity(model, "EntryLevelSystemCallcall2checkLogin1", new Immediate()); A31.on(T8); A31.synchCall(E11,1);
        Activity A32 = new Activity(model, "EntryLevelSystemCallcall2checkLogin2", new Immediate()); A32.on(T8); A32.synchCall(E11,1);
        Activity A33 = new Activity(model, "EntryLevelSystemCallcall2main2", new Immediate()); A33.on(T8); A33.synchCall(E9,1);
        Activity A34 = new Activity(model, "EntryLevelSystemCallcall2quickadd", new Immediate()); A34.on(T8); A34.synchCall(E14,1);
        Activity A35 = new Activity(model, "EntryLevelSystemCallcall2logout", new Immediate()); A35.on(T8); A35.synchCall(E12,1);
        Activity A36 = new Activity(model, "EntryLevelSystemCallcall2main3", new Immediate()); A36.on(T8); A36.synchCall(E9,1);
        Activity A37 = new Activity(model, "Stop9", new Immediate()); A37.on(T8);
        Activity A38 = new Activity(model, "StartAction_start__EkBkIMhoEeKON4DtRoKCMw_34_5", new Immediate()); A38.on(T9); A38.boundTo(E14);
        Activity A39 = new Activity(model, "InternalAction_addcart__5JEHQMhoEeKON4DtRoKCMw_34_50", new Immediate()); A39.on(T9); A39.synchCall(E6,1);
        Activity A40 = new Activity(model, "StopAction_stop__EkCLMMhoEeKON4DtRoKCMw_34_5", new Immediate()); A40.on(T9); A40.repliesTo(E14);

        T3.addPrecedence(ActivityPrecedence.Sequence("Start2", "EntryLevelSystemCallcall1main1"));
        T3.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall1main1", "EntryLevelSystemCallcall1login"));
        T3.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall1login", "EntryLevelSystemCallcall1checkLogin1"));
        T3.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall1checkLogin1", "EntryLevelSystemCallcall1checkLogin2"));
        T3.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall1checkLogin2", "EntryLevelSystemCallcall1logout"));
        T3.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall1logout", "EntryLevelSystemCallcall1main2"));
        T3.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall1main2", "Stop6"));
        T4.addPrecedence(ActivityPrecedence.Sequence("StartAction_start__EkLVIMhoEeKON4DtRoKCMw_34_5", "InternalAction_main__Iu1-wMhoEeKON4DtRoKCMw_34_50"));
        T4.addPrecedence(ActivityPrecedence.Sequence("InternalAction_main__Iu1-wMhoEeKON4DtRoKCMw_34_50", "StopAction_stop__EkL8MMhoEeKON4DtRoKCMw_34_5"));
        T5.addPrecedence(ActivityPrecedence.Sequence("StartAction_aName__XJKk8g26EeSPwb7XgvxhWQ_34_5", "InternalAction_login__YgxRGw26EeSPwb7XgvxhWQ_34_50"));
        T5.addPrecedence(ActivityPrecedence.Sequence("InternalAction_login__YgxRGw26EeSPwb7XgvxhWQ_34_50", "StopAction_aName__ae0SIA26EeSPwb7XgvxhWQ_34_5"));
        T6.addPrecedence(ActivityPrecedence.Sequence("StartAction_start__EkI44MhoEeKON4DtRoKCMw_34_5", "InternalAction_CheckLogin__vHMZsMhoEeKON4DtRoKCMw_34_50"));
        T6.addPrecedence(ActivityPrecedence.Sequence("InternalAction_CheckLogin__vHMZsMhoEeKON4DtRoKCMw_34_50", "StopAction_stop__EkJf8MhoEeKON4DtRoKCMw_34_5"));
        T7.addPrecedence(ActivityPrecedence.Sequence("StartAction_start__EkVtMMhoEeKON4DtRoKCMw_34_5", "InternalAction_logout__j-2E4MhrEeKON4DtRoKCMw_34_50"));
        T7.addPrecedence(ActivityPrecedence.Sequence("InternalAction_logout__j-2E4MhrEeKON4DtRoKCMw_34_50", "StopAction_stop__EkVtMchoEeKON4DtRoKCMw_34_5"));
        T8.addPrecedence(ActivityPrecedence.Sequence("Start8", "EntryLevelSystemCallcall2main1"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2main1", "EntryLevelSystemCallcall2login"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2login", "EntryLevelSystemCallcall2checkLogin1"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2checkLogin1", "EntryLevelSystemCallcall2checkLogin2"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2checkLogin2", "EntryLevelSystemCallcall2main2"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2main2", "EntryLevelSystemCallcall2quickadd"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2quickadd", "EntryLevelSystemCallcall2logout"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2logout", "EntryLevelSystemCallcall2main3"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2main3", "Stop9"));
        T9.addPrecedence(ActivityPrecedence.Sequence("StartAction_start__EkBkIMhoEeKON4DtRoKCMw_34_5", "InternalAction_addcart__5JEHQMhoEeKON4DtRoKCMw_34_50"));
        T9.addPrecedence(ActivityPrecedence.Sequence("InternalAction_addcart__5JEHQMhoEeKON4DtRoKCMw_34_50", "StopAction_stop__EkCLMMhoEeKON4DtRoKCMw_34_5"));

        return model;
    }

    public static LayeredNetwork  testOfBizFCFS() throws Exception {

        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "FrontEnd_CPU_Processor", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "USAGE_DELAY_Processor", 1, SchedStrategy.PS);
        Processor P3 = new Processor(model, "UsageScenario_userType1_1_Processor", 1, SchedStrategy.PS);
        Processor P4 = new Processor(model, "RequestHandler_HandlerIF_main_345_Processor", 1, SchedStrategy.PS);
        Processor P5 = new Processor(model, "RequestHandler_HandlerIF_login_345_Processor", 1, SchedStrategy.PS);
        Processor P6 = new Processor(model, "RequestHandler_HandlerIF_checkLogin_345_Processor", 1, SchedStrategy.PS);
        Processor P7 = new Processor(model, "RequestHandler_HandlerIF_logout_345_Processor", 1, SchedStrategy.PS);
        Processor P8 = new Processor(model, "UsageScenario_userType2_7_Processor", 1, SchedStrategy.PS);
        Processor P9 = new Processor(model, "RequestHandler_HandlerIF_quickadd_345_Processor", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "FrontEnd_CPU_Task", 1, SchedStrategy.FCFS); T1.on(P1);
        T1.setThinkTime(new Exp(1/5.0));
        Task T2 = new Task(model, "USAGE_DELAY_Task", 1, SchedStrategy.FCFS); T2.on(P2);
        T2.setThinkTime(new Exp(1/5.0));
        Task T3 = new Task(model, "UsageScenario_userType1_1_Task", 10, SchedStrategy.REF); T3.on(P3);
        T3.setThinkTime(new Exp(0.1));
        Task T4 = new Task(model, "RequestHandler_HandlerIF_main_345_Task", 1, SchedStrategy.FCFS); T4.on(P4);
        T4.setThinkTime(new Exp(1.0/2));
        Task T5 = new Task(model, "RequestHandler_HandlerIF_login_345_Task", 1, SchedStrategy.FCFS); T5.on(P5);
        T5.setThinkTime(new Exp(1));
        Task T6 = new Task(model, "RequestHandler_HandlerIF_checkLogin_345_Task", 1, SchedStrategy.FCFS); T6.on(P6);
        T6.setThinkTime(new Exp(1.0/2));
        Task T7 = new Task(model, "RequestHandler_HandlerIF_logout_345_Task", 1, SchedStrategy.FCFS); T7.on(P7);
        T7.setThinkTime(new Exp(1.0/5));
        Task T8 = new Task(model, "UsageScenario_userType2_7_Task", 10, SchedStrategy.REF); T8.on(P8);
        T8.setThinkTime(new Exp(0.0714286));
        Task T9 = new Task(model, "RequestHandler_HandlerIF_quickadd_345_Task", 1, SchedStrategy.FCFS); T9.on(P9);
        T9.setThinkTime(new Exp(1.0/10));

        Entry E1 = new Entry(model, "FrontEnd_CPU_Entry"); E1.on(T1);
        Entry E2 = new Entry(model, "InternalAction_main__Iu1-wMhoEeKON4DtRoKCMw_34_50_Entry"); E2.on(T1);
        Entry E3 = new Entry(model, "InternalAction_login__YgxRGw26EeSPwb7XgvxhWQ_34_50_Entry"); E3.on(T1);
        Entry E4 = new Entry(model, "InternalAction_CheckLogin__vHMZsMhoEeKON4DtRoKCMw_34_50_Entry"); E4.on(T1);
        Entry E5 = new Entry(model, "InternalAction_logout__j-2E4MhrEeKON4DtRoKCMw_34_50_Entry"); E5.on(T1);
        Entry E6 = new Entry(model, "InternalAction_addcart__5JEHQMhoEeKON4DtRoKCMw_34_50_Entry"); E6.on(T1);
        Entry E7 = new Entry(model, "USAGE_DELAY0_Entry"); E7.on(T2);
        Entry E8 = new Entry(model, "UsageScenario_userType1_1_Entry"); E8.on(T3);
        Entry E9 = new Entry(model, "RequestHandler_HandlerIF_main_345_Entry"); E9.on(T4);
        Entry E10 = new Entry(model, "RequestHandler_HandlerIF_login_345_Entry"); E10.on(T5);
        Entry E11 = new Entry(model, "RequestHandler_HandlerIF_checkLogin_345_Entry"); E11.on(T6);
        Entry E12 = new Entry(model, "RequestHandler_HandlerIF_logout_345_Entry"); E12.on(T7);
        Entry E13 = new Entry(model, "UsageScenario_userType2_7_Entry"); E13.on(T8);
        Entry E14 = new Entry(model, "RequestHandler_HandlerIF_quickadd_345_Entry"); E14.on(T9);

        Activity A1 = new Activity(model, "FrontEnd_CPU_Activity", new Exp(1.0)); A1.on(T1); A1.boundTo(E1); A1.repliesTo(E1);
        Activity A2 = new Activity(model, "InternalAction_main__Iu1-wMhoEeKON4DtRoKCMw_34_50_Activity", new Exp(100)); A2.on(T1); A2.boundTo(E2); A2.repliesTo(E2);
        Activity A3 = new Activity(model, "InternalAction_login__YgxRGw26EeSPwb7XgvxhWQ_34_50_Activity", new Exp(100)); A3.on(T1); A3.boundTo(E3); A3.repliesTo(E3);
        Activity A4 = new Activity(model, "InternalAction_CheckLogin__vHMZsMhoEeKON4DtRoKCMw_34_50_Activity", new Exp(100)); A4.on(T1); A4.boundTo(E4); A4.repliesTo(E4);
        Activity A5 = new Activity(model, "InternalAction_logout__j-2E4MhrEeKON4DtRoKCMw_34_50_Activity", new Exp(100)); A5.on(T1); A5.boundTo(E5); A5.repliesTo(E5);
        Activity A6 = new Activity(model, "InternalAction_addcart__5JEHQMhoEeKON4DtRoKCMw_34_50_Activity", new Exp(100)); A6.on(T1); A6.boundTo(E6); A6.repliesTo(E6);
        Activity A7 = new Activity(model, "USAGE_DELAY0_Activity", new Exp(1.0/5)); A7.on(T2); A7.boundTo(E7); A7.repliesTo(E7);
        Activity A8 = new Activity(model, "Start2", new Exp(1.0/3)); A8.on(T3); A8.boundTo(E8);
        Activity A9 = new Activity(model, "EntryLevelSystemCallcall1main1", new Exp(1.0/2)); A9.on(T3); A9.synchCall(E9,1);
        Activity A10 = new Activity(model, "EntryLevelSystemCallcall1login", new Exp(1.0/5)); A10.on(T3); A10.synchCall(E10,1);
        Activity A11 = new Activity(model, "EntryLevelSystemCallcall1checkLogin1", new Exp(1.0/5)); A11.on(T3); A11.synchCall(E11,1);
        Activity A12 = new Activity(model, "EntryLevelSystemCallcall1checkLogin2", new Exp(1.0/5)); A12.on(T3); A12.synchCall(E11,1);
        Activity A13 = new Activity(model, "EntryLevelSystemCallcall1logout", new Exp(1.0/5)); A13.on(T3); A13.synchCall(E12,1);
        Activity A14 = new Activity(model, "EntryLevelSystemCallcall1main2", new Exp(1.0/5)); A14.on(T3); A14.synchCall(E9,1);
        Activity A15 = new Activity(model, "Stop6", new Exp(1.0/5)); A15.on(T3);
        Activity A16 = new Activity(model, "StartAction_start__EkLVIMhoEeKON4DtRoKCMw_34_5", new Exp(1.0/5)); A16.on(T4); A16.boundTo(E9);
        Activity A17 = new Activity(model, "InternalAction_main__Iu1-wMhoEeKON4DtRoKCMw_34_50", new Exp(1.0/5)); A17.on(T4); A17.synchCall(E2,1);
        Activity A18 = new Activity(model, "StopAction_stop__EkL8MMhoEeKON4DtRoKCMw_34_5", new Exp(1.0/5)); A18.on(T4); A18.repliesTo(E9);
        Activity A19 = new Activity(model, "StartAction_aName__XJKk8g26EeSPwb7XgvxhWQ_34_5", new Exp(1.0/5)); A19.on(T5); A19.boundTo(E10);
        Activity A20 = new Activity(model, "InternalAction_login__YgxRGw26EeSPwb7XgvxhWQ_34_50", new Exp(1.0/5)); A20.on(T5); A20.synchCall(E3,1);
        Activity A21 = new Activity(model, "StopAction_aName__ae0SIA26EeSPwb7XgvxhWQ_34_5", new Exp(1.0/5)); A21.on(T5); A21.repliesTo(E10);
        Activity A22 = new Activity(model, "StartAction_start__EkI44MhoEeKON4DtRoKCMw_34_5", new Exp(1.0/5)); A22.on(T6); A22.boundTo(E11);
        Activity A23 = new Activity(model, "InternalAction_CheckLogin__vHMZsMhoEeKON4DtRoKCMw_34_50", new Exp(1.0/5)); A23.on(T6); A23.synchCall(E4,1);
        Activity A24 = new Activity(model, "StopAction_stop__EkJf8MhoEeKON4DtRoKCMw_34_5", new Exp(1.0/5)); A24.on(T6); A24.repliesTo(E11);
        Activity A25 = new Activity(model, "StartAction_start__EkVtMMhoEeKON4DtRoKCMw_34_5", new Exp(1.0/5)); A25.on(T7); A25.boundTo(E12);
        Activity A26 = new Activity(model, "InternalAction_logout__j-2E4MhrEeKON4DtRoKCMw_34_50", new Exp(1.0/5)); A26.on(T7); A26.synchCall(E5,1);
        Activity A27 = new Activity(model, "StopAction_stop__EkVtMchoEeKON4DtRoKCMw_34_5", new Exp(1.0/5)); A27.on(T7); A27.repliesTo(E12);
        Activity A28 = new Activity(model, "Start8", new Exp(1.0/5)); A28.on(T8); A28.boundTo(E13);
        Activity A29 = new Activity(model, "EntryLevelSystemCallcall2main1", new Exp(1.0/5)); A29.on(T8); A29.synchCall(E9,1);
        Activity A30 = new Activity(model, "EntryLevelSystemCallcall2login", new Exp(1.0/5)); A30.on(T8); A30.synchCall(E10,1);
        Activity A31 = new Activity(model, "EntryLevelSystemCallcall2checkLogin1", new Exp(1.0/5)); A31.on(T8); A31.synchCall(E11,1);
        Activity A32 = new Activity(model, "EntryLevelSystemCallcall2checkLogin2", new Exp(1.0/5)); A32.on(T8); A32.synchCall(E11,1);
        Activity A33 = new Activity(model, "EntryLevelSystemCallcall2main2", new Exp(1.0/5)); A33.on(T8); A33.synchCall(E9,1);
        Activity A34 = new Activity(model, "EntryLevelSystemCallcall2quickadd", new Exp(1.0/5)); A34.on(T8); A34.synchCall(E14,1);
        Activity A35 = new Activity(model, "EntryLevelSystemCallcall2logout", new Exp(1.0/5)); A35.on(T8); A35.synchCall(E12,1);
        Activity A36 = new Activity(model, "EntryLevelSystemCallcall2main3", new Exp(1.0/5)); A36.on(T8); A36.synchCall(E9,1);
        Activity A37 = new Activity(model, "Stop9", new Exp(1.0/5)); A37.on(T8);
        Activity A38 = new Activity(model, "StartAction_start__EkBkIMhoEeKON4DtRoKCMw_34_5", new Exp(1.0/5)); A38.on(T9); A38.boundTo(E14);
        Activity A39 = new Activity(model, "InternalAction_addcart__5JEHQMhoEeKON4DtRoKCMw_34_50", new Exp(1.0/5)); A39.on(T9); A39.synchCall(E6,1);
        Activity A40 = new Activity(model, "StopAction_stop__EkCLMMhoEeKON4DtRoKCMw_34_5", new Exp(1.0/5)); A40.on(T9); A40.repliesTo(E14);

        T3.addPrecedence(ActivityPrecedence.Sequence("Start2", "EntryLevelSystemCallcall1main1"));
        T3.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall1main1", "EntryLevelSystemCallcall1login"));
        T3.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall1login", "EntryLevelSystemCallcall1checkLogin1"));
        T3.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall1checkLogin1", "EntryLevelSystemCallcall1checkLogin2"));
        T3.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall1checkLogin2", "EntryLevelSystemCallcall1logout"));
        T3.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall1logout", "EntryLevelSystemCallcall1main2"));
        T3.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall1main2", "Stop6"));
        T4.addPrecedence(ActivityPrecedence.Sequence("StartAction_start__EkLVIMhoEeKON4DtRoKCMw_34_5", "InternalAction_main__Iu1-wMhoEeKON4DtRoKCMw_34_50"));
        T4.addPrecedence(ActivityPrecedence.Sequence("InternalAction_main__Iu1-wMhoEeKON4DtRoKCMw_34_50", "StopAction_stop__EkL8MMhoEeKON4DtRoKCMw_34_5"));
        T5.addPrecedence(ActivityPrecedence.Sequence("StartAction_aName__XJKk8g26EeSPwb7XgvxhWQ_34_5", "InternalAction_login__YgxRGw26EeSPwb7XgvxhWQ_34_50"));
        T5.addPrecedence(ActivityPrecedence.Sequence("InternalAction_login__YgxRGw26EeSPwb7XgvxhWQ_34_50", "StopAction_aName__ae0SIA26EeSPwb7XgvxhWQ_34_5"));
        T6.addPrecedence(ActivityPrecedence.Sequence("StartAction_start__EkI44MhoEeKON4DtRoKCMw_34_5", "InternalAction_CheckLogin__vHMZsMhoEeKON4DtRoKCMw_34_50"));
        T6.addPrecedence(ActivityPrecedence.Sequence("InternalAction_CheckLogin__vHMZsMhoEeKON4DtRoKCMw_34_50", "StopAction_stop__EkJf8MhoEeKON4DtRoKCMw_34_5"));
        T7.addPrecedence(ActivityPrecedence.Sequence("StartAction_start__EkVtMMhoEeKON4DtRoKCMw_34_5", "InternalAction_logout__j-2E4MhrEeKON4DtRoKCMw_34_50"));
        T7.addPrecedence(ActivityPrecedence.Sequence("InternalAction_logout__j-2E4MhrEeKON4DtRoKCMw_34_50", "StopAction_stop__EkVtMchoEeKON4DtRoKCMw_34_5"));
        T8.addPrecedence(ActivityPrecedence.Sequence("Start8", "EntryLevelSystemCallcall2main1"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2main1", "EntryLevelSystemCallcall2login"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2login", "EntryLevelSystemCallcall2checkLogin1"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2checkLogin1", "EntryLevelSystemCallcall2checkLogin2"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2checkLogin2", "EntryLevelSystemCallcall2main2"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2main2", "EntryLevelSystemCallcall2quickadd"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2quickadd", "EntryLevelSystemCallcall2logout"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2logout", "EntryLevelSystemCallcall2main3"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2main3", "Stop9"));
        T9.addPrecedence(ActivityPrecedence.Sequence("StartAction_start__EkBkIMhoEeKON4DtRoKCMw_34_5", "InternalAction_addcart__5JEHQMhoEeKON4DtRoKCMw_34_50"));
        T9.addPrecedence(ActivityPrecedence.Sequence("InternalAction_addcart__5JEHQMhoEeKON4DtRoKCMw_34_50", "StopAction_stop__EkCLMMhoEeKON4DtRoKCMw_34_5"));

//        SolverLN solver = new SolverLN(model);
//        solver.getEnsemble().get(7).getConnectionMatrix().print();
//        for(int i=0;i<solver.getEnsemble().size();i++) {
//            System.out.println("Solver result "+i);
//            Network layer = solver.getEnsemble().get(i);
//            SolverMVA layersolver = new SolverMVA(layer);
//            layersolver.getAvgTable();
//        }
//        solver.getEnsembleAvg();
        return model;
    }

    public static void main(String[] args) throws Exception{
        LayeredNetwork model = test35();
        SolverOptions solverOptions= new SolverOptions(SolverType.LN);
        SolverLN solver = new SolverLN(model, solverOptions);

//        // Use this block to test single layer.
//        // Layer 1
//        System.out.println("The MVA solver result for Layer 1");
//        Network network1 = solver.getEnsemble().get(0);
//        SolverMVA layersolver1 = new SolverMVA(network1);
//        layersolver1.getAvgTable();
//        System.out.println("----------------------------------------------------------------------------------------");
//
//        //layer 2
//        System.out.println("The MVA solver result for Layer 2");
//        Network network2 = solver.getEnsemble().get(1);
//        SolverMVA layersolver2 = new SolverMVA(network2);
//        layersolver2.getAvgTable();
//        System.out.println("----------------------------------------------------------------------------------------");

        // Use this Line to test the solverLN result
        System.out.println("The LN solver result ");
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        System.out.println(1);
    }
}
