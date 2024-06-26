package jline.examples;

import jline.lang.constant.*;
import jline.lang.distributions.*;
import jline.lang.layered.*;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.ln.LNOptions;
import jline.solvers.ln.SolverLN;
import jline.util.Matrix;

import java.util.ArrayList;

/**
 * Examples of layered queueing network models
 */
public class LayeredModel {

    public static LayeredNetwork example_layeredModel_1() throws Exception {
        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 10, SchedStrategy.REF).on(P1).setThinkTime(new Exp(0.01));
        Task T2 = new Task(model, "T2", 1, SchedStrategy.FCFS).on(P2).setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);

        Activity A1 = new Activity(model, "AS1", new Exp(0.625)).on(T1);
        A1.boundTo(E1);
        Activity A2 = new Activity(model, "AS2", new Immediate()).on(T1);
        A2.synchCall(E2, 1);
        Activity A3 = new Activity(model, "AS3", new Exp(0.2)).on(T2);
        A3.boundTo(E2);
        Activity A4 = new Activity(model, "AS4", new Exp(1)).on(T2);
        A4.repliesTo(E2);

        T1.addPrecedence(ActivityPrecedence.Sequence("AS1", "AS2"));
        T2.addPrecedence(ActivityPrecedence.Sequence("AS3", "AS4"));

        // Model solution
        SolverLN solver = new SolverLN(model);
        solver.getEnsembleAvg();
        return model;
    }

    public static LayeredNetwork example_layeredModel_2() throws Exception {

        LayeredNetwork model = new LayeredNetwork("LQN1");

        Processor P1 = new Processor(model, "P1", Integer.MAX_VALUE, SchedStrategy.INF);
        Processor P2 = new Processor(model, "P2", Integer.MAX_VALUE, SchedStrategy.INF);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF).on(P1);
        T1.setThinkTime(Erlang.fitMeanAndSCV(0.0001, 0.5));
        Task T2 = new Task(model, "T2", Integer.MAX_VALUE, SchedStrategy.INF).on(P2);
        T2.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);

        Activity A1 = new Activity(model, "A1", new Exp(1)).on(T1);
        A1.boundTo(E1);
        A1.synchCall(E2, 3);
        Activity A2 = new Activity(model, "A2", APH.fitMeanAndSCV(1, 10)).on(T2);
        A2.boundTo(E2);
        A2.repliesTo(E2);


        // Model solution
        SolverLN solver = new SolverLN(model);
        solver.getEnsembleAvg();
        return model;
    }

    public static LayeredNetwork example_layeredModel_3() throws Exception {

        LayeredNetwork model = new LayeredNetwork("example_layeredModel_1.xml");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 10, SchedStrategy.REF).on(P1);
        T1.setThinkTime(new Exp(0.01));
        Task T2 = new Task(model, "T2", 1, SchedStrategy.FCFS).on(P2);
        T2.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);

        Activity A1 = new Activity(model, "AS1", new Exp(0.625)).on(T1);
        A1.boundTo(E1);
        Activity A2 = new Activity(model, "AS2", new Immediate()).on(T1);
        A2.synchCall(E2, 1);
        Activity A3 = new Activity(model, "AS3", new Exp(0.2)).on(T2);
        A3.boundTo(E2);
        Activity A4 = new Activity(model, "AS4", new Exp(1)).on(T2);
        A4.repliesTo(E2);

        T1.addPrecedence(ActivityPrecedence.Sequence("AS1", "AS2"));
        T2.addPrecedence(ActivityPrecedence.Sequence("AS3", "AS4"));

        // Model solution
        SolverLN solver = new SolverLN(model);
        solver.getEnsembleAvg();
        return model;
    }

    public static LayeredNetwork example_layeredModel_4() throws Exception {

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

        Task T1 = new Task(model, "FrontEnd_CPU_Task", 1000, SchedStrategy.INF);
        T1.on(P1);
        T1.setThinkTime(new Immediate());
        Task T2 = new Task(model, "USAGE_DELAY_Task", 1000, SchedStrategy.INF);
        T2.on(P2);
        T2.setThinkTime(new Immediate());
        Task T3 = new Task(model, "UsageScenario_userType1_1_Task", 10, SchedStrategy.REF);
        T3.on(P3);
        T3.setThinkTime(new Exp(0.1));
        Task T4 = new Task(model, "RequestHandler_HandlerIF_main_345_Task", 1000, SchedStrategy.INF);
        T4.on(P4);
        T4.setThinkTime(new Immediate());
        Task T5 = new Task(model, "RequestHandler_HandlerIF_login_345_Task", 1000, SchedStrategy.INF);
        T5.on(P5);
        T5.setThinkTime(new Immediate());
        Task T6 = new Task(model, "RequestHandler_HandlerIF_checkLogin_345_Task", 1000, SchedStrategy.INF);
        T6.on(P6);
        T6.setThinkTime(new Immediate());
        Task T7 = new Task(model, "RequestHandler_HandlerIF_logout_345_Task", 1000, SchedStrategy.INF);
        T7.on(P7);
        T7.setThinkTime(new Immediate());
        Task T8 = new Task(model, "UsageScenario_userType2_7_Task", 10, SchedStrategy.REF);
        T8.on(P8);
        T8.setThinkTime(new Exp(0.0714286));
        Task T9 = new Task(model, "RequestHandler_HandlerIF_quickadd_345_Task", 1000, SchedStrategy.INF);
        T9.on(P9);
        T9.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "FrontEnd_CPU_Entry");
        E1.on(T1);
        Entry E2 = new Entry(model, "InternalAction_main__Iu1-wMhoEeKON4DtRoKCMw_34_50_Entry");
        E2.on(T1);
        Entry E3 = new Entry(model, "InternalAction_login__YgxRGw26EeSPwb7XgvxhWQ_34_50_Entry");
        E3.on(T1);
        Entry E4 = new Entry(model, "InternalAction_CheckLogin__vHMZsMhoEeKON4DtRoKCMw_34_50_Entry");
        E4.on(T1);
        Entry E5 = new Entry(model, "InternalAction_logout__j-2E4MhrEeKON4DtRoKCMw_34_50_Entry");
        E5.on(T1);
        Entry E6 = new Entry(model, "InternalAction_addcart__5JEHQMhoEeKON4DtRoKCMw_34_50_Entry");
        E6.on(T1);
        Entry E7 = new Entry(model, "USAGE_DELAY0_Entry");
        E7.on(T2);
        Entry E8 = new Entry(model, "UsageScenario_userType1_1_Entry");
        E8.on(T3);
        Entry E9 = new Entry(model, "RequestHandler_HandlerIF_main_345_Entry");
        E9.on(T4);
        Entry E10 = new Entry(model, "RequestHandler_HandlerIF_login_345_Entry");
        E10.on(T5);
        Entry E11 = new Entry(model, "RequestHandler_HandlerIF_checkLogin_345_Entry");
        E11.on(T6);
        Entry E12 = new Entry(model, "RequestHandler_HandlerIF_logout_345_Entry");
        E12.on(T7);
        Entry E13 = new Entry(model, "UsageScenario_userType2_7_Entry");
        E13.on(T8);
        Entry E14 = new Entry(model, "RequestHandler_HandlerIF_quickadd_345_Entry");
        E14.on(T9);

        Activity A1 = new Activity(model, "FrontEnd_CPU_Activity", new Immediate());
        A1.on(T1);
        A1.boundTo(E1);
        A1.repliesTo(E1);
        Activity A2 = new Activity(model, "InternalAction_main__Iu1-wMhoEeKON4DtRoKCMw_34_50_Activity", new Exp(100));
        A2.on(T1);
        A2.boundTo(E2);
        A2.repliesTo(E2);
        Activity A3 = new Activity(model, "InternalAction_login__YgxRGw26EeSPwb7XgvxhWQ_34_50_Activity", new Exp(100));
        A3.on(T1);
        A3.boundTo(E3);
        A3.repliesTo(E3);
        Activity A4 = new Activity(model, "InternalAction_CheckLogin__vHMZsMhoEeKON4DtRoKCMw_34_50_Activity", new Exp(100));
        A4.on(T1);
        A4.boundTo(E4);
        A4.repliesTo(E4);
        Activity A5 = new Activity(model, "InternalAction_logout__j-2E4MhrEeKON4DtRoKCMw_34_50_Activity", new Exp(100));
        A5.on(T1);
        A5.boundTo(E5);
        A5.repliesTo(E5);
        Activity A6 = new Activity(model, "InternalAction_addcart__5JEHQMhoEeKON4DtRoKCMw_34_50_Activity", new Exp(100));
        A6.on(T1);
        A6.boundTo(E6);
        A6.repliesTo(E6);
        Activity A7 = new Activity(model, "USAGE_DELAY0_Activity", new Immediate());
        A7.on(T2);
        A7.boundTo(E7);
        A7.repliesTo(E7);
        Activity A8 = new Activity(model, "Start2", new Immediate());
        A8.on(T3);
        A8.boundTo(E8);
        Activity A9 = new Activity(model, "EntryLevelSystemCallcall1main1", new Immediate());
        A9.on(T3);
        A9.synchCall(E9, 1);
        Activity A10 = new Activity(model, "EntryLevelSystemCallcall1login", new Immediate());
        A10.on(T3);
        A10.synchCall(E10, 1);
        Activity A11 = new Activity(model, "EntryLevelSystemCallcall1checkLogin1", new Immediate());
        A11.on(T3);
        A11.synchCall(E11, 1);
        Activity A12 = new Activity(model, "EntryLevelSystemCallcall1checkLogin2", new Immediate());
        A12.on(T3);
        A12.synchCall(E11, 1);
        Activity A13 = new Activity(model, "EntryLevelSystemCallcall1logout", new Immediate());
        A13.on(T3);
        A13.synchCall(E12, 1);
        Activity A14 = new Activity(model, "EntryLevelSystemCallcall1main2", new Immediate());
        A14.on(T3);
        A14.synchCall(E9, 1);
        Activity A15 = new Activity(model, "Stop6", new Immediate());
        A15.on(T3);
        Activity A16 = new Activity(model, "StartAction_start__EkLVIMhoEeKON4DtRoKCMw_34_5", new Immediate());
        A16.on(T4);
        A16.boundTo(E9);
        Activity A17 = new Activity(model, "InternalAction_main__Iu1-wMhoEeKON4DtRoKCMw_34_50", new Immediate());
        A17.on(T4);
        A17.synchCall(E2, 1);
        Activity A18 = new Activity(model, "StopAction_stop__EkL8MMhoEeKON4DtRoKCMw_34_5", new Immediate());
        A18.on(T4);
        A18.repliesTo(E9);
        Activity A19 = new Activity(model, "StartAction_aName__XJKk8g26EeSPwb7XgvxhWQ_34_5", new Immediate());
        A19.on(T5);
        A19.boundTo(E10);
        Activity A20 = new Activity(model, "InternalAction_login__YgxRGw26EeSPwb7XgvxhWQ_34_50", new Immediate());
        A20.on(T5);
        A20.synchCall(E3, 1);
        Activity A21 = new Activity(model, "StopAction_aName__ae0SIA26EeSPwb7XgvxhWQ_34_5", new Immediate());
        A21.on(T5);
        A21.repliesTo(E10);
        Activity A22 = new Activity(model, "StartAction_start__EkI44MhoEeKON4DtRoKCMw_34_5", new Immediate());
        A22.on(T6);
        A22.boundTo(E11);
        Activity A23 = new Activity(model, "InternalAction_CheckLogin__vHMZsMhoEeKON4DtRoKCMw_34_50", new Immediate());
        A23.on(T6);
        A23.synchCall(E4, 1);
        Activity A24 = new Activity(model, "StopAction_stop__EkJf8MhoEeKON4DtRoKCMw_34_5", new Immediate());
        A24.on(T6);
        A24.repliesTo(E11);
        Activity A25 = new Activity(model, "StartAction_start__EkVtMMhoEeKON4DtRoKCMw_34_5", new Immediate());
        A25.on(T7);
        A25.boundTo(E12);
        Activity A26 = new Activity(model, "InternalAction_logout__j-2E4MhrEeKON4DtRoKCMw_34_50", new Immediate());
        A26.on(T7);
        A26.synchCall(E5, 1);
        Activity A27 = new Activity(model, "StopAction_stop__EkVtMchoEeKON4DtRoKCMw_34_5", new Immediate());
        A27.on(T7);
        A27.repliesTo(E12);
        Activity A28 = new Activity(model, "Start8", new Immediate());
        A28.on(T8);
        A28.boundTo(E13);
        Activity A29 = new Activity(model, "EntryLevelSystemCallcall2main1", new Immediate());
        A29.on(T8);
        A29.synchCall(E9, 1);
        Activity A30 = new Activity(model, "EntryLevelSystemCallcall2login", new Immediate());
        A30.on(T8);
        A30.synchCall(E10, 1);
        Activity A31 = new Activity(model, "EntryLevelSystemCallcall2checkLogin1", new Immediate());
        A31.on(T8);
        A31.synchCall(E11, 1);
        Activity A32 = new Activity(model, "EntryLevelSystemCallcall2checkLogin2", new Immediate());
        A32.on(T8);
        A32.synchCall(E11, 1);
        Activity A33 = new Activity(model, "EntryLevelSystemCallcall2main2", new Immediate());
        A33.on(T8);
        A33.synchCall(E9, 1);
        Activity A34 = new Activity(model, "EntryLevelSystemCallcall2quickadd", new Immediate());
        A34.on(T8);
        A34.synchCall(E14, 1);
        Activity A35 = new Activity(model, "EntryLevelSystemCallcall2logout", new Immediate());
        A35.on(T8);
        A35.synchCall(E12, 1);
        Activity A36 = new Activity(model, "EntryLevelSystemCallcall2main3", new Immediate());
        A36.on(T8);
        A36.synchCall(E9, 1);
        Activity A37 = new Activity(model, "Stop9", new Immediate());
        A37.on(T8);
        Activity A38 = new Activity(model, "StartAction_start__EkBkIMhoEeKON4DtRoKCMw_34_5", new Immediate());
        A38.on(T9);
        A38.boundTo(E14);
        Activity A39 = new Activity(model, "InternalAction_addcart__5JEHQMhoEeKON4DtRoKCMw_34_50", new Immediate());
        A39.on(T9);
        A39.synchCall(E6, 1);
        Activity A40 = new Activity(model, "StopAction_stop__EkCLMMhoEeKON4DtRoKCMw_34_5", new Immediate());
        A40.on(T9);
        A40.repliesTo(E14);

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

    public static LayeredNetwork example_layeredModel_5() throws Exception {

        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 100, SchedStrategy.REF).on(P1);
        T1.setThinkTime(Erlang.fitMeanAndSCV(10, 1));
        Task T2 = new Task(model, "T2", 100, SchedStrategy.INF).on(P2);
        T2.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);
        Entry E3 = new Entry(model, "E3").on(T2);

        Activity A1 = new Activity(model, "A1", new Exp(1)).on(T1);
        A1.boundTo(E1);
        A1.synchCall(E2, 1).synchCall(E3, 1);
        Activity A2 = new Activity(model, "A20", new Exp(1)).on(T2);
        A2.boundTo(E2);
        Activity A3 = new Activity(model, "A21", new Exp(1)).on(T2);
        Activity A4 = new Activity(model, "A22", new Exp(1)).on(T2);
        A4.repliesTo(E2);
        Activity A5 = new Activity(model, "A3", new Exp(1)).on(T2);
        A5.boundTo(E3);
        A5.repliesTo(E3);

        T2.addPrecedence(ActivityPrecedence.Sequence("A20", "A21"));
        T2.addPrecedence(ActivityPrecedence.Sequence("A21", "A22"));

        // Model solution
        SolverLN solver = new SolverLN(model);
        solver.getEnsembleAvg();
        return model;
    }

    public static LayeredNetwork example_layeredModel_6() throws Exception {

        LayeredNetwork model = new LayeredNetwork("example_layeredModel_6.xml");

        Processor P1 = new Processor(model, "R1_Processor", 100, SchedStrategy.FCFS);
        Processor P2 = new Processor(model, "R2_Processor", Integer.MAX_VALUE, SchedStrategy.INF);
        Processor P3 = new Processor(model, "R3_Processor", 2, SchedStrategy.FCFS);
        Processor P4 = new Processor(model, "R1A_Processor", 7, SchedStrategy.FCFS);
        Processor P5 = new Processor(model, "R1B_Processor", 3, SchedStrategy.FCFS);
        Processor P6 = new Processor(model, "R2A_Processor", 4, SchedStrategy.FCFS);
        Processor P7 = new Processor(model, "R2B_Processor", 5, SchedStrategy.FCFS);

        Task T1 = new Task(model, "R1_Task", 100, SchedStrategy.REF).on(P1);
        T1.setThinkTime(new Exp(0.05));
        Task T2 = new Task(model, "R2_Task", Integer.MAX_VALUE, SchedStrategy.INF).on(P2);
        T2.setThinkTime(new Immediate());
        Task T3 = new Task(model, "R3_Task", 2, SchedStrategy.FCFS).on(P3);
        T3.setThinkTime(new Immediate());
        Task T4 = new Task(model, "R1A_Task", 7, SchedStrategy.FCFS).on(P4);
        T4.setThinkTime(new Immediate());
        Task T5 = new Task(model, "R1B_Task", 3, SchedStrategy.FCFS).on(P5);
        T5.setThinkTime(new Immediate());
        Task T6 = new Task(model, "R2A_Task", 4, SchedStrategy.FCFS).on(P6);
        T6.setThinkTime(new Immediate());
        Task T7 = new Task(model, "R2B_Task", 5, SchedStrategy.FCFS).on(P7);
        T7.setThinkTime(new Immediate());

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

        Activity A1 = new Activity(model, "A1_Empty", new Immediate()).on(T1);
        A1.boundTo(E1);
        A1.synchCall(E5, 1);
        Activity A2 = new Activity(model, "A2_Empty", new Immediate()).on(T1);
        A2.synchCall(E6, 1);
        Activity A3 = new Activity(model, "A5_Empty", new Immediate()).on(T1);
        A3.synchCall(E9, 1);
        Activity A4 = new Activity(model, "A6_Empty", new Immediate()).on(T1);
        A4.synchCall(E10, 1);
        Activity A5 = new Activity(model, "A3_Empty", new Immediate()).on(T1);
        A5.synchCall(E7, 1);
        Activity A6 = new Activity(model, "A4_Empty", new Immediate()).on(T1);
        A6.synchCall(E8, 1);
        Activity A7 = new Activity(model, "E4_Empty", new Immediate()).on(T2);
        A7.boundTo(E2);
        Activity A8 = new Activity(model, "A7_Empty", new Immediate()).on(T2);
        A8.synchCall(E11, 1);
        Activity A9 = new Activity(model, "A8_Empty", new Immediate()).on(T2);
        A9.synchCall(E12, 1);
        Activity A10 = new Activity(model, "A9_Empty", new Immediate()).on(T2);
        A10.synchCall(E14, 1);
        Activity A11 = new Activity(model, "A11_Empty", new Immediate()).on(T2);
        A11.synchCall(E13, 1);
        A11.repliesTo(E2);
        Activity A12 = new Activity(model, "A12_Empty", new Immediate()).on(T2);
        A12.boundTo(E3);
        A12.synchCall(E16, 1);
        A12.repliesTo(E3);
        Activity A13 = new Activity(model, "A10_Empty", new Immediate()).on(T2);
        A13.synchCall(E15, 1);
        Activity A14 = new Activity(model, "A13", new Exp(0.1)).on(T3);
        A14.boundTo(E4);
        A14.repliesTo(E4);
        Activity A15 = new Activity(model, "A1", new Exp(0.142857)).on(T4);
        A15.boundTo(E5);
        A15.repliesTo(E5);
        Activity A16 = new Activity(model, "A2", new Exp(0.25)).on(T4);
        A16.boundTo(E6);
        Activity A17 = new Activity(model, "A3", new Exp(0.2)).on(T4);
        A17.boundTo(E7);
        A17.repliesTo(E7);
        Activity A18 = new Activity(model, "A2_Res_Empty", new Immediate()).on(T4);
        A18.synchCall(E2, 1);
        A18.repliesTo(E6);
        Activity A19 = new Activity(model, "A4", new Exp(0.125)).on(T5);
        A19.boundTo(E8);
        A19.repliesTo(E8);
        Activity A20 = new Activity(model, "A5", new Exp(0.25)).on(T5);
        A20.boundTo(E9);
        Activity A21 = new Activity(model, "A6", new Exp(0.166667)).on(T5);
        A21.boundTo(E10);
        A21.repliesTo(E10);
        Activity A22 = new Activity(model, "A5_Res_Empty", new Immediate()).on(T5);
        A22.synchCall(E3, 1);
        A22.repliesTo(E9);
        Activity A23 = new Activity(model, "A7", new Exp(0.166667)).on(T6);
        A23.boundTo(E11);
        A23.repliesTo(E11);
        Activity A24 = new Activity(model, "A8", new Exp(0.125)).on(T6);
        A24.boundTo(E12);
        A24.repliesTo(E12);
        Activity A25 = new Activity(model, "A11", new Exp(0.25)).on(T6);
        A25.boundTo(E13);
        A25.repliesTo(E13);
        Activity A26 = new Activity(model, "A9", new Exp(0.25)).on(T7);
        A26.boundTo(E14);
        Activity A27 = new Activity(model, "A10", new Exp(0.166667)).on(T7);
        A27.boundTo(E15);
        A27.repliesTo(E15);
        Activity A28 = new Activity(model, "A12", new Exp(0.125)).on(T7);
        A28.boundTo(E16);
        A28.repliesTo(E16);
        Activity A29 = new Activity(model, "A9_Res_Empty", new Immediate()).on(T7);
        A29.synchCall(E4, 1);
        A29.repliesTo(E14);

        T1.addPrecedence(ActivityPrecedence.Sequence("A1_Empty", "A2_Empty"));
        T1.addPrecedence(ActivityPrecedence.Sequence("A5_Empty", "A6_Empty"));
        T2.addPrecedence(ActivityPrecedence.Sequence("E4_Empty", "A7_Empty"));
        T2.addPrecedence(ActivityPrecedence.Sequence("A9_Empty", "A10_Empty"));
        T4.addPrecedence(ActivityPrecedence.Sequence("A2", "A2_Res_Empty"));
        T5.addPrecedence(ActivityPrecedence.Sequence("A5", "A5_Res_Empty"));
        T7.addPrecedence(ActivityPrecedence.Sequence("A9", "A9_Res_Empty"));

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
        //T2.addPrecedence(ActivityPrecedence.AndFork("A7_Empty", postActs));

        // OrJoin Activity Precedence
        precActs.clear();
        precActs.add("A3_Empty");
        precActs.add("A4_Empty");
        T1.addPrecedence(ActivityPrecedence.OrJoin(precActs, "A5_Empty"));

        // AndJoin Activity Precedence
        precActs.clear();
        precActs.add("A8_Empty");
        precActs.add("A10_Empty");
        //T2.addPrecedence(ActivityPrecedence.AndJoin(precActs, "A11_Empty"));

        // Model solution
        SolverLN solver = new SolverLN(model);
        solver.getEnsembleAvg();
        return model;
    }

    public static LayeredNetwork example_layeredModel_7() throws Exception {
        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", Integer.MAX_VALUE, SchedStrategy.INF);
        Processor P2 = new Processor(model, "P2", Integer.MAX_VALUE, SchedStrategy.INF);
        Processor P3 = new Processor(model, "P3", 5, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF).on(P1);
        T1.setThinkTime(new Immediate());
        Task T2 = new Task(model, "T2", Integer.MAX_VALUE, SchedStrategy.INF).on(P2);
        T2.setThinkTime(new Immediate());
        Task T3 = new Task(model, "T3", Integer.MAX_VALUE, SchedStrategy.INF).on(P3);
        T3.setThinkTime(new Exp(0.1));

        Entry E1 = new Entry(model, "Entry").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);
        Entry E3 = new Entry(model, "E1").on(T3);

        Activity A1 = new Activity(model, "A1", new Exp(1)).on(T1);
        A1.boundTo(E1);
        Activity A2 = new Activity(model, "A2", new Exp(0.5)).on(T1);
        Activity A3 = new Activity(model, "A3", new Exp(0.333333)).on(T1);
        A3.synchCall(E2, 1);
        Activity A4 = new Activity(model, "B1", new Exp(10)).on(T2);
        A4.boundTo(E2);
        Activity A5 = new Activity(model, "B2", new Exp(5)).on(T2);
        Activity A6 = new Activity(model, "B3", new Exp(3.33333)).on(T2);
        Activity A7 = new Activity(model, "B4", new Exp(2.5)).on(T2);
        Activity A8 = new Activity(model, "B5", new Exp(2)).on(T2);
        Activity A9 = new Activity(model, "B6", new Exp(1.66667)).on(T2);
        A9.synchCall(E3, 1);
        A9.repliesTo(E2);
        Activity A10 = new Activity(model, "C1", new Exp(10)).on(T3);
        A10.boundTo(E3);
        Activity A11 = new Activity(model, "C2", new Exp(5)).on(T3);
        Activity A12 = new Activity(model, "C3", new Exp(3.33333)).on(T3);
        Activity A13 = new Activity(model, "C4", new Exp(2.5)).on(T3);
        Activity A14 = new Activity(model, "C5", new Exp(2)).on(T3);
        A14.repliesTo(E3);

        T2.addPrecedence(ActivityPrecedence.Sequence("B4", "B5"));

        // Loop Activity Precedence
        ArrayList<String> precActs = new ArrayList<String>();
        precActs.add("A2");
        precActs.add("A3");
        T1.addPrecedence(ActivityPrecedence.Loop("A1", precActs, new Matrix(3)));

        // OrFork Activity Precedence
        precActs.clear();
        Matrix probs = new Matrix(1, 3);
        precActs.add("C2");
        precActs.add("C3");
        precActs.add("C4");
        probs.set(0, 0, 0.3);
        probs.set(0, 1, 0.3);
        probs.set(0, 2, 0.4);
        T3.addPrecedence(ActivityPrecedence.OrFork("C1", precActs, probs));

        // AndFork Activity Precedence
        ArrayList<String> postActs = new ArrayList<String>();
        postActs.add("B2");
        postActs.add("B3");
        postActs.add("B4");
        //T2.addPrecedence(ActivityPrecedence.AndFork("B1", postActs));

        // OrJoin Activity Precedence
        precActs.clear();
        precActs.add("C2");
        precActs.add("C3");
        precActs.add("C4");
        //T3.addPrecedence(ActivityPrecedence.OrJoin(precActs, "C5"));

        // AndJoin Activity Precedence
        precActs.clear();
        precActs.add("B2");
        precActs.add("B3");
        precActs.add("B5");
        //T2.addPrecedence(ActivityPrecedence.AndJoin(precActs, "B6"));

        // Model solution
        SolverLN solver = new SolverLN(model);
        solver.getEnsembleAvg();
        return model;
    }


    public static void main(String[] args) throws Exception {
        LayeredNetwork model = example_layeredModel_7();
        new SolverLN(model, new LNOptions()).printEnsembleAvgTables();
    }
}
