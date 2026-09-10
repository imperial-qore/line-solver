"""
Layered Queueing Network - Sock Shop Microservice Application

This example demonstrates:
- LQN model of the Sock Shop microservice benchmark (7 processors, 7 tasks)
- Processor replication (P2_1 with replication=2)
- Fan-in and fan-out for replicated task communication
- Processor sharing (PS) scheduling with quantum
- Multiple entries per task with serial activity precedence
- Synchronous calls across a multi-tier architecture

Based on: atom2021/sockshop-line-opt-ord1000-rep1.lqnx
"""

from line_solver import *
import numpy as np


def lqn_sockshop():


    model = LayeredNetwork('sockshop')

    # Processors
    P = np.empty(7, dtype=object)
    P[0] = Processor(model, 'P1', 1, SchedStrategy.INF)
    P[1] = Processor(model, 'P2_1', 1, SchedStrategy.PS)
    P[1].setQuantum(0.1)
    P[1].setReplication(2)
    P[2] = Processor(model, 'P2_2', 1, SchedStrategy.PS)
    P[2].setQuantum(0.1)
    P[3] = Processor(model, 'P2_3', 1, SchedStrategy.PS)
    P[3].setQuantum(0.1)
    P[4] = Processor(model, 'P3_1', 1, SchedStrategy.PS)
    P[4].setQuantum(0.1)
    P[5] = Processor(model, 'P3_2', 1, SchedStrategy.PS)
    P[5].setQuantum(0.1)
    P[6] = Processor(model, 'P3_3', 1, SchedStrategy.PS)
    P[6].setQuantum(0.1)

    # Tasks
    T = np.empty(7, dtype=object)
    T[0] = Task(model, 'T0', 1000, SchedStrategy.REF).on(P[0]).set_think_time(Exp.fit_mean(7.0))  # reference task
    T[1] = Task(model, 'T1', 24, SchedStrategy.FCFS).on(P[1]).set_think_time(Immediate())   # edge router
    T[1].setFanOut('T2', 1)
    T[2] = Task(model, 'T2', 21, SchedStrategy.FCFS).on(P[2]).set_think_time(Immediate())   # front end
    T[2].setFanOut('T3', 1)
    T[2].setFanOut('T4', 1)
    T[2].setFanIn('T1', 1)
    T[3] = Task(model, 'T6', 100, SchedStrategy.FCFS).on(P[3]).set_think_time(Immediate())  # cartdb
    T[3].setFanIn('T3', 1)
    T[4] = Task(model, 'T3', 139, SchedStrategy.FCFS).on(P[4]).set_think_time(Immediate())  # cart
    T[4].setFanOut('T6', 1)
    T[4].setFanIn('T2', 1)
    T[5] = Task(model, 'T4', 16, SchedStrategy.FCFS).on(P[5]).set_think_time(Immediate())   # catalog
    T[5].setFanOut('T5', 1)
    T[5].setFanIn('T2', 1)
    T[6] = Task(model, 'T5', 151, SchedStrategy.FCFS).on(P[6]).set_think_time(Immediate())  # catalogdb
    T[6].setFanIn('T4', 1)

    # Entries
    E = np.empty(12, dtype=object)
    E[0] = Entry(model, 'E0').on(T[0])
    E[1] = Entry(model, 'E1').on(T[1])
    E[2] = Entry(model, 'E2').on(T[2])
    E[3] = Entry(model, 'E3').on(T[2])
    E[4] = Entry(model, 'E4').on(T[2])
    E[5] = Entry(model, 'E11').on(T[3])
    E[6] = Entry(model, 'E5').on(T[4])
    E[7] = Entry(model, 'E6').on(T[4])
    E[8] = Entry(model, 'E7').on(T[4])
    E[9] = Entry(model, 'E8').on(T[5])
    E[10] = Entry(model, 'E9').on(T[5])
    E[11] = Entry(model, 'E10').on(T[6])

    # Activities
    A = np.empty(24, dtype=object)

    # T0: reference task (AS -> AS0 -> synchCall E1)
    A[0] = Activity(model, 'AS', Exp.fit_mean(0.00000005)).on(T[0]).bound_to(E[0])
    A[1] = Activity(model, 'AS0', Exp.fit_mean(0.00000005)).on(T[0]).synch_call(E[1], 1.0)

    # T1: edge router (AS1 -> AS2 -> synchCall E2,E3,E4)
    A[2] = Activity(model, 'AS1', Exp.fit_mean(0.00000005)).on(T[1]).bound_to(E[1])
    A[3] = Activity(model, 'AS2', Exp.fit_mean(0.0022216)).on(T[1]).synch_call(E[2], 0.33).synch_call(E[3], 0.17).synch_call(E[4], 0.50).replies_to(E[1])

    # T2: front end - Entry E2 (AH1 -> AH2)
    A[4] = Activity(model, 'AH1', Exp.fit_mean(0.00000005)).on(T[2]).bound_to(E[2])
    A[5] = Activity(model, 'AH2', Exp.fit_mean(0.0021319)).on(T[2]).replies_to(E[2])
    # T2: front end - Entry E3 (AH3 -> AH4 -> synchCall E8,E9)
    A[6] = Activity(model, 'AH3', Exp.fit_mean(0.00000005)).on(T[2]).bound_to(E[3])
    A[7] = Activity(model, 'AH4', Exp.fit_mean(0.0037561)).on(T[2]).synch_call(E[9], 0.5).synch_call(E[10], 0.5).replies_to(E[3])
    # T2: front end - Entry E4 (AH5 -> AH6 -> synchCall E5,E6,E7)
    A[8] = Activity(model, 'AH5', Exp.fit_mean(0.00000005)).on(T[2]).bound_to(E[4])
    A[9] = Activity(model, 'AH6', Exp.fit_mean(0.0051774)).on(T[2]).synch_call(E[6], 0.33).synch_call(E[7], 0.33).synch_call(E[8], 0.33).replies_to(E[4])

    # T6: cartdb - Entry E11 (AH15 -> AH16)
    A[10] = Activity(model, 'AH15', Exp.fit_mean(0.0000000005)).on(T[3]).bound_to(E[5])
    A[11] = Activity(model, 'AH16', Exp.fit_mean(0.0040355)).on(T[3]).replies_to(E[5])

    # T3: cart - Entry E5 (AH7 -> AH8 -> synchCall E11)
    A[12] = Activity(model, 'AH7', Exp.fit_mean(0.0000000005)).on(T[4]).bound_to(E[6])
    A[13] = Activity(model, 'AH8', Exp.fit_mean(0.0029469)).on(T[4]).synch_call(E[5], 1).replies_to(E[6])
    # T3: cart - Entry E6 (AH9 -> AH10 -> synchCall E11)
    A[14] = Activity(model, 'AH9', Exp.fit_mean(0.0000000005)).on(T[4]).bound_to(E[7])
    A[15] = Activity(model, 'AH10', Exp.fit_mean(0.012323)).on(T[4]).synch_call(E[5], 1).replies_to(E[7])
    # T3: cart - Entry E7 (AH11 -> AH12 -> synchCall E11)
    A[16] = Activity(model, 'AH11', Exp.fit_mean(0.0000000005)).on(T[4]).bound_to(E[8])
    A[17] = Activity(model, 'AH12', Exp.fit_mean(0.0033488)).on(T[4]).synch_call(E[5], 1).replies_to(E[8])

    # T4: catalog - Entry E8 (AS3 -> AS4 -> synchCall E10)
    A[18] = Activity(model, 'AS3', Exp.fit_mean(0.0000000005)).on(T[5]).bound_to(E[9])
    A[19] = Activity(model, 'AS4', Exp.fit_mean(0.0034925)).on(T[5]).synch_call(E[11], 1).replies_to(E[9])
    # T4: catalog - Entry E9 (AS5 -> AS6 -> synchCall E10)
    A[20] = Activity(model, 'AS5', Exp.fit_mean(0.0000000005)).on(T[5]).bound_to(E[10])
    A[21] = Activity(model, 'AS6', Exp.fit_mean(0.0030162)).on(T[5]).synch_call(E[11], 1).replies_to(E[10])

    # T5: catalogdb - Entry E10 (AH13 -> AH14)
    A[22] = Activity(model, 'AH13', Exp.fit_mean(0.0000000005)).on(T[6]).bound_to(E[11])
    A[23] = Activity(model, 'AH14', Exp.fit_mean(0.0032434)).on(T[6]).replies_to(E[11])

    # Serial precedence for all tasks
    T[0].add_precedence(ActivityPrecedence.serial([A[0], A[1]]))
    T[1].add_precedence(ActivityPrecedence.serial([A[2], A[3]]))
    T[2].add_precedence(ActivityPrecedence.serial([A[4], A[5]]))
    T[2].add_precedence(ActivityPrecedence.serial([A[6], A[7]]))
    T[2].add_precedence(ActivityPrecedence.serial([A[8], A[9]]))
    T[3].add_precedence(ActivityPrecedence.serial([A[10], A[11]]))
    T[4].add_precedence(ActivityPrecedence.serial([A[12], A[13]]))
    T[4].add_precedence(ActivityPrecedence.serial([A[14], A[15]]))
    T[4].add_precedence(ActivityPrecedence.serial([A[16], A[17]]))
    T[5].add_precedence(ActivityPrecedence.serial([A[18], A[19]]))
    T[5].add_precedence(ActivityPrecedence.serial([A[20], A[21]]))
    T[6].add_precedence(ActivityPrecedence.serial([A[22], A[23]]))

    return model


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.STD)
    model = lqn_sockshop()

    # Solve using LN with MVA
    avg_table_ln = LN(model).get_avg_table()
    print('\nLN(MVA) Results:')
    print(avg_table_ln)
