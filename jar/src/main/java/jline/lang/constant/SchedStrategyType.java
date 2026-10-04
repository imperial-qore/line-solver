/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

/**
 * Constants for specifying a scheduling strategy type at stations.
 *
 * <p>Declaration order is MATLAB's numbering in
 * {@code matlab/src/lang/constant/SchedStrategyType.m} (NP = 0), so an ordinal
 * means the same thing in all three codebases. It used to start at PR, and the
 * two Prio members below did not exist here at all.</p>
 */
public enum SchedStrategyType {
    NP,      // non-preemptive
    PR,      // preemptive resume
    PNR,     // preemptive non-resume
    NPPrio,  // non-preemptive priority
    PRPrio,  // preemptive resume priority
    PNRPrio; // preemptive non-resume priority

    /**
     * Classifies a scheduling discipline by preemption behaviour and priority.
     *
     * <p>This IS the classification {@code Queue} stores in {@code schedPolicy}:
     * the constructor calls this and keeps no table of its own, as MATLAB's
     * {@code Queue.m} and native python's {@code Queue} do. The three were
     * written out separately until 2026-09-17 and disagreed on eleven
     * disciplines: this table answered PNR for the preemptive-identical family,
     * PR for PSJF/FB/LRPT and NPPrio for LCFSPRIO where MATLAB answered PR, NP
     * and NP. Nothing in any codebase READS schedPolicy, so the disagreement
     * was never a wrong number, only a wrong label, which is why it survived.</p>
     *
     * <p>The preemption axis follows each discipline's own definition in
     * {@link SchedStrategy}. PI is "preemptive independent", i.e. the preempted
     * job RESTARTS, so it is non-resume and not resume; SRPT, PSJF, FB, LRPT,
     * FSP and EDF preempt, against SJF, LJF, SEPT, LEPT, EDD and SETF, which
     * rank the same jobs without preempting.</p>
     *
     * @param strategy the scheduling discipline
     * @return the scheduling strategy type
     */
    public static SchedStrategyType getTypeId(SchedStrategy strategy) {
        switch (strategy) {
            case INF:
            case FCFS:
            case LCFS:
            case SIRO:
            case SJF:
            case LJF:
            case SEPT:
            case LEPT:
            case EDD:
            case SETF:
            case PAS:
            case OI:
            case POLLING:
                return NP;
            case PS:
            case DPS:
            case GPS:
            case LPS:
            case LCFSPR:
            case FCFSPR:
            case EDF:
            case SRPT:
            case PSJF:
            case FB:
            case LRPT:
            case FSP:
                return PR;
            case FCFSPI:
            case LCFSPI:
                return PNR;
            case HOL:
            case FCFSPRIO:
            case LCFSPRIO:
                return NPPrio;
            case PSPRIO:
            case DPSPRIO:
            case GPSPRIO:
            case LCFSPRPRIO:
            case FCFSPRPRIO:
            case SRPTPRIO:
                return PRPrio;
            case FCFSPIPRIO:
            case LCFSPIPRIO:
                return PNRPrio;
            default:
                throw new RuntimeException("Unrecognized scheduling strategy type: " + strategy);
        }
    }
}
