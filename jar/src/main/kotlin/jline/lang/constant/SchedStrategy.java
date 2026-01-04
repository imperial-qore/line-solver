/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

/**
 * Enumeration of scheduling strategies available at queueing stations.
 * <p>
 * These constants specify how jobs are selected for service when multiple
 * jobs are waiting at a station. Each strategy implements different
 * prioritization and fairness policies.
 * <p>
 * Common strategies include:
 * <ul>
 * <li>FCFS - First Come First Serve</li>
 * <li>PS - Processor Sharing</li>
 * <li>LCFS - Last Come First Serve</li>
 * <li>SJF - Shortest Job First</li>
 * <li>INF - Infinite servers (delay station)</li>
 * </ul>
 */
public enum SchedStrategy {
    /** Infinite Server - each job gets its own server immediately (delay station) */
    INF,
    
    /** First Come First Serve - jobs served in order of arrival */
    FCFS,

    /** First Come First Serve Preemptive Resume - higher priority job preempts, preempted job resumes */
    FCFSPR,

    /** First Come First Serve Preemptive Independent - higher priority job preempts, preempted job restarts */
    FCFSPI,

    /** Last Come First Serve - most recent arrival served first (non-preemptive) */
    LCFS,
    
    /** Last Come First Serve Preemptive Resume - arriving job preempts current service */
    LCFSPR,
    
    /** Last Come First Serve Preemptive Independent - arriving job preempts and service restarts from beginning */
    LCFSPI,
    
    /** Service In Random Order - jobs selected randomly from queue */
    SIRO,
    
    /** Shortest Job First - jobs with smallest service requirement served first */
    SJF,
    
    /** Longest Job First - jobs with largest service requirement served first */
    LJF,
    
    /** Processor Sharing - all jobs receive equal share of service capacity */
    PS,
    
    /** Discriminatory Processor Sharing - jobs receive weighted shares based on class */
    DPS,
    
    /** Generalized Processor Sharing - minimum service guarantee with work-conserving */
    GPS,
    
    /** Shortest Expected Processing Time - based on expected remaining service time */
    SEPT,
    
    /** Longest Expected Processing Time - based on expected remaining service time */
    LEPT,

    /** Shortest Remaining Processing Time - preemptive scheduling based on remaining service time */
    SRPT,

    /** Shortest Remaining Processing Time with Priority - SRPT within priority levels */
    SRPTPRIO,

    /** Head Of Line priority (FCFS with priorities) - highest priority class served first */
    HOL, // FCFCPRIO
    
    /** First Come First Served with Priority - alias for HOL */
    FCFSPRIO,
    
    /** Fork node - jobs are split into parallel tasks */
    FORK,
    
    /** External routing - used for open classes at sources */
    EXT,
    
    /** Reference station - used for closed classes */
    REF,
    
    /** Polling - server visits queues in cyclic order */
    POLLING,
    
    /** Processor Sharing with Priorities - PS within priority levels */
    PSPRIO,
    
    /** Discriminatory PS with Priorities - DPS within priority levels */
    DPSPRIO,
    
    /** Generalized PS with Priorities - GPS within priority levels */
    GPSPRIO,
    
    /** Last Come First Served with Priority - LCFS with priority classes */
    LCFSPRIO,
    
    /** Last Come First Served Preemptive Resume with Priority - LCFSPR with priority classes */
    LCFSPRPRIO,

    /** Last Come First Served Preemptive Independent with Priority - LCFSPI with priority classes */
    LCFSPIPRIO,

    /** First Come First Served Preemptive Resume with Priority - FCFSPR with priority classes */
    FCFSPRPRIO,

    /** First Come First Served Preemptive Independent with Priority - FCFSPI with priority classes */
    FCFSPIPRIO,

    /** Earliest Due Date - jobs with earliest deadline served first (non-preemptive) */
    EDD,

    /** Earliest Deadline First - jobs with earliest deadline served first (preemptive) */
    EDF,

    /** Least Progress Scheduling - hybrid PS/FCFS with max concurrent jobs limit */
    LPS,

    /** Preemptive Shortest Job First - priority based on original job size (not remaining) */
    PSJF,

    /** Feedback / Least Attained Service - priority based on attained service (age) */
    FB,

    /** Longest Remaining Processing Time - priority to jobs with longest remaining time */
    LRPT,

    /** Shortest Elapsed Time First - non-preemptive FB/LAS (priority by attained service, no preemption) */
    SETF;

    /**
     * Converts a LINE string representation to a SchedStrategy enum value.
     *
     * @param string the string representation (e.g., "fcfs", "ps", "inf")
     * @return the corresponding SchedStrategy enum value
     * @throws RuntimeException if the string is not recognized
     */
    public static SchedStrategy fromText(String string) {
        switch (string) {
            case "inf":
                return INF;
            case "fcfs":
                return FCFS;
            case "fcfspr":
                return FCFSPR;
            case "fcfspi":
                return FCFSPI;
            case "lcfs":
                return LCFS;
            case "siro":
                return SIRO;
            case "sjf":
                return SJF;
            case "ljf":
                return LJF;
            case "ps":
                return PS;
            case "dps":
                return DPS;
            case "gps":
                return GPS;
            case "sept":
                return SEPT;
            case "lept":
                return LEPT;
            case "srpt":
                return SRPT;
            case "srptprio":
                return SRPTPRIO;
            case "hol":
                return HOL;
            case "fcfsprio":
                return FCFSPRIO;
            case "fork":
                return FORK;
            case "ext":
                return EXT;
            case "ref":
                return REF;
            case "lcfspr":
                return LCFSPR;
            case "lcfspi":
                return LCFSPI;
            case "polling":
                return POLLING;
            case "psprio":
                return PSPRIO;
            case "dpsprio":
                return DPSPRIO;
            case "gpsprio":
                return GPSPRIO;
            case "lcfsprio":
                return LCFSPRIO;
            case "lcfsprprio":
                return LCFSPRPRIO;
            case "lcfspiprio":
                return LCFSPIPRIO;
            case "fcfsprprio":
                return FCFSPRPRIO;
            case "fcfspiprio":
                return FCFSPIPRIO;
            case "edd":
                return EDD;
            case "edf":
                return EDF;
            case "lps":
                return LPS;
            case "psjf":
                return PSJF;
            case "fb":
                return FB;
            case "las":
                return FB;
            case "lrpt":
                return LRPT;
            case "setf":
                return SETF;
            case "set":
                return SETF;
            case "pp":
                // LQNS preemptive priority - maps to FCFS with preemptive resume priority
                return FCFSPRPRIO;
            case "cfs":
                // LQNS completely fair scheduling - maps to Generalized Processor Sharing
                return GPS;
            default:
                throw new RuntimeException("Unable to return a SchedStrategy, check string and try again.");
        }
    }

    /**
     * Converts a SchedStrategy enum value to a feature string for analysis.
     *
     * @param scheduling the scheduling strategy to convert
     * @return the feature string representation
     */
    public static String toFeature(SchedStrategy scheduling) {
        switch (scheduling) {
            case INF:
                return "SchedStrategy_INF";
            case FCFS:
                return "SchedStrategy_FCFS";
            case FCFSPR:
                return "SchedStrategy_FCFSPR";
            case FCFSPI:
                return "SchedStrategy_FCFSPI";
            case LCFS:
                return "SchedStrategy_LCFS";
            case LCFSPR:
                return "SchedStrategy_LCFSPR";
            case LCFSPI:
                return "SchedStrategy_LCFSPI";
            case POLLING:
                return "SchedStrategy_POLLING";
            case SIRO:
                return "SchedStrategy_SIRO";
            case SJF:
                return "SchedStrategy_SJF";
            case LJF:
                return "SchedStrategy_LJF";
            case PS:
                return "SchedStrategy_PS";
            case DPS:
                return "SchedStrategy_DPS";
            case GPS:
                return "SchedStrategy_GPS";
            case SEPT:
                return "SchedStrategy_SEPT";
            case LEPT:
                return "SchedStrategy_LEPT";
            case SRPT:
                return "SchedStrategy_SRPT";
            case SRPTPRIO:
                return "SchedStrategy_SRPTPRIO";
            case HOL:
                return "SchedStrategy_HOL";
            case FCFSPRIO:
                return "SchedStrategy_HOL";
            case EXT:
                return "SchedStrategy_EXT";
            case PSPRIO:
                return "SchedStrategy_PSPRIO";
            case DPSPRIO:
                return "SchedStrategy_DPSPRIO";
            case GPSPRIO:
                return "SchedStrategy_GPSPRIO";
            case LCFSPRIO:
                return "SchedStrategy_LCFSPRIO";
            case LCFSPRPRIO:
                return "SchedStrategy_LCFSPRPRIO";
            case LCFSPIPRIO:
                return "SchedStrategy_LCFSPIPRIO";
            case FCFSPRPRIO:
                return "SchedStrategy_FCFSPRPRIO";
            case FCFSPIPRIO:
                return "SchedStrategy_FCFSPIPRIO";
            case EDD:
                return "SchedStrategy_EDD";
            case EDF:
                return "SchedStrategy_EDF";
            case LPS:
                return "SchedStrategy_LPS";
            case PSJF:
                return "SchedStrategy_PSJF";
            case FB:
                return "SchedStrategy_FB";
            case LRPT:
                return "SchedStrategy_LRPT";
            case SETF:
                return "SchedStrategy_SETF";
            default:
                return "";
        }
    }

    /**
     * Converts a SchedStrategy enum value to its text representation.
     *
     * @param scheduling the scheduling strategy to convert
     * @return the text representation of the strategy
     */
    public static String toText(SchedStrategy scheduling) {
        switch (scheduling) {
            case INF:
                return "inf";
            case FCFS:
                return "fcfs";
            case FCFSPR:
                return "fcfspr";
            case FCFSPI:
                return "fcfspi";
            case LCFS:
                return "lcfs";
            case LCFSPR:
                return "lcfspr";
            case LCFSPI:
                return "lcfspi";
            case POLLING:
                return "polling";
            case SIRO:
                return "siro";
            case PS:
                return "ps";
            case DPS:
                return "dps";
            case GPS:
                return "gps";
            case SEPT:
                return "sept";
            case LEPT:
                return "lept";
            case SRPT:
                return "srpt";
            case SRPTPRIO:
                return "srptprio";
            case HOL:
                return "hol";
            case FCFSPRIO:
                return "hol";
            case FORK:
                return "fork";
            case EXT:
                return "ext";
            case REF:
                return "ref";
            case PSPRIO:
                return "psprio";
            case DPSPRIO:
                return "dpsprio";
            case GPSPRIO:
                return "gpsprio";
            case LCFSPRIO:
                return "lcfsprio";
            case LCFSPRPRIO:
                return "lcfsprprio";
            case LCFSPIPRIO:
                return "lcfspiprio";
            case FCFSPRPRIO:
                return "fcfsprprio";
            case FCFSPIPRIO:
                return "fcfspiprio";
            case EDD:
                return "edd";
            case EDF:
                return "edf";
            case LPS:
                return "lps";
            case PSJF:
                return "psjf";
            case FB:
                return "fb";
            case LRPT:
                return "lrpt";
            case SETF:
                return "setf";
            default:
                return "";
        }
    }

}
