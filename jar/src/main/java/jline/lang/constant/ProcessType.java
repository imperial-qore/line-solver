/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

import jline.lang.processes.Distribution;

import java.io.Serializable;

/**
 * Constants for specifying a point process type
 */
public enum ProcessType implements Serializable {
    EXP,
    ERLANG,
    DISABLED,
    IMMEDIATE,
    HYPEREXP,
    APH,
    COXIAN,
    PH,
    MAP,
    DMAP,
    UNIFORM,
    DET,
    GAMMA,
    PARETO,
    WEIBULL,
    LOGNORMAL,
    MMPP2,
    BMAP,
    MMAP,
    ME,
    RAP,
    REPLAYER,
    TRACE,
    COX2,
    BINOMIAL,
    POISSON,
    GEOMETRIC,
    DUNIFORM,
    BERNOULLI,
    PRIOR,
    DISCRETESAMPLER,
    ZIPF,
    NHPP,
    MAPT,
    PHT,
    /**
     * The MARKED families. A mark is a label carried by an event, and at a Source it
     * selects the class of the arriving job (Source.setMarkedArrival, sn.markidx).
     * MPH is the renewal special case of MMAP, obtained from a PH with K marked exits
     * by D0 = S and D1k = s_k*alpha; MMAPT is MMAP with a piecewise-constant matrix
     * schedule, and MPHT is the same lowering applied segment by segment. MPH keeps an
     * id of its own rather than aliasing MMAP, so a solver that cannot honour a marked
     * renewal process refuses it explicitly instead of inheriting MMAP's support.
     */
    MPH,
    MMAPT,
    MPHT,
    /**
     * BMMAPT crosses the BATCH axis with the two above: a block is indexed by segment,
     * mark and batch size, and one epoch releases a batch of jobs that all carry the same
     * mark. It reduces to MMAPT when every batch size is 1 and to BMAP when the schedule is
     * flat, and like MPH it keeps an id of its own so a solver that cannot release batches
     * refuses it by name instead of inheriting MMAPT's support.
     */
    BMMAPT;

    public static ProcessType fromDistribution(Distribution d) {
        return fromText(d.getName());
    }

    /**
     * True when sn.proc carries an exact matrix representation of a process of
     * this type: a genuine (D0,D1) pair, or its matrix-exponential analogue for
     * ME/RAP.
     *
     * <p>The distinction is about sn.proc, NOT about what
     * {@link jline.lang.processes.Distribution#getProcess()} returns.
     * getProcess hands back raw distribution PARAMETERS for several
     * non-Markovian families -- Gamma, Weibull, Lognormal, Pareto and Uniform
     * return two scalars (Pareto {alpha,k}, Uniform {min,max}) -- and the
     * network refresh replaces those with map_erlang(mean, n) before storing
     * them in sn.proc, where n = ceil(1/SCV) capped at 100 (n = 20 when
     * SCV &lt; CoarseTol). That fit matches the mean, and matches the SCV only
     * when SCV &lt;= 1: Pareto with SCV 64 gives n = 1, a single exponential of
     * SCV 1. So for these types sn.proc is an approximation, not the law that
     * was requested, and nothing on the cell says so -- the only signal is
     * sn.procid.</p>
     *
     * <p>Solvers that read sn.proc as if it were the exact law must gate on this
     * predicate. It is the procid-level counterpart of
     * Distribution.isMarkovian(), i.e. of the Markovian class hierarchy, so the
     * two must stay in step.</p>
     *
     * @param t process type
     * @return true if sn.proc holds the exact law for this type
     */
    public static boolean isMarkovian(ProcessType t) {
        if (t == null) {
            return false;
        }
        switch (t) {
            case EXP:
            case ERLANG:
            case HYPEREXP:
            case PH:
            case APH:
            case MAP:
            case COXIAN:
            case COX2:
            case MMPP2:
            case ME:
            case RAP:
            case DMAP:
            case BMAP:
            case MMAP:
            case MPH:
            case MMAPT:
            case MPHT:
            case BMMAPT:
                return true;
            default:
                return false;
        }
    }

    /**
     * True when the type carries PER-MARK arrival blocks, i.e. an event of this process
     * is labelled and the label is meaningful to the model. At a Source the label selects
     * the class of the arriving job (Source.setMarkedArrival, sn.markidx).
     *
     * @param t process type
     * @return true for MMAP, MPH, MMAPT, MPHT and BMMAPT
     */
    public static boolean isMarked(ProcessType t) {
        return isMarkedStationary(t) || isMarkedSchedule(t);
    }

    /**
     * True when sn.proc holds the STATIONARY marked cell, the M3A layout
     * {D0, D1_agg, D11, ..., D1K}.
     *
     * <p>MPH is the renewal special case of MMAP and lowers to exactly that cell
     * (D0 = S, D1k = s_k*alpha), so every consumer that reads the M3A layout serves
     * both and must gate on this predicate rather than on equality with MMAP.
     *
     * @param t process type
     * @return true for MMAP and MPH
     */
    public static boolean isMarkedStationary(ProcessType t) {
        return t == ProcessType.MMAP || t == ProcessType.MPH;
    }

    /**
     * True when sn.proc holds the MARKED SCHEDULE slot,
     * {breakpoints, D0 segments, D1_agg segments, cyclic, per-mark segments}.
     *
     * <p>MPHT is stored lowered to MMAPT form segment by segment, so one walk serves
     * both, exactly as one MAPt walk serves MAPt and PHt.
     *
     * <p>BMMAPT IS INCLUDED, and its slot is that one with the batch blocks appended, so a
     * consumer gated on this predicate reads a BMMAPT as the MMAPT it aggregates down to.
     * That is right for anything time-blind or batch-blind and WRONG for anything that
     * releases jobs: an arrival or service sampler must branch on {@link #isBatch} as well,
     * or it silently delivers one job per epoch.
     *
     * @param t process type
     * @return true for MMAPT, MPHT and BMMAPT
     */
    public static boolean isMarkedSchedule(ProcessType t) {
        return t == ProcessType.MMAPT || t == ProcessType.MPHT || t == ProcessType.BMMAPT;
    }

    /**
     * True when an EVENT of this process releases (or, as a service process, completes) a
     * BATCH of jobs whose size the process itself carries in its blocks.
     *
     * <p>This is the batch twin of {@link #isMarkedStationary} and {@link #isMarkedSchedule},
     * and the same rule applies: a procid test that should serve every batch family is a
     * MEMBERSHIP test, never equality with BMAP.
     *
     * <p>It is disjoint from {@code sn.arrivalbatch}, which is a SEPARATE batch-size law
     * bolted onto a renewal stream by {@code Source.setArrivalBatch}. A process that is
     * isBatch already carries its own sizes, so the two are mutually exclusive by
     * construction.
     *
     * @param t process type
     * @return true for BMAP and BMMAPT
     */
    public static boolean isBatch(ProcessType t) {
        return t == ProcessType.BMAP || t == ProcessType.BMMAPT;
    }

    public static ProcessType fromText(String name) {
        switch (name) {
            case "Exp":
                return ProcessType.EXP;
            case "Erlang":
                return ProcessType.ERLANG;
            case "HyperExp":
                return ProcessType.HYPEREXP;
            case "PH":
                return ProcessType.PH;
            case "APH":
                return ProcessType.APH;
            case "MAP":
                return ProcessType.MAP;
            case "DMAP":
                return ProcessType.DMAP;
            case "Uniform":
                return ProcessType.UNIFORM;
            case "Det":
                return ProcessType.DET;
            case "Coxian":
                return ProcessType.COXIAN;
            case "Gamma":
                return ProcessType.GAMMA;
            case "Pareto":
                return ProcessType.PARETO;
            case "MMPP2":
                return ProcessType.MMPP2;
            case "BMAP":
                return ProcessType.BMAP;
            case "MMAP":
            case "MarkedMAP":
                return ProcessType.MMAP;
            case "MPH":
            case "MarkedPH":
                return ProcessType.MPH;
            case "BMMAPt":
                return ProcessType.BMMAPT;
            case "MMAPt":
                return ProcessType.MMAPT;
            case "MPHt":
                return ProcessType.MPHT;
            case "ME":
                return ProcessType.ME;
            case "RAP":
                return ProcessType.RAP;
            case "Replayer":
            case "Trace":
                return ProcessType.REPLAYER;
            case "Immediate":
                return ProcessType.IMMEDIATE;
            case "Disabled":
                return ProcessType.DISABLED;
            case "Cox2":
                return ProcessType.COX2;
            case "Weibull":
                return ProcessType.WEIBULL;
            case "Lognormal":
                return ProcessType.LOGNORMAL;
            case "Poisson":
                return ProcessType.POISSON;
            case "Binomial":
                return ProcessType.BINOMIAL;
            case "Geometric":
                return ProcessType.GEOMETRIC;
            case "DiscreteUniform":
                return ProcessType.DUNIFORM;
            case "Bernoulli":
                return ProcessType.BERNOULLI;
            case "Prior":
                return ProcessType.PRIOR;
            case "DiscreteSampler":
                return ProcessType.DISCRETESAMPLER;
            case "Zipf":
                return ProcessType.ZIPF;
            case "NHPP":
                return ProcessType.NHPP;
            case "MAPt":
                return ProcessType.MAPT;
            case "PHt":
                return ProcessType.PHT;
            default:
                throw new IllegalArgumentException("Unknown ProcessType: " + name);
        }
    }

    public static String toText(ProcessType type) {
        switch (type) {
            case EXP:
                return "Exp";
            case ERLANG:
                return "Erlang";
            case HYPEREXP:
                return "HyperExp";
            case PH:
                return "PH";
            case APH:
                return "APH";
            case MAP:
                return "MAP";
            case DMAP:
                return "DMAP";
            case UNIFORM:
                return "Uniform";
            case DET:
                return "Det";
            case COXIAN:
                return "Coxian";
            case GAMMA:
                return "Gamma";
            case PARETO:
                return "Pareto";
            case MMPP2:
                return "MMPP2";
            case BMAP:
                return "BMAP";
            case MMAP:
                return "MMAP";
            case MPH:
                return "MPH";
            case MMAPT:
                return "MMAPt";
            case BMMAPT:
                return "BMMAPt";
            case MPHT:
                return "MPHt";
            case ME:
                return "ME";
            case RAP:
                return "RAP";
            case REPLAYER:
            case TRACE:
                return "Replayer";
            case IMMEDIATE:
                return "Immediate";
            case DISABLED:
                return "Disabled";
            case COX2:
                return "Cox2";
            case WEIBULL:
                return "Weibull";
            case LOGNORMAL:
                return "Lognormal";
            case POISSON:
                return "Poisson";
            case BINOMIAL:
                return "Binomial";
            case GEOMETRIC:
                return "Geometric";
            case DUNIFORM:
                return "DiscreteUniform";
            case BERNOULLI:
                return "Bernoulli";
            case PRIOR:
                return "Prior";
            case DISCRETESAMPLER:
                return "DiscreteSampler";
            case ZIPF:
                return "Zipf";
            case NHPP:
                return "NHPP";
            case MAPT:
                return "MAPt";
            case PHT:
                return "PHt";
            default:
                return type.name();
        }
    }
}