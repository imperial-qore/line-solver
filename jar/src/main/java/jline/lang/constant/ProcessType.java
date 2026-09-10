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
    PHT;

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
                return true;
            default:
                return false;
        }
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