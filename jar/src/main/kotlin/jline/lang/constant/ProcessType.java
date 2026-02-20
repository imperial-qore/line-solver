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
    ZIPF;

    public static ProcessType fromDistribution(Distribution d) {
        return fromText(d.getName());
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
            default:
                return type.name();
        }
    }
}